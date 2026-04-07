# -*- coding: utf-8 -*-
"""
Intra-tumor Heterogeneity (ITH) Extraction and ML Pipeline
Includes:
1. Unsupervised ITH Feature Extraction (Pre, Post, and Delta)
2. Spatial Cluster Visualization
3. Machine Learning Evaluation Pipeline (SMOTE + Classifiers)
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import SimpleITK as sitk
import matplotlib.pyplot as plt
import joblib

from sklearn.cluster import KMeans
from scipy.ndimage import label as cc_label, generate_binary_structure
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.metrics import (
    roc_auc_score, roc_curve, accuracy_score, f1_score, precision_score, 
    recall_score, confusion_matrix, balanced_accuracy_score
)
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier
from imblearn.pipeline import Pipeline as ImbPipeline
from imblearn.over_sampling import BorderlineSMOTE

warnings.filterwarnings("ignore")

# ==========================================================
# 0. CONFIGURATION (Please update paths before running)
# ==========================================================
class Config:
    # --- Paths ---
    PRE_IMG_DIR   = "./data/First-image"
    PRE_MASK_DIR  = "./data/First-mask"
    POST_IMG_DIR  = "./data/Second-image"
    POST_MASK_DIR = "./data/Second-mask"
    CLINICAL_CSV  = "./data/clinical_data.csv"
    OUTPUT_DIR    = "./output"
    
    # --- ITH Parameters ---
    N_CLUSTERS = 3
    RANDOM_STATE = 42
    EPS = 1e-8
    
    # --- ML Parameters ---
    TEST_SIZE = 0.30
    USE_SMOTE = True
    SMOTE_K_NEIGHBORS = 5
    
    # --- Clinical Column Matchers ---
    ID_COL_CANDIDATES = ["ID", "id", "Case", "case", "case_id", "CaseID", "PatientID"]
    LABEL_COL_CANDIDATES = ["group", "label", "Label", "y", "response"]

os.makedirs(Config.OUTPUT_DIR, exist_ok=True)

# ==========================================================
# 1. UTILITIES
# ==========================================================
def list_nii_files(folder):
    if not os.path.isdir(folder): return []
    return sorted([f for f in os.listdir(folder) if f.lower().endswith(('.nii', '.nii.gz'))])

def case_id_from_filename(fname: str) -> str:
    base = fname
    if base.lower().endswith(".nii.gz"): return base[:-7]
    elif base.lower().endswith(".nii"): return base[:-4]
    return base

def extract_int_id(x):
    if pd.isna(x): return np.nan
    s = str(x).strip()
    if not s: return np.nan
    try:
        f = float(s)
        if np.isfinite(f): return int(round(f))
    except Exception: pass
    m = re.findall(r"\d+", s)
    if not m: return np.nan
    return int(m[-1])

def percent_delta(pre_val, post_val, eps=1e-8):
    pre = float(pre_val)
    post = float(post_val)
    denom = pre if abs(pre) > eps else (eps if pre >= 0 else -eps)
    return (post - pre) / denom * 100.0

def load_clinical_data(csv_path):
    encodings = ["utf-8-sig", "utf-8", "gbk", "latin1"]
    for enc in encodings:
        try:
            df = pd.read_csv(csv_path, encoding=enc)
            # Find ID and Label columns
            id_col = next((c for c in df.columns if any(cand.lower() in c.lower() for cand in Config.ID_COL_CANDIDATES)), None)
            lbl_col = next((c for c in df.columns if any(cand.lower() in c.lower() for cand in Config.LABEL_COL_CANDIDATES)), None)
            
            if not id_col or not lbl_col:
                raise ValueError(f"Could not find ID or Label column in {csv_path}")
                
            df["_merge_id_int_"] = df[id_col].astype(str).apply(extract_int_id)
            return df, id_col, lbl_col
        except UnicodeDecodeError:
            continue
    raise RuntimeError("Failed to read clinical CSV with common encodings.")

# ==========================================================
# 2. ITH EXTRACTION (PURE UNSUPERVISED)
# ==========================================================
def calculate_connectivity_metric(cluster_label_map, spacing_xyz, n_clusters, is_3d=True):
    sx, sy, sz = spacing_xyz
    unit = float(sx) * float(sy) * float(sz) if is_3d else float(sx) * float(sy)
    struct = generate_binary_structure(3, 3) if is_3d else generate_binary_structure(2, 2)
    
    fg = (cluster_label_map > 0)
    total_size = np.sum(fg) * unit
    if total_size <= 0: return 0.0

    sum_term = 0.0
    for l in range(1, n_clusters + 1):
        cm = (cluster_label_map == l).astype(np.uint8)
        if cm.sum() == 0: continue
        labeled_arr, num_cc = cc_label(cm, structure=struct)
        if num_cc <= 0: continue
        
        sizes = [(np.sum(labeled_arr == cc_id) * unit) for cc_id in range(1, num_cc + 1)]
        max_size = float(max(sizes)) if sizes else 0.0
        sum_term += (max_size / float(num_cc))

    ith = 1.0 - (1.0 / total_size) * sum_term
    return float(max(min(ith, 1.0), 0.0))

def calc_ith_unsupervised(image_sitk, mask_sitk, n_clusters=3):
    img = sitk.GetArrayFromImage(image_sitk)
    msk = sitk.GetArrayFromImage(mask_sitk)
    spacing = image_sitk.GetSpacing()

    roi = (msk > 0)
    voxels = img[roi].reshape(-1, 1).astype(np.float64)
    if voxels.shape[0] < max(10, n_clusters * 2):
        return 0.0, 0.0

    # PURE UNSUPERVISED KMEANS
    km = KMeans(n_clusters=n_clusters, random_state=Config.RANDOM_STATE, n_init=5)
    labels = km.fit_predict(voxels)

    cluster_map = np.zeros_like(msk, dtype=np.int16)
    cluster_map[roi] = (labels + 1).astype(np.int16)

    # 2D TLCRH (Max Slice)
    slice_areas = np.sum(roi, axis=(1, 2))
    tlcrh = 0.0
    if slice_areas.max() > 0:
        z = int(np.argmax(slice_areas))
        tlcrh = calculate_connectivity_metric(cluster_map[z, :, :], spacing, n_clusters, is_3d=False)

    # 3D TSH
    tsh = calculate_connectivity_metric(cluster_map, spacing, n_clusters, is_3d=True)
    return tsh, tlcrh, cluster_map

def run_extraction():
    print("--- Starting ITH Feature Extraction ---")
    files = list_nii_files(Config.PRE_IMG_DIR)
    clin_df, _, _ = load_clinical_data(Config.CLINICAL_CSV)
    
    results = []
    for f in files:
        cid = case_id_from_filename(f)
        pre_img = sitk.ReadImage(os.path.join(Config.PRE_IMG_DIR, f))
        pre_msk = sitk.ReadImage(os.path.join(Config.PRE_MASK_DIR, f))
        post_img = sitk.ReadImage(os.path.join(Config.POST_IMG_DIR, f))
        post_msk = sitk.ReadImage(os.path.join(Config.POST_MASK_DIR, f))
        
        pre_tsh, pre_tlcrh, _ = calc_ith_unsupervised(pre_img, pre_msk, Config.N_CLUSTERS)
        post_tsh, post_tlcrh, _ = calc_ith_unsupervised(post_img, post_msk, Config.N_CLUSTERS)
        
        results.append({
            "ID": cid,
            "_merge_id_int_": extract_int_id(cid),
            "pre_ITH_TSH": pre_tsh,
            "post_ITH_TSH": post_tsh,
            "pctDelta_ITH_TSH": percent_delta(pre_tsh, post_tsh),
            "pre_ITH_TLCRH": pre_tlcrh,
            "post_ITH_TLCRH": post_tlcrh,
            "pctDelta_ITH_TLCRH": percent_delta(pre_tlcrh, post_tlcrh),
        })
        print(f"Extracted features for {cid}")

    df_feats = pd.DataFrame(results)
    df_merged = df_feats.merge(clin_df, on="_merge_id_int_", how="left")
    
    out_path = os.path.join(Config.OUTPUT_DIR, "Extracted_ITH_Features.csv")
    df_merged.to_csv(out_path, index=False)
    print(f"Extraction complete. Saved to {out_path}\n")
    return out_path

# ==========================================================
# 3. VISUALIZATION
# ==========================================================
def create_overlay(img_gray, mask_bin, alpha=0.4, color=(1, 0, 0)):
    img_norm = img_gray - np.min(img_gray)
    max_val = np.max(img_norm)
    if max_val > 0: img_norm = img_norm / max_val
    overlay = np.stack([img_norm, img_norm, img_norm], axis=-1)
    mask_indices = mask_bin > 0
    for c in range(3):
        overlay[mask_indices, c] = overlay[mask_indices, c] * (1 - alpha) + color[c] * alpha
    return overlay

def run_visualization():
    print("--- Starting ITH Visualization ---")
    vis_dir = os.path.join(Config.OUTPUT_DIR, "Visualizations")
    os.makedirs(vis_dir, exist_ok=True)
    files = list_nii_files(Config.PRE_IMG_DIR)
    
    for f in files[:5]: # Process first 5 for demonstration
        cid = case_id_from_filename(f)
        try:
            pre_img_sitk = sitk.ReadImage(os.path.join(Config.PRE_IMG_DIR, f))
            pre_msk_sitk = sitk.ReadImage(os.path.join(Config.PRE_MASK_DIR, f))
            post_img_sitk = sitk.ReadImage(os.path.join(Config.POST_IMG_DIR, f))
            post_msk_sitk = sitk.ReadImage(os.path.join(Config.POST_MASK_DIR, f))
            
            _, _, pre_cmap = calc_ith_unsupervised(pre_img_sitk, pre_msk_sitk, Config.N_CLUSTERS)
            _, _, post_cmap = calc_ith_unsupervised(post_img_sitk, post_msk_sitk, Config.N_CLUSTERS)
            
            # Crop to ROI Z-slice
            def extract_roi_slice(img_sitk, msk_sitk, cmap):
                img = sitk.GetArrayFromImage(img_sitk)
                msk = sitk.GetArrayFromImage(msk_sitk)
                z = int(np.argmax(np.sum(msk > 0, axis=(1,2))))
                return img[z], msk[z], cmap[z]

            pre_img, pre_msk, pre_c = extract_roi_slice(pre_img_sitk, pre_msk_sitk, pre_cmap)
            post_img, post_msk, post_c = extract_roi_slice(post_img_sitk, post_msk_sitk, post_cmap)
            
            fig, axes = plt.subplots(2, 3, figsize=(15, 9))
            
            axes[0,0].imshow(pre_img, cmap='gray'); axes[0,0].set_title(f"PRE Original ({cid})")
            axes[0,1].imshow(create_overlay(pre_img, pre_msk)); axes[0,1].set_title("PRE Mask Overlay")
            axes[0,2].imshow(pre_c, cmap='viridis'); axes[0,2].set_title("PRE ITH Cluster Map")
            
            axes[1,0].imshow(post_img, cmap='gray'); axes[1,0].set_title("POST Original")
            axes[1,1].imshow(create_overlay(post_img, post_msk)); axes[1,1].set_title("POST Mask Overlay")
            axes[1,2].imshow(post_c, cmap='viridis'); axes[1,2].set_title("POST ITH Cluster Map")
            
            for ax in axes.ravel(): ax.axis('off')
            plt.tight_layout()
            plt.savefig(os.path.join(vis_dir, f"{cid}_ITH_Vis.pdf"))
            plt.close()
            print(f"Saved visualization for {cid}")
        except Exception as e:
            print(f"Failed to visualize {cid}: {e}")
    print("Visualization complete.\n")

# ==========================================================
# 4. MACHINE LEARNING PIPELINE
# ==========================================================
def youden_threshold(y_true, y_prob):
    fpr, tpr, thr = roc_curve(y_true, y_prob)
    idx = int(np.argmax(tpr - fpr))
    return float(thr[idx])

def build_model(name):
    steps = [
        ("imputer", SimpleImputer(strategy="median")),
        ("poly", PolynomialFeatures(degree=2, include_bias=False)),
        ("scaler", StandardScaler())
    ]
    if Config.USE_SMOTE:
        steps.append(("smote", BorderlineSMOTE(random_state=Config.RANDOM_STATE, k_neighbors=Config.SMOTE_K_NEIGHBORS)))
        
    models = {
        "LogReg": LogisticRegression(C=1.0, solver="liblinear", class_weight="balanced", random_state=Config.RANDOM_STATE),
        "SVM-RBF": SVC(C=1.0, kernel="rbf", probability=True, class_weight="balanced", random_state=Config.RANDOM_STATE),
        "RF": RandomForestClassifier(n_estimators=500, max_depth=5, class_weight="balanced_subsample", random_state=Config.RANDOM_STATE),
    }
    steps.append(("clf", models[name]))
    return ImbPipeline(steps) if Config.USE_SMOTE else Pipeline(steps)

def run_ml_pipeline(csv_path):
    print("--- Starting Machine Learning Pipeline ---")
    df = pd.read_csv(csv_path)
    
    # Locate label
    lbl_col = next((c for c in df.columns if any(cand.lower() in c.lower() for cand in Config.LABEL_COL_CANDIDATES)), None)
    if not lbl_col: raise ValueError("Label column not found in extracted data.")
    
    df = df.dropna(subset=[lbl_col]).reset_index(drop=True)
    y = df[lbl_col].astype(int)
    
    tasks = {
        "PRE": ["pre_ITH_TSH", "pre_ITH_TLCRH"],
        "POST": ["post_ITH_TSH", "post_ITH_TLCRH"],
        "DELTA": ["pctDelta_ITH_TSH", "pctDelta_ITH_TLCRH"]
    }
    
    ml_out = os.path.join(Config.OUTPUT_DIR, "ML_Results")
    os.makedirs(ml_out, exist_ok=True)
    
    overall_results = []
    
    for task_name, features in tasks.items():
        X = df[features].apply(pd.to_numeric, errors='coerce')
        valid = ~X.isna().all(axis=1)
        X, y_valid = X[valid], y[valid]
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_valid, test_size=Config.TEST_SIZE, stratify=y_valid, random_state=Config.RANDOM_STATE
        )
        
        plt.figure(figsize=(8, 6))
        for model_name in ["LogReg", "SVM-RBF", "RF"]:
            pipe = build_model(model_name)
            pipe.fit(X_train, y_train)
            
            # Predict
            train_prob = pipe.predict_proba(X_train)[:, 1]
            test_prob = pipe.predict_proba(X_test)[:, 1]
            
            thr = youden_threshold(y_train, train_prob)
            y_pred = (test_prob >= thr).astype(int)
            
            auc = roc_auc_score(y_test, test_prob)
            overall_results.append({
                "Task": task_name, "Model": model_name,
                "AUC": auc, "Accuracy": accuracy_score(y_test, y_pred),
                "F1": f1_score(y_test, y_pred)
            })
            
            fpr, tpr, _ = roc_curve(y_test, test_prob)
            plt.plot(fpr, tpr, label=f'{model_name} (AUC={auc:.3f})')
        
        plt.plot([0, 1], [0, 1], 'k--')
        plt.title(f"ROC Curves - Task: {task_name}")
        plt.legend()
        plt.savefig(os.path.join(ml_out, f"ROC_{task_name}.png"))
        plt.close()
        print(f"Evaluated ML models for Task: {task_name}")

    pd.DataFrame(overall_results).to_csv(os.path.join(ml_out, "ML_Summary.csv"), index=False)
    print("Machine learning pipeline complete. Summary saved.\n")

# ==========================================================
# 5. MAIN EXECUTION
# ==========================================================
if __name__ == "__main__":
    # You can comment/uncomment these depending on what you want to run
    extracted_csv = run_extraction()
    run_visualization()
    run_ml_pipeline(extracted_csv)