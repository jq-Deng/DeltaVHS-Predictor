#!/usr/bin/env python3


import os
import warnings
import numpy as np
import pandas as pd
from scipy.ndimage import label as cc_label, generate_binary_structure
from scipy.ndimage import sobel, laplace
from sklearn.cluster import KMeans

warnings.filterwarnings("ignore")


class Config:


    TARGET_SPACING = None
    IMAGE_INTERPOLATOR = "bspline"
    MASK_INTERPOLATOR = "nearest"


    INTENSITY_CLIP_LOWER_PERCENTILE = 1
    INTENSITY_CLIP_UPPER_PERCENTILE = 99
    BIN_WIDTH_HU = 25


    WINDOW_SIZE = 3
    N_STANDARD_FEATURES = 93
    N_CUSTOM_FEATURES = 11
    N_TOTAL_FEATURES = 104


    N_CLUSTERS = 3
    KMEANS_INIT = "k-means++"
    KMEANS_N_INIT = 5
    KMEANS_RANDOM_STATE = 42


    CONNECTIVITY = 26


    VHS_CLIP_LOWER = 1e-4
    VHS_CLIP_UPPER = 1.0


def load_nifti(image_path, mask_path):

    import SimpleITK as sitk
    image = sitk.ReadImage(image_path)
    mask = sitk.ReadImage(mask_path)
    return image, mask


def get_native_min_spacing(image_path):


    import SimpleITK as sitk
    img = sitk.ReadImage(image_path)
    return float(min(img.GetSpacing()))


def resample_to_isotropic(image_sitk, mask_sitk, target_spacing=None):


    import SimpleITK as sitk

    original_spacing = image_sitk.GetSpacing()

    if target_spacing is None:
        target_spacing = min(original_spacing)
    if not np.isfinite(target_spacing) or target_spacing <= 0:
        raise ValueError("target_spacing must be a positive finite number")

    target_spacing_3d = (target_spacing, target_spacing, target_spacing)

    original_size = image_sitk.GetSize()
    scaling_factor = [original_spacing[i] / target_spacing_3d[i] for i in range(3)]
    target_size = [int(round(original_size[i] * scaling_factor[i])) for i in range(3)]


    resampler = sitk.ResampleImageFilter()
    resampler.SetOutputSpacing(target_spacing_3d)
    resampler.SetSize(target_size)
    resampler.SetOutputDirection(image_sitk.GetDirection())
    resampler.SetOutputOrigin(image_sitk.GetOrigin())
    resampler.SetTransform(sitk.Transform())
    resampler.SetDefaultPixelValue(image_sitk.GetPixelIDValue())

    resampler.SetInterpolator(sitk.sitkBSpline)
    image_resampled = resampler.Execute(image_sitk)


    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    mask_resampled = resampler.Execute(mask_sitk)

    return image_resampled, mask_resampled


def retain_largest_connected_component(mask_arr):


    structure = generate_binary_structure(3, 3)
    labeled, num_features = cc_label(mask_arr > 0, structure=structure)

    if num_features == 0:
        return mask_arr


    component_sizes = np.bincount(labeled.ravel())
    component_sizes[0] = 0
    largest_component = np.argmax(component_sizes)

    mask_largest = (labeled == largest_component).astype(np.uint8)
    return mask_largest


def clip_intensities(img_arr, mask_arr, lower_percentile=1, upper_percentile=99):


    voi_intensities = img_arr[mask_arr > 0]
    if voi_intensities.size == 0:
        raise ValueError("VAT mask contains no valid voxels after preprocessing")
    lower = np.percentile(voi_intensities, lower_percentile)
    upper = np.percentile(voi_intensities, upper_percentile)
    img_clipped = np.clip(img_arr, lower, upper)
    return img_clipped


def discretize_to_bins(img_arr, mask_arr, bin_width_hu=25):


    voi_intensities = img_arr[mask_arr > 0]
    if voi_intensities.size == 0:
        raise ValueError("VAT mask contains no valid voxels for discretization")
    if bin_width_hu <= 0:
        raise ValueError("bin_width_hu must be positive")
    min_voi = np.min(voi_intensities)


    img_discrete = np.floor((img_arr - min_voi) / bin_width_hu).astype(np.int32) + 1
    img_discrete[mask_arr == 0] = -1

    n_levels = int(np.max(img_discrete[mask_arr > 0])) + 1
    return img_discrete, n_levels


def compute_first_order_features(values, discrete_values=None, voxel_volume=1.0):


    values = np.asarray(values, dtype=np.float64)
    n = len(values)
    if n == 0:
        return [0.0] * 18

    mean = np.mean(values)
    variance = np.var(values, ddof=0)
    std = np.sqrt(variance) if variance > 0 else 0.0


    if discrete_values is None:
        raise ValueError("discrete_values (scan-level 25-HU bins) are required")
    discrete_values = np.asarray(discrete_values, dtype=np.int64)
    if discrete_values.shape != values.shape:
        raise ValueError("Raw and discretized first-order input shapes must match")
    _, hist = np.unique(discrete_values, return_counts=True)
    prob = hist.astype(np.float64) / n
    prob_nonzero = prob[prob > 0]

    entropy = -np.sum(prob_nonzero * np.log2(prob_nonzero))
    uniformity = np.sum(prob ** 2)


    p10 = np.percentile(values, 10)
    p90 = np.percentile(values, 90)
    p75 = np.percentile(values, 75)
    p25 = np.percentile(values, 25)

    energy = np.sum(values ** 2)
    total_energy = energy * float(voxel_volume)

    mad = np.mean(np.abs(values - mean))
    rmad = np.mean(np.abs(values[(values >= p10) & (values <= p90)] -
                            np.mean(values[(values >= p10) & (values <= p90)]))) if np.sum((values >= p10) & (values <= p90)) > 0 else 0.0
    rms = np.sqrt(np.mean(values ** 2))


    if std > 0:
        skewness = np.mean(((values - mean) / std) ** 3)
        kurtosis = np.mean(((values - mean) / std) ** 4) - 3.0
    else:
        skewness = 0.0
        kurtosis = 0.0

    features = [
        energy,
        total_energy,
        entropy,
        np.min(values),
        p10,
        p90,
        np.max(values),
        mean,
        np.median(values),
        p75 - p25,
        np.max(values) - np.min(values),
        mad,
        rmad,
        rms,
        skewness,
        kurtosis,
        variance,
        uniformity,
    ]
    return features


def compute_glcm_features(patch_discrete, n_levels, valid_mask=None, distances=[1]):


    patch = np.asarray(patch_discrete, dtype=np.int32)


    if valid_mask is None:
        valid_mask = np.ones_like(patch, dtype=bool)
    else:
        valid_mask = np.asarray(valid_mask, dtype=bool)

    if np.sum(valid_mask) < 2 or n_levels < 2:
        return [0.0] * 24


    offsets = [
        (0, 0, 1), (0, 1, 0), (1, 0, 0),
        (0, 1, 1), (0, 1, -1), (1, 0, 1), (1, 0, -1), (1, 1, 0), (1, -1, 0),
        (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)
    ]


    glcm_avg = np.zeros((n_levels, n_levels), dtype=np.float64)
    sz, sy, sx = patch.shape
    count = 0

    for dz, dy, dx in offsets:


        for z in range(sz):
            z2 = z + dz
            if z2 < 0 or z2 >= sz:
                continue
            for y in range(sy):
                y2 = y + dy
                if y2 < 0 or y2 >= sy:
                    continue
                for x in range(sx):
                    x2 = x + dx
                    if x2 < 0 or x2 >= sx:
                        continue
                    if not (valid_mask[z, y, x] and valid_mask[z2, y2, x2]):
                        continue
                    i_val = patch[z, y, x]
                    j_val = patch[z2, y2, x2]
                    if 0 <= i_val < n_levels and 0 <= j_val < n_levels:
                        glcm_avg[i_val, j_val] += 1
                        glcm_avg[j_val, i_val] += 1
        count += 1


    glcm_sum = np.sum(glcm_avg)
    if glcm_sum == 0:
        return [0.0] * 24
    glcm_avg = glcm_avg / glcm_sum


    px = np.sum(glcm_avg, axis=1)
    py = np.sum(glcm_avg, axis=0)


    i_indices = np.arange(n_levels)
    j_indices = np.arange(n_levels)
    mu_x = np.sum(i_indices * px)
    mu_y = np.sum(j_indices * py)
    var_x = np.sum((i_indices - mu_x) ** 2 * px)
    var_y = np.sum((j_indices - mu_y) ** 2 * py)


    joint_energy = np.sum(glcm_avg ** 2)
    joint_entropy = -np.sum(glcm_avg[glcm_avg > 0] * np.log2(glcm_avg[glcm_avg > 0]))


    diff_matrix = np.zeros((2 * n_levels - 1,), dtype=np.float64)
    sum_matrix = np.zeros((2 * n_levels - 1,), dtype=np.float64)
    for i in range(n_levels):
        for j in range(n_levels):
            if glcm_avg[i, j] > 0:
                diff_matrix[abs(i - j)] += glcm_avg[i, j]
                sum_matrix[i + j] += glcm_avg[i, j]

    diff_avg = np.sum(np.arange(len(diff_matrix)) * diff_matrix)
    diff_entropy = -np.sum(diff_matrix[diff_matrix > 0] * np.log2(diff_matrix[diff_matrix > 0]))
    diff_var = np.sum((np.arange(len(diff_matrix)) - diff_avg) ** 2 * diff_matrix)

    sum_avg = np.sum(np.arange(len(sum_matrix)) * sum_matrix)
    sum_entropy = -np.sum(sum_matrix[sum_matrix > 0] * np.log2(sum_matrix[sum_matrix > 0]))
    sum_var = np.sum((np.arange(len(sum_matrix)) - sum_avg) ** 2 * sum_matrix)


    contrast = np.sum(glcm_avg * (i_indices[:, None] - j_indices[None, :]) ** 2)


    if var_x > 0 and var_y > 0:
        correlation = np.sum(glcm_avg * (i_indices[:, None] - mu_x) * (j_indices[None, :] - mu_y)) / np.sqrt(var_x * var_y)
    else:
        correlation = 0.0


    autocorrelation = np.sum(glcm_avg * i_indices[:, None] * j_indices[None, :])


    cluster_prominence = np.sum(glcm_avg * ((i_indices[:, None] + j_indices[None, :] - mu_x - mu_y) ** 4))
    cluster_shade = np.sum(glcm_avg * ((i_indices[:, None] + j_indices[None, :] - mu_x - mu_y) ** 3))
    cluster_tendency = np.sum(glcm_avg * ((i_indices[:, None] + j_indices[None, :] - mu_x - mu_y) ** 2))


    id_val = np.sum(glcm_avg / (1 + np.abs(i_indices[:, None] - j_indices[None, :])))
    idm = np.sum(glcm_avg / (1 + (i_indices[:, None] - j_indices[None, :]) ** 2))


    diff_abs = np.abs(i_indices[:, None] - j_indices[None, :])
    idmn = np.sum(glcm_avg / (1.0 + diff_abs ** 2 / float(n_levels - 1) ** 2)) if n_levels > 1 else 0.0
    idn = np.sum(glcm_avg / (1.0 + diff_abs / float(n_levels - 1))) if n_levels > 1 else 0.0


    hx = -np.sum(px[px > 0] * np.log2(px[px > 0]))
    hy = -np.sum(py[py > 0] * np.log2(py[py > 0]))
    hxy = joint_entropy


    pxpy_outer = px[:, None] * py[None, :]
    valid_mask_glcm = glcm_avg > 0

    hxy1 = -np.sum(glcm_avg[valid_mask_glcm] * np.log2(pxpy_outer[valid_mask_glcm]))

    if hx > 0 and hy > 0:
        imc1 = (hxy - hxy1) / max(hx, hy)
    else:
        imc1 = 0.0


    imc2 = np.sqrt(max(0, 1 - np.exp(-2 * (hxy1 - hxy)))) if (hxy1 - hxy) > 0 else 0.0


    sum_squares = var_x


    max_prob = np.max(glcm_avg)


    joint_average = mu_x


    diff_sq = (i_indices[:, None] - j_indices[None, :]) ** 2
    inv_var_mask = diff_sq > 0
    inverse_variance = np.sum(glcm_avg[inv_var_mask] / diff_sq[inv_var_mask]) if np.any(inv_var_mask) else 0.0

    features = [
        autocorrelation,
        cluster_prominence,
        cluster_shade,
        cluster_tendency,
        contrast,
        correlation,
        diff_avg,
        diff_entropy,
        diff_var,
        joint_average,
        joint_energy,
        joint_entropy,
        imc1,
        imc2,
        id_val,
        idm,
        idmn,
        idn,
        sum_squares,
        max_prob,
        sum_avg,
        sum_entropy,
        sum_var,
        inverse_variance,
    ]

    return features


def compute_glrlm_features(patch_discrete, n_levels, valid_mask=None):


    patch = np.asarray(patch_discrete, dtype=np.int32)
    if valid_mask is None:
        valid_mask = patch >= 1
    else:
        valid_mask = np.asarray(valid_mask, dtype=bool) & (patch >= 1)
    if np.sum(valid_mask) < 2 or n_levels < 2:
        return [0.0] * 16

    offsets = [
        (0, 0, 1), (0, 1, 0), (1, 0, 0),
        (0, 1, 1), (0, 1, -1), (1, 0, 1), (1, 0, -1), (1, 1, 0), (1, -1, 0),
        (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)
    ]

    max_run_length = np.sum(valid_mask)
    glrlm_avg = np.zeros((n_levels, max_run_length + 1), dtype=np.float64)
    n_voxels_total = 0


    sz, sy, sx = patch.shape
    for dz, dy, dx in offsets:
        for z in range(sz):
            for y in range(sy):
                for x in range(sx):
                    if not valid_mask[z, y, x]:
                        continue
                    g = int(patch[z, y, x])
                    pz, py, px = z - dz, y - dy, x - dx
                    if (0 <= pz < sz and 0 <= py < sy and 0 <= px < sx and
                        valid_mask[pz, py, px] and patch[pz, py, px] == g):
                        continue
                    length = 0
                    cz, cy, cx = z, y, x
                    while (0 <= cz < sz and 0 <= cy < sy and 0 <= cx < sx and
                           valid_mask[cz, cy, cx] and patch[cz, cy, cx] == g):
                        length += 1
                        cz, cy, cx = cz + dz, cy + dy, cx + dx
                    if length:
                        glrlm_avg[g, length] += 1
        n_voxels_total += int(np.sum(valid_mask))

    total_runs = np.sum(glrlm_avg)
    if total_runs == 0:
        return [0.0] * 16


    glrlm = glrlm_avg / total_runs


    gray_level_marginal = np.sum(glrlm, axis=1)
    run_length_marginal = np.sum(glrlm, axis=0)

    gray_indices = np.arange(n_levels)
    run_indices = np.arange(max_run_length + 1)


    gray_level_nonuniformity = np.sum(gray_level_marginal ** 2) * total_runs
    gray_level_nonuniformity_norm = np.sum(gray_level_marginal ** 2)
    gray_level_variance = np.sum((gray_indices - np.sum(gray_indices * gray_level_marginal)) ** 2 * gray_level_marginal)

    high_gray_level_run_emphasis = np.sum(gray_indices[:, None] ** 2 * glrlm)
    long_run_emphasis = np.sum(run_indices[None, :] ** 2 * glrlm)
    long_run_high_gray = np.sum(gray_indices[:, None] ** 2 * run_indices[None, :] ** 2 * glrlm)
    long_run_low_gray = np.sum(run_indices[None, :] ** 2 * glrlm / (gray_indices[:, None] ** 2 + 1e-10))
    low_gray_level_run_emphasis = np.sum(glrlm / (gray_indices[:, None] ** 2 + 1e-10))

    run_entropy = -np.sum(glrlm[glrlm > 0] * np.log2(glrlm[glrlm > 0]))
    run_length_nonuniformity = np.sum(run_length_marginal ** 2) * total_runs
    run_length_nonuniformity_norm = np.sum(run_length_marginal ** 2)
    run_percentage = total_runs / n_voxels_total if n_voxels_total > 0 else 0.0
    run_variance = np.sum((run_indices - np.sum(run_indices * run_length_marginal)) ** 2 * run_length_marginal)

    short_run_emphasis = np.sum(glrlm / (run_indices[None, :] ** 2 + 1e-10))
    short_run_high_gray = np.sum(gray_indices[:, None] ** 2 * glrlm / (run_indices[None, :] ** 2 + 1e-10))
    short_run_low_gray = np.sum(glrlm / ((gray_indices[:, None] ** 2 + 1e-10) * (run_indices[None, :] ** 2 + 1e-10)))

    features = [
        gray_level_nonuniformity,
        gray_level_nonuniformity_norm,
        gray_level_variance,
        high_gray_level_run_emphasis,
        long_run_emphasis,
        long_run_high_gray,
        long_run_low_gray,
        low_gray_level_run_emphasis,
        run_entropy,
        run_length_nonuniformity,
        run_length_nonuniformity_norm,
        run_percentage,
        run_variance,
        short_run_emphasis,
        short_run_high_gray,
        short_run_low_gray,
    ]
    return features


def compute_glszm_features(patch_discrete, n_levels, valid_mask=None):


    patch = np.asarray(patch_discrete, dtype=np.int32)
    valid_mask = patch >= 0
    if np.sum(valid_mask) < 2 or n_levels < 2:
        return [0.0] * 16

    structure = generate_binary_structure(3, 3)
    max_zone_size = np.sum(valid_mask)
    glszm = np.zeros((n_levels, max_zone_size + 1), dtype=np.float64)

    for level in range(n_levels):
        level_mask = (patch == level) & valid_mask
        if np.sum(level_mask) == 0:
            continue
        labeled, num_zones = cc_label(level_mask, structure=structure)
        if num_zones == 0:
            continue


        zone_sizes = np.bincount(labeled.ravel())[1:]
        for zone_size in zone_sizes:
            if 1 <= zone_size <= max_zone_size:
                glszm[level, zone_size] += 1

    total_zones = np.sum(glszm)
    n_voxels = np.sum(valid_mask)
    if total_zones == 0:
        return [0.0] * 16

    glszm_norm = glszm / total_zones
    gray_marginal = np.sum(glszm_norm, axis=1)
    size_marginal = np.sum(glszm_norm, axis=0)

    gray_indices = np.arange(n_levels)
    size_indices = np.arange(max_zone_size + 1)

    gray_level_nonuniformity = np.sum(gray_marginal ** 2) * total_zones
    gray_level_nonuniformity_norm = np.sum(gray_marginal ** 2)
    gray_level_variance = np.sum((gray_indices - np.sum(gray_indices * gray_marginal)) ** 2 * gray_marginal)

    high_gray_zone_emphasis = np.sum(gray_indices[:, None] ** 2 * glszm_norm)
    large_area_emphasis = np.sum(size_indices[None, :] ** 2 * glszm_norm)
    large_area_high_gray = np.sum(gray_indices[:, None] ** 2 * size_indices[None, :] ** 2 * glszm_norm)
    large_area_low_gray = np.sum(size_indices[None, :] ** 2 * glszm_norm / (gray_indices[:, None] ** 2 + 1e-10))
    low_gray_zone_emphasis = np.sum(glszm_norm / (gray_indices[:, None] ** 2 + 1e-10))

    size_zone_nonuniformity = np.sum(size_marginal ** 2) * total_zones
    size_zone_nonuniformity_norm = np.sum(size_marginal ** 2)
    zone_entropy = -np.sum(glszm_norm[glszm_norm > 0] * np.log2(glszm_norm[glszm_norm > 0]))
    zone_percentage = total_zones / n_voxels if n_voxels > 0 else 0.0
    zone_variance = np.sum((size_indices - np.sum(size_indices * size_marginal)) ** 2 * size_marginal)

    small_area_emphasis = np.sum(glszm_norm / (size_indices[None, :] ** 2 + 1e-10))
    small_area_high_gray = np.sum(gray_indices[:, None] ** 2 * glszm_norm / (size_indices[None, :] ** 2 + 1e-10))
    small_area_low_gray = np.sum(glszm_norm / ((gray_indices[:, None] ** 2 + 1e-10) * (size_indices[None, :] ** 2 + 1e-10)))

    features = [
        gray_level_nonuniformity,
        gray_level_nonuniformity_norm,
        gray_level_variance,
        high_gray_zone_emphasis,
        large_area_emphasis,
        large_area_high_gray,
        large_area_low_gray,
        low_gray_zone_emphasis,
        size_zone_nonuniformity,
        size_zone_nonuniformity_norm,
        zone_entropy,
        zone_percentage,
        zone_variance,
        small_area_emphasis,
        small_area_high_gray,
        small_area_low_gray,
    ]
    return features


def compute_gldm_features(patch_discrete, n_levels, valid_mask=None):


    patch = np.asarray(patch_discrete, dtype=np.int32)
    valid_mask = patch >= 0
    if np.sum(valid_mask) < 2 or n_levels < 2:
        return [0.0] * 14


    offsets = []
    for dz in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dx in [-1, 0, 1]:
                if dz == 0 and dy == 0 and dx == 0:
                    continue
                offsets.append((dz, dy, dx))

    max_dependence = 27
    gldm = np.zeros((n_levels, max_dependence + 1), dtype=np.float64)


    padded = np.pad(patch, 1, constant_values=-1)
    padded_valid = np.pad(valid_mask, 1, constant_values=False)

    for z in range(patch.shape[0]):
        for y in range(patch.shape[1]):
            for x in range(patch.shape[2]):
                if not valid_mask[z, y, x]:
                    continue
                gray_level = patch[z, y, x]
                dependence_count = 0
                for dz, dy, dx in offsets:
                    nz, ny, nx = z + 1 + dz, y + 1 + dy, x + 1 + dx
                    if padded_valid[nz, ny, nx] and padded[nz, ny, nx] == gray_level:
                        dependence_count += 1
                gldm[gray_level, dependence_count + 1] += 1

    total_dependencies = np.sum(gldm)
    if total_dependencies == 0:
        return [0.0] * 14

    gldm_norm = gldm / total_dependencies
    gray_marginal = np.sum(gldm_norm, axis=1)
    dep_marginal = np.sum(gldm_norm, axis=0)

    gray_indices = np.arange(n_levels)
    dep_indices = np.arange(max_dependence + 1)

    dependence_entropy = -np.sum(gldm_norm[gldm_norm > 0] * np.log2(gldm_norm[gldm_norm > 0]))
    dependence_nonuniformity = np.sum(dep_marginal ** 2) * total_dependencies
    dependence_nonuniformity_norm = np.sum(dep_marginal ** 2)
    dependence_variance = np.sum((dep_indices - np.sum(dep_indices * dep_marginal)) ** 2 * dep_marginal)

    gray_level_nonuniformity = np.sum(gray_marginal ** 2) * total_dependencies
    gray_level_variance = np.sum((gray_indices - np.sum(gray_indices * gray_marginal)) ** 2 * gray_marginal)

    high_gray_emphasis = np.sum(gray_indices[:, None] ** 2 * gldm_norm)
    large_dependence_emphasis = np.sum(dep_indices[None, :] ** 2 * gldm_norm)
    large_dependence_high_gray = np.sum(gray_indices[:, None] ** 2 * dep_indices[None, :] ** 2 * gldm_norm)
    large_dependence_low_gray = np.sum(dep_indices[None, :] ** 2 * gldm_norm / (gray_indices[:, None] ** 2 + 1e-10))
    low_gray_emphasis = np.sum(gldm_norm / (gray_indices[:, None] ** 2 + 1e-10))

    small_dependence_emphasis = np.sum(gldm_norm / (dep_indices[None, :] ** 2 + 1e-10))
    small_dependence_high_gray = np.sum(gray_indices[:, None] ** 2 * gldm_norm / (dep_indices[None, :] ** 2 + 1e-10))
    small_dependence_low_gray = np.sum(gldm_norm / ((gray_indices[:, None] ** 2 + 1e-10) * (dep_indices[None, :] ** 2 + 1e-10)))

    features = [
        dependence_entropy,
        dependence_nonuniformity,
        dependence_nonuniformity_norm,
        dependence_variance,
        gray_level_nonuniformity,
        gray_level_variance,
        high_gray_emphasis,
        large_dependence_emphasis,
        large_dependence_high_gray,
        large_dependence_low_gray,
        low_gray_emphasis,
        small_dependence_emphasis,
        small_dependence_high_gray,
        small_dependence_low_gray,
    ]
    return features


def compute_ngtdm_features(patch_discrete, n_levels, valid_mask=None):


    patch = np.asarray(patch_discrete, dtype=np.int32)
    valid_mask = patch >= 0
    if np.sum(valid_mask) < 2 or n_levels < 2:
        return [0.0] * 5


    offsets = []
    for dz in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dx in [-1, 0, 1]:
                if dz == 0 and dy == 0 and dx == 0:
                    continue
                offsets.append((dz, dy, dx))

    padded = np.pad(patch, 1, constant_values=-1)
    padded_valid = np.pad(valid_mask, 1, constant_values=False)


    n_i = np.zeros(n_levels, dtype=np.float64)
    p_i = np.zeros(n_levels, dtype=np.float64)
    s_i = np.zeros(n_levels, dtype=np.float64)

    for z in range(patch.shape[0]):
        for y in range(patch.shape[1]):
            for x in range(patch.shape[2]):
                if not valid_mask[z, y, x]:
                    continue
                gray_level = patch[z, y, x]
                n_i[gray_level] += 1


                neighbor_sum = 0.0
                neighbor_count = 0
                for dz, dy, dx in offsets:
                    nz, ny, nx = z + 1 + dz, y + 1 + dy, x + 1 + dx
                    if padded_valid[nz, ny, nx]:
                        neighbor_sum += padded[nz, ny, nx]
                        neighbor_count += 1

                if neighbor_count > 0:
                    avg_neighbor = neighbor_sum / neighbor_count
                    s_i[gray_level] += abs(gray_level - avg_neighbor)

    n_voxels = np.sum(n_i)
    if n_voxels == 0:
        return [0.0] * 5

    p_i = n_i / n_voxels


    active = n_i > 0
    g = np.arange(n_levels, dtype=np.float64)[active]
    p = p_i[active]
    d = s_i[active]
    n_active = len(g)
    eps = 1e-10
    weighted_difference = np.sum(p * d)
    coarseness = 1.0 / (weighted_difference + eps) if weighted_difference > 0 else 0.0


    g_diff2 = (g[:, None] - g[None, :]) ** 2
    pair_product = p[:, None] * p[None, :]
    if n_active > 1:
        contrast = (np.sum(pair_product * g_diff2) /
                    (n_active * (n_active - 1)) * (np.sum(d) / n_voxels))
    else:
        contrast = 0.0


    weighted_levels = g * p
    busy_denom = np.sum(np.abs(weighted_levels[:, None] - weighted_levels[None, :]))
    busyness = weighted_difference / (busy_denom + eps) if busy_denom > 0 else 0.0


    denom = p[:, None] + p[None, :]
    numer = (p * d)[:, None] + (p * d)[None, :]
    mask_off_diag = ~np.eye(n_active, dtype=bool)
    complexity = (np.sum((np.abs(g[:, None] - g[None, :]) * numer / (denom + eps))
                         [mask_off_diag]) / n_voxels)
    strength_numer = np.sum(((p[:, None] + p[None, :]) * g_diff2)[mask_off_diag])
    strength = strength_numer / (np.sum(d) + eps) if np.sum(d) > 0 else 0.0

    features = [
        busyness,
        coarseness,
        complexity,
        contrast,
        strength,
    ]
    return features


def _masked_gradient_and_laplacian(intensities, valid_mask):


    arr = np.asarray(intensities, dtype=np.float64)
    valid = np.asarray(valid_mask, dtype=bool)
    grad2 = np.zeros_like(arr)
    lap = np.zeros_like(arr)
    for axis in range(3):
        before = np.zeros_like(arr)
        after = np.zeros_like(arr)
        has_before = np.zeros_like(valid)
        has_after = np.zeros_like(valid)
        cur = [slice(None)] * 3
        neigh = [slice(None)] * 3
        cur[axis] = slice(1, None)
        neigh[axis] = slice(None, -1)
        before[tuple(cur)] = arr[tuple(neigh)]
        has_before[tuple(cur)] = valid[tuple(neigh)]
        cur[axis] = slice(None, -1)
        neigh[axis] = slice(1, None)
        after[tuple(cur)] = arr[tuple(neigh)]
        has_after[tuple(cur)] = valid[tuple(neigh)]
        both = valid & has_before & has_after
        forward = valid & (~has_before) & has_after
        backward = valid & has_before & (~has_after)
        derivative = np.zeros_like(arr)
        derivative[both] = (after[both] - before[both]) / 2.0
        derivative[forward] = after[forward] - arr[forward]
        derivative[backward] = arr[backward] - before[backward]
        grad2 += derivative ** 2
        lap[both] += after[both] - 2 * arr[both] + before[both]
    return np.sqrt(grad2), lap


def compute_custom_local_features(patch_raw, patch_discrete, valid_mask_patch):


    valid = np.asarray(valid_mask_patch, dtype=bool)
    values = np.asarray(patch_raw, dtype=np.float64)[valid]
    n = len(values)
    if n == 0:
        return [0.0] * 11
    mean = float(np.mean(values))
    variance = float(np.var(values, ddof=0))
    std = float(np.sqrt(variance))
    disc_values = np.asarray(patch_discrete)[valid]
    _, counts = np.unique(disc_values, return_counts=True)
    prob = counts.astype(float) / n
    local_entropy = float(-np.sum(prob * np.log2(prob)))
    local_range = float(np.max(values) - np.min(values))
    local_contrast = std
    gradients, laplacian = _masked_gradient_and_laplacian(patch_raw, valid)
    g = gradients[valid]
    lv = laplacian[valid]
    grad_mean = float(np.mean(g))
    grad_var = float(np.var(g, ddof=0))
    lap_mean = float(np.mean(lv))
    lap_var = float(np.var(lv, ddof=0))
    tau_grad = float(np.mean(g) + np.std(g, ddof=0))
    edge_density = float(np.mean(g > tau_grad))

    lihi = float(variance / (mean ** 2 + 1e-6))
    return [mean, variance, local_entropy, local_range, local_contrast,
            grad_mean, grad_var, lap_mean, lap_var, edge_density, lihi]


def extract_local_feature_matrix(img_arr_hu, img_discrete, mask_arr, n_levels,
                                 window_size=3, voxel_volume=1.0):


    pad = window_size // 2


    img_padded = np.pad(img_arr_hu, pad, mode='constant', constant_values=0)
    disc_padded = np.pad(img_discrete, pad, mode='constant', constant_values=-1)
    mask_padded = np.pad(mask_arr, pad, mode='constant', constant_values=0)


    voxel_coords = np.argwhere(mask_arr > 0)
    n_voxels = len(voxel_coords)

    if n_voxels == 0:
        return np.zeros((0, 104)), np.zeros((0, 3), dtype=int)

    feature_matrix = np.zeros((n_voxels, 104), dtype=np.float64)

    for idx, (z, y, x) in enumerate(voxel_coords):

        z_p, y_p, x_p = z, y, x
        patch_raw = img_padded[z_p:z_p + window_size,
                               y_p:y_p + window_size,
                               x_p:x_p + window_size]
        patch_disc = disc_padded[z_p:z_p + window_size,
                                  y_p:y_p + window_size,
                                  x_p:x_p + window_size]
        patch_mask = mask_padded[z_p:z_p + window_size,
                                  y_p:y_p + window_size,
                                  x_p:x_p + window_size]


        valid_mask_patch = patch_mask > 0

        if np.sum(valid_mask_patch) == 0:
            continue


        values_raw = patch_raw[valid_mask_patch]


        values_disc = patch_disc[valid_mask_patch]


        first_order = compute_first_order_features(values_raw, values_disc,
                                                   voxel_volume=voxel_volume)


        glcm = compute_glcm_features(patch_disc, n_levels, valid_mask_patch)


        glrlm = compute_glrlm_features(patch_disc, n_levels, valid_mask_patch)


        glszm = compute_glszm_features(patch_disc, n_levels, valid_mask_patch)


        gldm = compute_gldm_features(patch_disc, n_levels, valid_mask_patch)


        ngtdm = compute_ngtdm_features(patch_disc, n_levels, valid_mask_patch)


        custom = compute_custom_local_features(patch_raw, patch_disc, valid_mask_patch)


        all_features = first_order + glcm + glrlm + glszm + gldm + ngtdm + custom
        if len(all_features) != 104:
            raise RuntimeError(f"Expected 104 features, found {len(all_features)}")
        feature_matrix[idx, :] = all_features

    return feature_matrix, voxel_coords


def standardize_features_within_scan(feature_matrix):


    standardized = feature_matrix.copy()
    for j in range(feature_matrix.shape[1]):
        col = feature_matrix[:, j]
        mean = np.mean(col)
        std = np.std(col, ddof=0)
        if std > 1e-10:
            standardized[:, j] = (col - mean) / std
        else:
            standardized[:, j] = 0.0
    return standardized


def impute_missing_values(feature_matrix):


    imputed = feature_matrix.copy()
    for j in range(feature_matrix.shape[1]):
        col = feature_matrix[:, j]
        finite_mask = np.isfinite(col)
        if np.sum(finite_mask) > 0:
            median_val = np.median(col[finite_mask])
        else:
            median_val = 0.0
        imputed[~finite_mask, j] = median_val
    return imputed


def perform_kmeans_clustering(feature_matrix, n_clusters=3, init='k-means++',
                               n_init=5, random_state=42):


    kmeans = KMeans(n_clusters=n_clusters, init=init, n_init=n_init,
                    random_state=random_state)
    labels = kmeans.fit_predict(feature_matrix)
    return labels, kmeans


def reconstruct_habitat_map(labels, voxel_coords, image_shape):


    habitat_map = np.zeros(image_shape, dtype=np.int32)
    for idx, (z, y, x) in enumerate(voxel_coords):
        habitat_map[z, y, x] = labels[idx] + 1
    return habitat_map


def compute_vhs(habitat_map, mask_arr, connectivity=26):


    structure = generate_binary_structure(3, 3)

    total_vat_volume = np.sum(mask_arr > 0)
    if total_vat_volume == 0:
        return 0.0, 0


    habitat_labels = np.unique(habitat_map[habitat_map > 0])
    if len(habitat_labels) == 0:
        return 0.0, 0

    sum_vmax_over_nk = 0.0
    n_components_total = 0


    for hab_label in habitat_labels:

        hab_mask = (habitat_map == hab_label) & (mask_arr > 0)

        if np.sum(hab_mask) == 0:
            continue


        labeled, n_components = cc_label(hab_mask, structure=structure)
        n_components_total += n_components

        if n_components == 0:
            continue


        component_sizes = np.bincount(labeled.ravel())
        component_sizes[0] = 0
        largest_component_size = np.max(component_sizes)


        sum_vmax_over_nk += largest_component_size / n_components


    vhs = 1.0 - (sum_vmax_over_nk / total_vat_volume)


    vhs = float(np.clip(vhs, 0.0, 1.0))

    return vhs, n_components_total


def compute_vhs_for_single_scan(image_path, mask_path, config=None,
                               target_spacing_override=None):


    if config is None:
        config = Config()


    eff_target = target_spacing_override if target_spacing_override is not None else config.TARGET_SPACING
    image_sitk, mask_sitk = load_nifti(image_path, mask_path)
    if image_sitk.GetDimension() != 3 or mask_sitk.GetDimension() != 3:
        raise ValueError("CT image and VAT mask must both be 3-D")
    if image_sitk.GetSize() != mask_sitk.GetSize():
        raise ValueError("CT image and VAT mask have different voxel-grid sizes")
    image_resampled, mask_resampled = resample_to_isotropic(
        image_sitk, mask_sitk, target_spacing=eff_target)

    import SimpleITK as sitk
    img_arr = sitk.GetArrayFromImage(image_resampled)
    mask_arr = sitk.GetArrayFromImage(mask_resampled)


    mask_arr = retain_largest_connected_component(mask_arr)


    img_clipped = clip_intensities(
        img_arr, mask_arr,
        lower_percentile=config.INTENSITY_CLIP_LOWER_PERCENTILE,
        upper_percentile=config.INTENSITY_CLIP_UPPER_PERCENTILE)


    img_discrete, n_levels = discretize_to_bins(
        img_clipped, mask_arr, bin_width_hu=config.BIN_WIDTH_HU)


    feature_matrix, voxel_coords = extract_local_feature_matrix(
        img_clipped, img_discrete, mask_arr, n_levels,
        window_size=config.WINDOW_SIZE,
        voxel_volume=float(np.prod(image_resampled.GetSpacing())))


    feature_matrix = impute_missing_values(feature_matrix)


    feature_matrix_std = standardize_features_within_scan(feature_matrix)
    if len(feature_matrix_std) < config.N_CLUSTERS:
        raise ValueError("Fewer valid VAT voxels than requested clusters")


    labels, kmeans = perform_kmeans_clustering(
        feature_matrix_std,
        n_clusters=config.N_CLUSTERS,
        init=config.KMEANS_INIT,
        n_init=config.KMEANS_N_INIT,
        random_state=config.KMEANS_RANDOM_STATE)


    habitat_map = reconstruct_habitat_map(labels, voxel_coords, img_arr.shape)


    vhs, n_components = compute_vhs(habitat_map, mask_arr, connectivity=config.CONNECTIVITY)

    return vhs, habitat_map, feature_matrix


def compute_percent_delta_vhs(vhs_pre, vhs_post):


    vhs_pre_effective = max(vhs_pre, Config.VHS_CLIP_LOWER)
    return (vhs_pre_effective - vhs_post) / vhs_pre_effective * 100.0


def process_cohort(pre_image_dir, post_image_dir, pre_mask_dir, post_mask_dir,
                    output_path=None, config=None):


    if config is None:
        config = Config()


    pre_images = sorted([f for f in os.listdir(pre_image_dir) if f.endswith(('.nii', '.nii.gz'))])

    results = []
    for case_number, pre_img_file in enumerate(pre_images, start=1):
        case_id = f"case_{case_number:04d}"


        pre_img_path = os.path.join(pre_image_dir, pre_img_file)
        post_img_path = os.path.join(post_image_dir, pre_img_file)
        pre_mask_path = os.path.join(pre_mask_dir, pre_img_file)
        post_mask_path = os.path.join(post_mask_dir, pre_img_file)


        if not all(os.path.exists(p) for p in [post_img_path, pre_mask_path, post_mask_path]):
            print(f"{case_id}: skipped (missing inputs)")
            continue

        try:
            print(f"Processing {case_id}")


            native_pre  = get_native_min_spacing(pre_img_path)
            native_post = get_native_min_spacing(post_img_path)
            shared_target = min(native_pre, native_post)

            vhs_pre, _, _  = compute_vhs_for_single_scan(pre_img_path, pre_mask_path, config,
                                                        target_spacing_override=shared_target)
            vhs_post, _, _ = compute_vhs_for_single_scan(post_img_path, post_mask_path, config,
                                                         target_spacing_override=shared_target)
            pct_delta = compute_percent_delta_vhs(vhs_pre, vhs_post)

            results.append({
                'CaseID': case_id,
                'VHSpre': vhs_pre,
                'VHSpost': vhs_post,
                'pctDelta_VHS': pct_delta,
            })

        except Exception as e:
            print(f"{case_id}: failed ({type(e).__name__})")
            continue

    results_df = pd.DataFrame(results)

    if output_path:
        if output_path.endswith('.xlsx'):
            results_df.to_excel(output_path, index=False)
        else:
            results_df.to_csv(output_path, index=False)
        print("Results saved")

    return results_df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="IBD-ITH VHS Pipeline")
    parser.add_argument("--pre-image-dir", required=True, help="Directory with pre-induction CT images")
    parser.add_argument("--post-image-dir", required=True, help="Directory with post-induction CT images")
    parser.add_argument("--pre-mask-dir", required=True, help="Directory with pre-induction VAT masks")
    parser.add_argument("--post-mask-dir", required=True, help="Directory with post-induction VAT masks")
    parser.add_argument("--output", default="vhs_results.xlsx", help="Output file path")
    parser.add_argument("--target-spacing", type=float, default=None, help="Target isotropic spacing (mm)")
    parser.add_argument("--bin-width", type=int, default=25, help="Gray-level bin width (HU)")
    parser.add_argument("--n-clusters", type=int, default=3, help="Number of K-means clusters")
    parser.add_argument("--random-state", type=int, default=42, help="K-means random seed")
    parser.add_argument("--n-init", type=int, default=5, help="K-means n_init")

    args = parser.parse_args()

    config = Config()
    config.TARGET_SPACING = args.target_spacing
    config.BIN_WIDTH_HU = args.bin_width
    config.N_CLUSTERS = args.n_clusters
    config.KMEANS_RANDOM_STATE = args.random_state
    config.KMEANS_N_INIT = args.n_init

    results = process_cohort(
        args.pre_image_dir, args.post_image_dir,
        args.pre_mask_dir, args.post_mask_dir,
        output_path=args.output, config=config)

    print(f"\nProcessed {len(results)} patients.")

