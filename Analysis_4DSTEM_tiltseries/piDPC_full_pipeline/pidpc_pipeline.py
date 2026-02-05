# pidpc_pipeline.py
# Author: Shahar Seifer (Elbaum Lab), Weizmann Institute of Science
# GPLv3 license
# Python translation supported by M365 Copilot (Feb 2026)
#
# Pipeline:
#   DEFAULT (legacy HDF5, from Arina camera):
#     load_legacy_h5_to_mat(...)  --> mat: (Nk_y, Nk_x, Ny*Nx)
#   RAW mode (when --raw is set):
#     load_4dstem_h5_to_mat_raw(Ny, Nx, Nky, Nkx)
#
#pip install numpy scipy h5py mrcfile scikit-image opencv
#Examples of use:
#python pidpc_pipeline.py your.h5
#python pidpc_pipeline.py raw_rows.h5 --raw
#python pidpc_pipeline.py data/scan_0102.h5 --out-base out/scan_0102 --nx 128 --ny 96
#python pidpc_pipeline.py your.h5 --out-base out/pi --px 0.5 --py 0.5  (if pixel size is 0.5 Angstrom)

from __future__ import annotations
import os
import math
from dataclasses import dataclass
import numpy as np
import h5py
import mrcfile
from scipy import fft as spfft
from scipy.ndimage import gaussian_filter, shift as nd_shift
from scipy.optimize import least_squares
from scipy.signal import fftconvolve
from scipy.signal import correlate2d

from numpy.linalg import norm
import cv2

# ============================================================
# Helpers (unambiguous, symmetric) — use everywhere for shifts
# ============================================================
def imtranslate_xy(image: np.ndarray, dx: float, dy: float, order: int = 1) -> np.ndarray:

    return nd_shift(image, shift=(float(dx), float(dy)), order=order,
                    mode='constant', cval=0.0, prefilter=False)

def shift_to_center(image: np.ndarray, Xc: float, Yc: float, order: int = 1) -> np.ndarray:
    """
    Shift the image so a feature at (Xc, Yc) moves to the origin.
    That is: apply (dx, dy) = (-Xc, -Yc) using imtranslate_xy.
    """
    return imtranslate_xy(image, dx=-float(Xc), dy=-float(Yc), order=order)

def compose_shifts(*pairs: tuple[float, float]) -> tuple[float, float]:
    """
    Vector-add multiple (dx, dy) shifts.
    Example: compose_shifts((-Xc, -Yc), dither, (-r[0], -r[1])) → total (dx, dy).
    """
    dx = sum(p[0] for p in pairs)
    dy = sum(p[1] for p in pairs)
    return (dx, dy)


# ============================================================
# Utilities
# ============================================================


def fft_interpolate(arr, new_shape):
    """N-D Fourier interpolation (zero-fill/up, truncate/down, global scale) – MATLAB-like."""
    arr = np.asarray(arr)
    old_shape = np.array(arr.shape, dtype=int)
    new_shape = np.array(new_shape, dtype=int)
    if np.array_equal(new_shape, old_shape):
        return arr.copy()
    if np.any(new_shape <= 0):
        raise ValueError("fft_interpolate :: bad size")

    k = np.fft.fftshift(np.fft.fftn(arr))
    up_axes = new_shape > old_shape

    center = np.floor_divide(old_shape, 2) + 1
    lo = np.floor(center - new_shape / 2).astype(int)
    hi = lo + new_shape - 1
    lo[up_axes] = 1
    hi[up_axes] = old_shape[up_axes]

    src_slices = tuple(slice(max(0, lo[d]-1), min(old_shape[d], hi[d])) for d in range(arr.ndim))
    k_cropped = k[src_slices]
    k_new = np.zeros(new_shape, dtype=k.dtype)

    cropped_shape = np.array(k_cropped.shape, dtype=int)
    centerd = np.floor_divide(new_shape, 2)
    lod = np.ceil(centerd - new_shape / 2).astype(int) + 1
    hid = lod + new_shape - 1
    centeru = np.floor_divide(new_shape, 2) + 1 + (new_shape % 2)
    lou = np.ceil(centeru - cropped_shape / 2).astype(int)
    hiu = lou + cropped_shape - 1

    dst_lo = lod.copy(); dst_hi = hid.copy()
    dst_lo[new_shape > old_shape] = lou[new_shape > old_shape]
    dst_hi[new_shape > old_shape] = hiu[new_shape > old_shape]
    dst_slices = tuple(slice(max(0, dst_lo[d]-1), min(new_shape[d], dst_hi[d])) for d in range(arr.ndim))
    k_new[dst_slices] = k_cropped

    scale = float(np.prod(new_shape / old_shape))
    out = scale * np.fft.ifftn(np.fft.ifftshift(k_new))
    return out.real if np.isrealobj(arr) else out


def _q_grids_96_matlab_exact(flip_kx: bool = False, flip_ky: bool = False):
    """
    MATLAB:
        [qY, qX] = meshgrid((1:96)-(1+96)/2, (1:96)-(1+96)/2) with 'xy' order.
    IMPORTANT: as in MATLAB naming, qY holds the X-grid (columns), qX holds the Y-grid (rows).
    """
    gv = np.arange(1, 96 + 1) - (1 + 96) / 2.0
    Xg, Yg = np.meshgrid(gv, gv, indexing='xy')
    qY, qX = Xg, Yg  # crossed names on purpose (MATLAB-compatible)
    if flip_kx: qY = -qY
    if flip_ky: qX = -qX
    return qY, qX


def _bin192_to_96_like_matlab(im: np.ndarray) -> np.ndarray:
    """MATLAB mask+reshape (column-major) downsampling 192→96 via odd rows/cols."""
    assert im.shape[0] == 192 and im.shape[1] == 192
    mask = np.zeros_like(im, dtype=bool)
    mask[0::2, 0::2] = True
    vec = np.ravel(im, order='F')[np.ravel(mask, order='F')]
    return np.reshape(vec, (96, 96), order='F')


# ============================================================
# H5 Loading
# ============================================================
@dataclass
class H5Meta:
    Ny: int; Nx: int; NKy: int; NKx: int; dtype: np.dtype; row_keys: list[str]


def inspect_h5_structure_raw(h5_path: str) -> H5Meta:
    with h5py.File(h5_path, 'r') as f:
        rows = list(f.keys())
        Ny = len(rows)
        sample = f[rows[0]]
        Nx, NKy, NKx = sample.shape
        dtype = sample.dtype
    return H5Meta(Ny=Ny, Nx=Nx, NKy=NKy, NKx=NKx, dtype=dtype, row_keys=rows)


def load_4dstem_h5_to_mat_raw(h5_path: str) -> tuple[np.ndarray, int, int]:
    """RAW reader: top-level row groups; return (H,W,T) in MATLAB scan order (y desc, x asc)."""
    meta = inspect_h5_structure_raw(h5_path)
    Ny, Nx, NKy, NKx = meta.Ny, meta.Nx, meta.NKy, meta.NKx
    mat = np.empty((NKy, NKx, Ny * Nx), dtype=np.float32)
    t = 0
    with h5py.File(h5_path, 'r') as f:
        for y_idx in range(Ny - 1, -1, -1):
            row = f[meta.row_keys[y_idx]][...]  # (Nx, NKy, NKx)
            for x_idx in range(Nx):
                frame = row[x_idx, :, :]
                mat[:, :, t] = frame
                t += 1
    return mat, Nx, Ny


def _process_dataset_legacy(dset, start_index: int, array: np.ndarray) -> int:
    image_index = start_index
    for i in range(dset.shape[0]):
        array[image_index] = dset[i].astype(array.dtype, copy=False)
        image_index += 1
    return image_index


def load_legacy_h5_to_mat(h5_path: str,
                          nx_hint: int | None = None,
                          ny_hint: int | None = None) -> tuple[np.ndarray, int, int]:
    """
    Legacy reader under /entry/data (single dataset /entry/data/data or multiple datasets).
    Returns (H,W,T) in MATLAB scan order (y desc, x asc) and (Nx,Ny).
    """
    with h5py.File(h5_path, "r") as f:
        entry = f["entry"]
        if isinstance(entry.get("data"), h5py.Dataset):
            dset = entry["data"]
            nimages, H, W = dset.shape[0], dset.shape[1], dset.shape[2]
            array = np.empty((nimages, H, W), dtype=dset.dtype)
            _process_dataset_legacy(dset, 0, array)
        else:
            grp = entry["data"]
            nimages, H, W, dtype = 0, None, None, None
            for name in grp:
                d = grp[name]
                nimages += d.shape[0]; H, W = d.shape[1], d.shape[2]; dtype = d.dtype
            array = np.empty((nimages, H, W), dtype=dtype)
            image_index = 0
            for name in grp:
                image_index = _process_dataset_legacy(grp[name], image_index, array)

        if array.dtype == np.uint32:
            print("Dataset is uint32 but will be converted to uint16")
            array = array.astype(np.uint16, copy=False)

        if nx_hint and ny_hint and nx_hint * ny_hint == array.shape[0]:
            Nx, Ny = int(nx_hint), int(ny_hint)
        else:
            Nx = int(math.sqrt(array.shape[0])); Ny = int(array.shape[0] / Nx)
        print(f"Dataset [{Nx} {Ny} {W} {H} @ {array.dtype}]")

        mat = np.transpose(array.astype(np.float32, copy=False), (1, 2, 0))  # (H,W,T)

        T = mat.shape[2]; assert T == Nx * Ny, f"Frame count {T} != Nx*Ny {Nx*Ny}"
        perm = []
        for y_desc in range(Ny - 1, -1, -1):
            y_rm = (Ny - 1) - y_desc
            for x in range(Nx):
                perm.append(y_rm * Nx + x)
        mat = mat[:, :, np.asarray(perm, dtype=int)]
    return mat, Nx, Ny


# ============================================================
# Image shift detection
# ============================================================

def _normxcorr2_full(template: np.ndarray, image: np.ndarray) -> np.ndarray:
    T = np.asarray(template, dtype=np.float64)
    I = np.asarray(image,    dtype=np.float64)
    Tm = T - T.mean()
    t_norm = np.sqrt(np.sum(Tm * Tm)) + np.finfo(np.float64).eps
    num = fftconvolve(I, np.flipud(np.fliplr(Tm)), mode='full')
    ones = np.ones(T.shape, dtype=np.float64)
    sum_I  = fftconvolve(I, ones, mode='full')
    sum_I2 = fftconvolve(I * I, ones, mode='full')
    n = T.size
    mean_I = sum_I / n
    var_I  = np.maximum(sum_I2 - (sum_I * sum_I) / n, 0.0)
    std_I  = np.sqrt(var_I) + np.finfo(np.float64).eps
    return ((num - mean_I * np.sum(Tm)) / (std_I * t_norm)).T

from scipy.signal import correlate2d


def r_mn(Imagem: np.ndarray, Imagen: np.ndarray, shift_limit: int = 50, do_filt: int = 1) -> np.ndarray:


    im1 = np.asarray(Imagem, dtype=np.float64)
    im2 = np.asarray(Imagen, dtype=np.float64)

    if do_filt == 1:
        # same band-pass you used (HP then LP)
        im1 = gaussian_filter(im1 - gaussian_filter(im1, 30), 3)
        im2 = gaussian_filter(im2 - gaussian_filter(im2, 30), 3)

    H, W = im1.shape
    # MATLAB template bounds (your code): 30% central crop
    tempx = int(np.floor(0.25 * H))
    tempy = int(np.floor(0.25 * W))
    tempux = H - tempx
    tempuy = W - tempy
    # This produces a (Th, Tw) template (rows: tempx..tempux-1, cols: tempy..tempuy-1)
    templ = im1[tempx:tempux, tempy:tempuy]

    #import matplotlib.pyplot as plt
    #plt.figure()
    #plt.subplot(1, 2, 1);
    #plt.imshow(templ);
    #plt.title("Template")
    #plt.subplot(1, 2, 2);
    #plt.imshow(im2);
    #plt.title("Image 2")
    #plt.show()

    # Full normalized cross-correlation, option without opencv (requires further modifications: yoffset = ypeak - tempuy; xoffset = xpeak - tempux)
    #corr = _normxcorr2_full(templ, im2)

    corr = cv2.matchTemplate(templ.astype(np.float32), im2.astype(np.float32), cv2.TM_CCOEFF_NORMED)

    xpeak, ypeak = np.unravel_index(np.argmax(np.abs(corr)), corr.shape)

    center_im2_y = (im2.shape[1] - 1) / 2.0
    center_im2_x = (im2.shape[0] - 1) / 2.0

    center_templ_on_im2_y = ypeak + (templ.shape[1] - 1) / 2.0
    center_templ_on_im2_x = xpeak + (templ.shape[0] - 1) / 2.0

    yoffset = center_templ_on_im2_y - center_im2_y
    xoffset = center_templ_on_im2_x - center_im2_x

    # Sub pixel correction
    ys0 = ypeak - 7; ys1 = ypeak + 8   # inclusive in MATLAB
    xs0 = xpeak - 7; xs1 = xpeak + 8
    sample16 = np.zeros((16, 16), dtype=np.float64)
    # compute valid source ranges
    sx0 = max(0, xs0); sx1 = min(corr.shape[0]-1, xs1)
    sy0 = max(0, ys0); sy1 = min(corr.shape[1]-1, ys1)
    # map to destination indices (remember: we swap to match MATLAB's 'corr(xrange, yrange)')
    # dst rows correspond to x-range; dst cols correspond to y-range
    dst_x0 = max(0, sx0 - xs0); dst_x1 = dst_x0 + (sx1 - sx0 + 1)
    dst_y0 = max(0, sy0 - ys0); dst_y1 = dst_y0 + (sy1 - sy0 + 1)

    block = corr[sx0:sx1 + 1, sy0:sy1 + 1]

    # place into 16x16 canvas
    sample16[dst_x0:dst_x1, dst_y0:dst_y1] = block

    # ---- MATLAB fftInterpolate(sample16, [512 512]) ----
    Intsample16 = fft_interpolate(sample16, (512, 512))

     # Find sub-pixel peak (MATLAB: [xpeak2,ypeak2] = ind2sub(size(...), maxIdx))
    maxIndex2 = int(np.argmax(np.abs(Intsample16)))
    xpeak2, ypeak2 = np.unravel_index(maxIndex2, Intsample16.shape)  # row, col
    yoffset2 = yoffset+ (ypeak2 - 256.0+32) * (16.0 / 512.0)
    xoffset2 = xoffset + (xpeak2 - 256.0+32) * (16.0 / 512.0)

    # Reject unreasonable estimates (MATLAB also guards by shift_limit)
    if abs(xoffset2) > shift_limit or abs(yoffset2) > shift_limit:
        return np.array([0.0, 0.0], dtype=float)

    # Return (dx, dy) in columns, rows
    return np.array([xoffset2, yoffset2], dtype=float)


def deshift(img1, img2, img3, img4):
    """Pairwise NCC + pseudoinverse; rx from second comps, ry from first comps (like MATLAB)."""
    if all(np.sum(np.abs(img)) > 0 for img in [img1, img2, img3, img4]):
        r21 = r_mn(img2, img1); r12 = -r21
        r31 = r_mn(img3, img1); r13 = -r31
        r41 = r_mn(img4, img1); r14 = -r41
        r32 = r_mn(img3, img2); r23 = -r32
        r42 = r_mn(img4, img2); r24 = -r42
        r43 = r_mn(img4, img3); r34 = -r43
        A = r21 + r31 + r41
        B = r12 + r32 + r42
        C = r13 + r23 + r43
        D = r14 + r24 + r34
        M = np.array([[ 3, -1, -1, -1],
                      [-1,  3, -1, -1],
                      [-1, -1,  3, -1],
                      [-1, -1, -1,  3]], dtype=float)
        invM = np.linalg.pinv(M)
        rx = invM @ np.array([A[1], B[1], C[1], D[1]])
        ry = invM @ np.array([A[0], B[0], C[0], D[0]])
        return np.array([rx[0], ry[0]]), np.array([rx[1], ry[1]]), np.array([rx[2], ry[2]]), np.array([rx[3], ry[3]])
    z = np.array([0.0, 0.0])
    return z, z, z, z


# ============================================================
# Beam tracking (fit disk shift vs scan position)
# ============================================================
def _moving_average_indices(ind: int, half: int, span_left: int, span_right: int, total: int) -> np.ndarray:
    a = max(max(ind - half, ind - span_left + 1), 1)
    b = min(min(ind + half, ind + span_right), total)
    return np.arange(a, b + 1, dtype=int)


def arina_trackbeam(mat: np.ndarray, nX: int, nY: int,
                    flip_kx: bool = False, flip_ky: bool = False) -> np.ndarray:
    """function_Arina_trackbeam (96x96 grid), COM swap (Xd from qY, Yd from qX). Returns [a0 a1 a2 b0 b1 b2]."""
    H, W, _ = mat.shape
    factor_bin = 2 if (H == 192 and W == 192) else 1
    qY, qX = _q_grids_96_matlab_exact(flip_kx=flip_kx, flip_ky=flip_ky)

    Xd_list, Yd_list, Xp_list, Yp_list = [], [], [], []
    movingavg_halfsize = int(round(20 * nX / 1024)) or 1

    Xp, Yp = 1, nY
    veclength = mat.shape[2]
    for ind in range(1, veclength + 1):
        if ind % movingavg_halfsize == 0:
            vavg = _moving_average_indices(ind, movingavg_halfsize, (Xp - 1), (nX - Xp), veclength)
            im_raw = mat[:, :, vavg - 1].mean(axis=2) if len(vavg) > 1 else mat[:, :, ind - 1]
        else:
            im_raw = mat[:, :, ind - 1]
        im = _bin192_to_96_like_matlab(im_raw) if factor_bin == 2 else im_raw

        meanv = im.mean()
        mask = im > meanv
        if np.any(mask):
            mw = (im * mask) / float((im * mask).sum())
            Xd = (qY[mask] * mw[mask]).sum()
            Yd = (qX[mask] * mw[mask]).sum()
            Xd_list.append(Xd); Yd_list.append(Yd)
            Xp_list.append(Xp); Yp_list.append(Yp)

        if Xp == nX:
            Xp = 1
            Yp -= 1
            if Yp == 0:
                break
        else:
            Xp += 1

    Xd_arr = np.array(Xd_list) - np.mean(Xd_list)
    Yd_arr = np.array(Yd_list) - np.mean(Yd_list)
    XY = np.column_stack([np.array(Xp_list), np.array(Yp_list)])

    def residuals(P):
        a0, a1, a2, b0, b1, b2 = P
        pred_x = a0 + a1 * XY[:, 0] + a2 * XY[:, 1]
        pred_y = b0 + b1 * XY[:, 0] + b2 * XY[:, 1]
        return np.concatenate([pred_x - Xd_arr, pred_y - Yd_arr])

    P0 = np.zeros(6, dtype=float)
    res = least_squares(residuals, P0, method='lm')
    return res.x  # [a0 a1 a2 b0 b1 b2]


# ============================================================
# Rings (sectorization; NO theta modulo; POST-REMAP for anchored reversal)
# ============================================================
def arina_tomo_rings_v9(mat: np.ndarray, out_base: str, nX: int, nY: int,
                        flip_kx: bool = False, flip_ky: bool = False,
                        sector_rot: int = 0,
                        theta_mirror: bool = False,
                        digest: bool = False) -> None:
    """
    Rings stage with MATLAB-parity geometry and **no theta modulo**.
    Writes 16 sector files:
      _sect1..8  = inside BF per wedge
      _sect9..16 = outside BF per wedge
    NEW: if write_sector_masks=True, also writes 16 binary masks (float32 0/1):
      _sectmask1..8  (inside BF), _sectmask9..16 (outside BF)
    """
    H, W, T = mat.shape
    assert H in (96, 192) and W in (96, 192), "Expected 96x96 or 192x192 detector."
    factor_bin = 2 if (H == 192 and W == 192) else 1

    qY, qX = _q_grids_96_matlab_exact(flip_kx=flip_kx, flip_ky=flip_ky)
    q = np.sqrt(qX**2 + qY**2)

    theta = np.angle(qX + 1j * qY) + np.pi

    if theta_mirror:
        theta = (2.0 * np.pi - theta)

        # Build forward wedges (k=1..8), shift by sector_rot; split if (high>2π)
    mask_section = np.zeros((96, 96, 8), dtype=bool)
    rot = (sector_rot % 8) * (2.0 * np.pi / 8.0)
    inrad = (q <= 50)
    for tsec in range(8):
        k = (tsec ) % 8
        low  = 2.0 * np.pi * (k) / 8.0 + rot
        high = 2.0 * np.pi * (k+1) / 8.0 + rot
        sel = (theta >= low) & (theta < high) & inrad
        mask_section[:, :, tsec] = sel


    # BF COM from average
    im_grand = np.zeros((96, 96), dtype=float)
    frames_considered = min(T, nX * nY)
    for ind in range(frames_considered):
        im_raw = mat[:, :, ind]
        im = _bin192_to_96_like_matlab(im_raw) if factor_bin == 2 else im_raw
        im_grand += im

    midv = 0.5 * (im_grand.max() + im_grand.min())
    mask_bf = im_grand > midv
    m_weight = (im_grand * mask_bf) / float((im_grand * mask_bf).sum())

    Yd0 = (qY[mask_bf] * m_weight[mask_bf]).sum()
    Xd0 = (qX[mask_bf] * m_weight[mask_bf]).sum()

    radius_disk = int(round(np.sqrt(mask_bf.sum() / np.pi)))
    mask_com = (q <= (radius_disk + 1))

    # Fit tilt surface; zero intercepts
    a0, a1, a2, b0, b1, b2 = arina_trackbeam(mat, nX=nX, nY=nY, flip_kx=flip_kx, flip_ky=flip_ky)
    a0 = 0.0; b0 = 0.0

    def surf_fit(xy: tuple[int, int]) -> tuple[float, float]:
        x, y = xy
        xd = a0 + a1 * x + a2 * y  # x/cols
        yd = b0 + b1 * x + b2 * y  # y/rows
        return xd, yd

    # Accumulate inside/outside sums AND (optionally) write masks
    scan_grandt = np.zeros((nX, nY, 16), dtype=np.float32)

    # ------------------------------------------------------------
    # Digest (first 200 frames; recenter by Xd0,Yd0 ONLY; q<=50 disk)
    # ------------------------------------------------------------
    if digest:
        frames_to_report = min(200, T)
        report_sums_in = np.zeros(8, dtype=np.float64)
        report_sums_out = np.zeros(8, dtype=np.float64)
        for ind in range(frames_to_report):
            im_raw = mat[:, :, ind]
            im = _bin192_to_96_like_matlab(im_raw) if factor_bin == 2 else im_raw
            # Recenter by (Xd0, Yd0) ONLY — no tilt terms for a fair DIGEST comparison
            im0 = shift_to_center(im, Xc=Xd0, Yc=Yd0, order=1)
            # Integrate per wedge (inside q<=50 vs outside)
            for tsec in range(8):
                inside_mask = (mask_section[:, :, tsec] & mask_com)
                outside_mask = (mask_section[:, :, tsec] & (~mask_com))
                report_sums_in[tsec] += im0[inside_mask].sum()
                report_sums_out[tsec] += im0[outside_mask].sum()
        area_in = [int(np.count_nonzero(mask_section[:, :, t] & mask_com)) for t in range(8)]
        area_out = [int(np.count_nonzero(mask_section[:, :, t] & (~mask_com))) for t in range(8)]
        print("\n=== MATLAB-PARITY SECTOR REPORT (first %d frames; recentered by Xd0,Yd0; q<=50) ===" % frames_to_report)
        print("Xd0: %.4f   Yd0: %.4f" % (Xd0, Yd0))
        for tsec in range(8):
            print("  Sector %d:  S_in=%10.3f   S_out=%10.3f   (area_in=%d, area_out=%d)"
                  % (tsec + 1, report_sums_in[tsec], report_sums_out[tsec], area_in[tsec], area_out[tsec]))
        print("=== END REPORT ===\n")
        # Save 16 masks: 1..8 inside, 9..16 outside
        for tsec in range(8):
            mask_in  = (mask_section[:, :, tsec] & mask_com).astype(np.float32)
            mask_out = (mask_section[:, :, tsec] & (~mask_com)).astype(np.float32)
            with mrcfile.new(f"{out_base}_sectmask{tsec+1}.mrc", overwrite=True) as m:
                m.set_data(mask_in.T)
            with mrcfile.new(f"{out_base}_sectmask{tsec+9}.mrc", overwrite=True) as m:
                m.set_data(mask_out.T)
        # Quick area printout
        print("SECTOR MASK AREAS (pixels):")
        for tsec in range(8):
            a_in  = int(np.count_nonzero(mask_section[:, :, tsec] & mask_com))
            a_out = int(np.count_nonzero(mask_section[:, :, tsec] & (~mask_com)))
            print(f"  sector {tsec+1}: inside={a_in}, outside={a_out}")


    posX, posY = 1, nY
    for ind in range(T):
        im_raw = mat[:, :, ind]
        im = _bin192_to_96_like_matlab(im_raw) if factor_bin == 2 else im_raw
        xd_pred, yd_pred = surf_fit((posX, posY))

        # Recentering (symmetric signs via helper)
        Xc = Xd0 + xd_pred
        Yc = Yd0 + yd_pred
        im_corr = shift_to_center(im, Xc, Yc, order=1)

        for tsec in range(8):
            inside  = (mask_section[:, :, tsec] & mask_com)
            outside = (mask_section[:, :, tsec] & (~mask_com))
            scan_grandt[posX - 1, posY - 1, tsec    ] += im_corr[inside].sum()   # 1..8 (inside)
            scan_grandt[posX - 1, posY - 1, tsec + 8] += im_corr[outside].sum()  # 9..16 (outside)

        if posX == nX:
            posX = 1
            if posY > 1: posY -= 1
            else:        posY = nY
        else:
            posX += 1

    # Write 16 sector MRCs as (Y, X)
    for k in range(16):
        path = f"{out_base}_sect{k+1}.mrc"
        with mrcfile.new(path, overwrite=True) as m:
            m.set_data(scan_grandt[:, :, k].T)



def _mk_synthetic(H=128, W=128):
    Y, X = np.mgrid[0:H, 0:W]
    img = np.exp(-((X - W/2)**2 + (Y - H/2)**2) / (2*3.0**2))    # Gaussian blob
    img += ((X > 57) & (X < 60) & (Y > 57) & (Y < 60)).astype(float)  # small square
    return img.astype(np.float32)

def _peak_offset_int(a, b):
    # integer-precision correlation peak (good enough as a residual sanity check)
    c = correlate2d(b, a, mode='same')
    y0, x0 = np.unravel_index(np.argmax(c), c.shape)
    return np.array([x0 - c.shape[1]//2, y0 - c.shape[0]//2], dtype=int)

def verify_deshift_subpixel():
    print("\n=== VERIFY DESHIFT ON SYNTHETIC IMAGES ===")
    base = _mk_synthetic()
    # --- Synthetic base image (your Gaussian or any high-contrast object) ---
    H = W = 128
    Y, X = np.mgrid[0:H, 0:W]
    base = np.exp(-((X - 64) ** 2 + (Y - 64) ** 2) / (2 * 3.0 ** 2))
    base += ((X > 55) & (X < 60) & (Y > 55) & (Y < 60)).astype(float)
    base = base.astype(np.float64)
    # --- True shifts: (dx, dy) in COLS, ROWS ---
    true_shifts = [
        np.array([0, 0]),
        np.array([-3.10, -0.90]),
        np.array([0.60, 3.20]),
        np.array([-1.80, 2.50]),
        np.array([2.4, -1.75]),
    ]
    # --- Generate shifted copies using YOUR imtranslate_xy (must use shift=(dy,dx)) ---
    imgs = [imtranslate_xy(base, dx=r[0], dy=r[1], order=3) for r in true_shifts]
    # --- Recover using your r_mn directly against img1 as reference ---
    rec_shifts = [np.array([0.0, 0.0])]
    for k in range(1, len(imgs)):
        rec = r_mn(imgs[0], imgs[k], do_filt=0)  # <-- IMPORTANT: no filtering for synthetic test
        rec_shifts.append(rec)
    # --- Print true vs recovered ---
    print("\nTrue shifts (dx,dy):")
    for i, r in enumerate(true_shifts, 1):
        print(f"  img{i}: {r}")
    print("\nRecovered shifts (dx,dy) from r_mn:")
    for i, r in enumerate(rec_shifts, 1):
        print(f"  img{i}: {r}")
    # --- Align and compute RMS differences ---
    aligned = [imtranslate_xy(imgs[k], dx=-rec_shifts[k][0], dy=-rec_shifts[k][1])
               for k in range(len(imgs))]
    rms = [np.sqrt(np.mean((aligned[0] - aligned[k]) ** 2)) for k in range(len(imgs))]
    print("\nRMS after alignment:")
    print(rms)
    print("=== END VERIFY DESHIFT ===\n")


# ============================================================
# piDPC reconstruction (ORIGINAL mapping only; quad_rot optional)
# ============================================================
def pidpc_from_sects(out_base: str,
                     pixel_size_Axy: tuple[float, float] | None = None,
                     flip_kx: bool = False, flip_ky: bool = False,
                     quad_rot: int = 0) -> str:
    """
    MATLAB-parity piDPC reconstruction from 8 sector images (inside BF).
    Matches the  MATLAB reference:
        [Y, X] = meshgrid((1:nY)-(1+nY)/2, (1:nX)-(1+nX)/2);
        kyp = Y/nY;  kxp = X/nX;
        kpnorm2 = kxp.^2 + kyp.^2;  kpnorm2(kpnorm2==0) = 1e-6;
        for tiltno=1:ntilts
            [corr_offset(1,:),...,corr_offset(4,:)] = deshift(img1,img2,img3,img4);
            ... de-shift → Gx,Gy ...
            iDPC1fft = (1/(1i*2*pi))*((kxp.*F(Gx)+kyp.*F(Gy)).*(1-1*(abs(kpnorm2)<1e-9)))./kpnorm2;
            iDPC1 = real(ifftshift(ifft2(fftshift(iDPC1fft))));
            ... repeat with ±1 px dither ...
            DiDPC1_BP = (iDPC11 - iDPC1) - imgaussfilt(iDPC11 - iDPC1, 50);
            iDPC11tilt(:,:,tiltno) = DiDPC1_BP;
        end
    """

    # --- Load 8 corrected sectors (inside) into channels_by_sector[0..7] ---
    channels_by_sector = None
    sizeX_A = sizeY_A = None
    nY = nX = ntilts = None
    for s in range(1, 9):
        p = f"{out_base}_sect{s}.mrc"
        with mrcfile.open(p, permissive=True) as m:
            vol = m.data
            if vol.ndim == 2:
                tilt = vol.astype(np.float32)
                nY, nX = tilt.shape
                tilt = tilt[..., np.newaxis]
                ntilts_local = 1
            elif vol.ndim == 3:
                tilt = np.moveaxis(vol.astype(np.float32), 0, -1)  # (Y, X, Z)
                nY, nX, ntilts_local = tilt.shape
            else:
                raise ValueError("Unsupported MRC dimensionality.")
            if channels_by_sector is None:
                ntilts = ntilts_local
                channels_by_sector = [np.zeros((nY, nX, ntilts), dtype=np.float32) for _ in range(8)]
            else:
                if ntilts_local == 1 and ntilts > 1:
                    tilt = np.repeat(tilt, ntilts, axis=2)
            channels_by_sector[s - 1] += tilt

            try:
                sizeX_A = float(m.voxel_size.x); sizeY_A = float(m.voxel_size.y)
            except Exception:
                pass

    # --- Build quadrants from sectors: MATLAB mapping [2,2,3,3,4,4,1,1] then rotate by quad_rot ---
    sect2quad = [2, 2, 3, 3, 4, 4, 1, 1]                             # sectors 1..8 → quadrants 1..4

    tilts_channels = np.zeros((nY, nX, ntilts, 4), dtype=np.float32)
    for s in range(8):
        qidx = sect2quad[s] - 1                                      # 0..3
        tilts_channels[:, :, :, qidx] += channels_by_sector[s]

    # --- MATLAB k-grid parity ---
    # [Y, X] = meshgrid((1:nY)-(1+nY)/2, (1:nX)-(1+nX)/2)  (note the order)
    X, Y = np.meshgrid(
        (np.arange(1, nY + 1) - (1 + nY) / 2.0),
        (np.arange(1, nX + 1) - (1 + nX) / 2.0),
        indexing='xy'
    )
    # kxp/kyp as in MATLAB:
    kxp = X / float(nX)
    kyp = Y / float(nY)
    if flip_kx: kxp = -kxp
    if flip_ky: kyp = -kyp

    kpnorm2 = kxp**2 + kyp**2
    kpnorm2 = kpnorm2.astype(np.float64, copy=False)
    kpnorm2[kpnorm2 == 0] = 1e-6
    # MATLAB mask (1 - 1*(abs(kpnorm2)<1e-9)) → zero DC exactly
    dc_mask = (np.abs(kpnorm2) >= 1e-9).astype(np.float64)

    def F(img):
         return spfft.fftshift(spfft.fft2(img))

    def Finv(spec):
        return np.real(spfft.ifft2(spfft.ifftshift(spec)))

    # Pre-allocate as MATLAB: iDPC11tilt is stored as (nX, nY, ntilts)
    iDPC11tilt_store = np.zeros((nX, nY, ntilts), dtype=np.float32)

    for t in range(ntilts):
        img1 = tilts_channels[:, :, t, 0]
        img2 = tilts_channels[:, :, t, 1]
        img3 = tilts_channels[:, :, t, 2]
        img4 = tilts_channels[:, :, t, 3]

        # deshift → corr_offset(1,:)..(4,:) (dx,dy) each
        r1, r2, r3, r4 = deshift(img1, img2, img3, img4)  # each r = [dx, dy]

        # (Optional) MATLAB's shift_avg_pix (not used afterwards)
        shift_avg_pix = (-r1[0] - r1[1] - r2[0] + r2[1] + r3[0] + r3[1] + r4[0] - r4[1]) / 8.0
        _ = shift_avg_pix  # keep for parity; not used later

        # ---- Pass 1: exact MATLAB de-shift and integration ----
        i1 = imtranslate_xy(img1, dx=-r1[0], dy=-r1[1])
        i2 = imtranslate_xy(img2, dx=-r2[0], dy=-r2[1])
        i3 = imtranslate_xy(img3, dx=-r3[0], dy=-r3[1])
        i4 = imtranslate_xy(img4, dx=-r4[0], dy=-r4[1])

        img_gradx = -(i1 + i2 - i3 - i4)
        img_grady = -(i2 - i1 + i3 - i4)

        GxF = F(img_gradx)
        GyF = F(img_grady)

        numer = (kxp * GxF + kyp * GyF) * dc_mask
        iDPC1fft = (1.0 / (1j * 2.0 * np.pi)) * (numer / kpnorm2)
        iDPC1 = Finv(iDPC1fft)

        # ---- Pass 2: ±1 px dither (same offsets as MATLAB) ----
        i1b = imtranslate_xy(img1, dx=-r1[0] + 1, dy=-r1[1] - 1)  # [ +1, -1 ]
        i2b = imtranslate_xy(img2, dx=-r2[0] - 1, dy=-r2[1] - 1)  # [ -1, -1 ]
        i3b = imtranslate_xy(img3, dx=-r3[0] - 1, dy=-r3[1] + 1)  # [ -1, +1 ]
        i4b = imtranslate_xy(img4, dx=-r4[0] + 1, dy=-r4[1] + 1)  # [ +1, +1 ]

        img_gradx2 = -(i1b + i2b - i3b - i4b)
        img_grady2 = -(i2b - i1b + i3b - i4b)

        GxF2 = F(img_gradx2)
        GyF2 = F(img_grady2)

        numer2 = (kxp * GxF2 + kyp * GyF2) * dc_mask
        iDPC11fft = (1.0 / (1j * 2.0 * np.pi)) * (numer2 / kpnorm2)
        iDPC11 = Finv(iDPC11fft)

        # Band-pass: (iDPC11 - iDPC1) minus its low-pass (σ=50)
        DiDPC1 = iDPC11 - iDPC1
        DiDPC1_LP = gaussian_filter(DiDPC1, sigma=50)
        DiDPC1_BP = DiDPC1 - DiDPC1_LP

        # Store as (nX, nY, ntilts)
        iDPC11tilt_store[:, :, t] = DiDPC1_BP.astype(np.float32)

    # --- Write MRC stack ---
    out_pidpc = f"{out_base}_piDPC.mrc"
    with mrcfile.new(out_pidpc, overwrite=True) as m:
        # MRC expects (Z, Y, X): we already have (nX, nY, ntilts) → transpose to (ntilts, nY, nX)
        m.set_data(np.moveaxis(iDPC11tilt_store, -1, 0))
        if pixel_size_Axy is not None:
            m.voxel_size = (float(pixel_size_Axy[0]), float(pixel_size_Axy[1]), 1.0)
        elif sizeX_A is not None and sizeY_A is not None:
            m.voxel_size = (float(sizeX_A), float(sizeY_A), 1.0)

    return out_pidpc

# ============================================================
# Driver
# ============================================================
def run_pidpc(h5_path: str, out_base: str,
              pixel_size_Axy: tuple[float, float] | None = None,
              use_raw: bool = False,
              nx_hint: int | None = None,
              ny_hint: int | None = None,
              flip_kx: bool = False, flip_ky: bool = False,
              sector_rot: int = 0,
              theta_mirror: bool = False,
              sector_reverse: bool = False,
              digest: bool = False,
              quad_rot: int = 0) -> str:

    #verify_deshift_subpixel()
    if use_raw:
        mat, Nx, Ny = load_4dstem_h5_to_mat_raw(h5_path)
    else:
        mat, Nx, Ny = load_legacy_h5_to_mat(h5_path, nx_hint=nx_hint, ny_hint=ny_hint)
        # >>> One-time XY swap of detector plane <<<
        # mat shape is (H, W, T). We want MATLAB-compatible orientation.
        mat = np.transpose(mat, (1, 0, 2))  # (W, H, T) → now columns/x match qY and rows/y match qX

    arina_tomo_rings_v9(mat, out_base=out_base, nX=Nx, nY=Ny,
                        flip_kx=flip_kx, flip_ky=flip_ky,
                        sector_rot=sector_rot,
                        theta_mirror=theta_mirror,
                        digest=digest)

    pidpc_path = pidpc_from_sects(out_base=out_base,
                                  pixel_size_Axy=pixel_size_Axy,
                                  flip_kx=flip_kx, flip_ky=flip_ky,
                                  quad_rot=quad_rot)
    return pidpc_path


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="piDPC pipeline: H5 → 16 sectors → piDPC (helpers-integrated, no theta modulo)")

    parser.add_argument("h5_path", type=str, help="Path to input 4D-STEM HDF5")
    parser.add_argument("--out-base", type=str, default=None, help="Output base path (prefix)")
    parser.add_argument("--px", type=float, default=None, help="Pixel size X in Å (optional)")
    parser.add_argument("--py", type=float, default=None, help="Pixel size Y in Å (optional)")

    # H5 mode / scan hints
    parser.add_argument("--raw", action="store_true", help="Use RAW reader (top-level row groups). Default: legacy /entry/data")
    parser.add_argument("--nx", type=int, default=None, help="Scan width override (optional)")
    parser.add_argument("--ny", type=int, default=None, help="Scan height override (optional)")

    # Axis flips & wedge rotation
    parser.add_argument("--flip-kx", action="store_true", help="Flip sign of kx axis (applied to rings & DPC)")
    parser.add_argument("--flip-ky", action="store_true", help="Flip sign of ky axis (applied to rings & DPC)")
    parser.add_argument("--sector-rot", type=int, default=0, help="Rotate 8-sector wedges by N⋅45° (N=0..7)")
    parser.add_argument("--theta-mirror", action="store_true", help="Reverse angular direction: theta := (2π - theta) (no modulo)")
    # Diagnostics: write 16 wedge masks (1..8 inside, 9..16 outside) as MRCs for verification
    parser.add_argument("--digest", action="store_true", help="Write average in sectors and images of the 16 sector masks as MRCs.")
    # Quadrant rotation for DPC (0,1,2,3 → 0/90/180/270°)
    parser.add_argument("--quad-rot", type=int, default=0, help="Rotate sector→quadrant mapping by Q⋅90° (Q=0..3)")

    args = parser.parse_args()
    base = args.out_base or os.path.splitext(args.h5_path)[0]
    px_xy = (args.px, args.py) if (args.px is not None and args.py is not None) else None

    out = run_pidpc(args.h5_path, base,
                    pixel_size_Axy=px_xy,
                    use_raw=args.raw,
                    nx_hint=args.nx, ny_hint=args.ny,
                    flip_kx=args.flip_kx, flip_ky=args.flip_ky,
                    sector_rot=args.sector_rot,
                    theta_mirror=args.theta_mirror,
                    digest=args.digest,
                    quad_rot=args.quad_rot)
    print("piDPC written to:", out)

