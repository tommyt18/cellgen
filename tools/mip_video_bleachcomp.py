"""Bleach-compensated XY MIP movies for the website (long CPE sequences).

Same layout as simcell/gen/mip_video.py (raw | clean, nuc magenta, cyto green), but the display
follows the simulated photobleaching.  The acquisition is  F_t = X_t * g_t + O  with g_t =
exp(-k t) (k = acquisition.bleach_per_frame from config.yaml), X the signal + autofluorescence and
O the static camera offset field.  Shown is (F_t - O) / (hi0 * g_t):
  * O is estimated from the clean movie's last frames by extrapolating the decay
    (O = F_T - (F_{T-D} - F_T) / (exp(k D) - 1), then smoothed; O is a smooth field);
  * hi0 is the 99.8 percentile of F - O over the first frames of the clean movie (both panels);
  * the raw panel also subtracts the per-frame median of raw MIP - clean MIP (the bias a
    max-projection picks up from noise), so its noise is centred on the clean panel.
With the fixed range of mip_video.py a 1440-frame sequence (k = 0.004 -> 0.3% signal at the end)
is black after ~10 h.  The data are untouched; the raw panel turns to noise as the signal bleaches.

  PYTHONPATH=<simcell> python mip_video_bleachcomp.py <seq_dir> ... --out <dir> [--stride 2]
"""
import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor

import imageio.v2 as imageio
import numpy as np
import yaml
import zarr
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter

from gen.mip_video import composite, pick_level


def mipf(arr, t):
    return np.asarray(arr[t]).max(axis=1).astype(np.float32)          # (c, y, x)


def offset_field(clean, k, D=300, n=5, sigma=6.0):
    T = clean.shape[0]
    late = np.mean([mipf(clean, t) for t in range(T - n, T)], axis=0)
    early = np.mean([mipf(clean, t) for t in range(T - n - D, T - D)], axis=0)
    O = late - (early - late) / (np.exp(k * D) - 1)
    return np.stack([gaussian_filter(o, sigma, mode="nearest") for o in O])


def render(seq, out, max_px, stride):
    meta = json.load(open(os.path.join(seq, "meta.json")))
    cfg = yaml.safe_load(open(os.path.join(seq, "config.yaml")))
    k = float(cfg["acquisition"]["bleach_per_frame"])
    chans = list(cfg["render"]["channels"])
    za = json.load(open(os.path.join(seq, "raw.ome.zarr", ".zattrs")))
    chans = [c["label"] for c in za.get("omero", {}).get("channels", [])] or chans
    name = os.path.basename(seq.rstrip("/"))
    src = {kk: zarr.open_group(os.path.join(seq, f"{kk}.ome.zarr"), mode="r") for kk in ("raw", "clean")}
    lvl = pick_level(src["raw"], max_px)
    raw, clean = src["raw"][lvl], src["clean"][lvl]
    T = raw.shape[0]
    O = offset_field(clean, k)
    hi0 = np.percentile(np.stack([mipf(clean, t) - O for t in range(6)]), 99.8, axis=(0, 2, 3))
    zero = np.zeros_like(hi0)
    dt = float(meta["frame_interval_min"])
    frames = range(0, T, stride)
    fps = int(np.clip(round(len(frames) / 24), 8, 30))
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, name + ".mp4")
    font = ImageFont.load_default()
    w = imageio.get_writer(path, fps=fps, codec="libx264", quality=None, macro_block_size=1,
                           ffmpeg_params=["-crf", "24", "-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    for t in frames:
        g = np.exp(-k * t)
        c = mipf(clean, t)
        r = mipf(raw, t)
        bias = np.median(r - c, axis=(1, 2), keepdims=True)
        panels = [composite(r - O - bias, zero, hi0 * g, chans), composite(c - O, zero, hi0 * g, chans)]
        gap = np.full((panels[0].shape[0], 4, 3), 40, np.uint8)
        frame = np.concatenate([panels[0], gap, panels[1]], axis=1)
        bar = np.zeros((30, frame.shape[1], 3), np.uint8)
        img = Image.fromarray(np.concatenate([bar, frame], axis=0))
        d = ImageDraw.Draw(img)
        d.text((4, 2), f"{name}   t={t}  {t * dt / 60:.1f} h   display gain x{1 / g:.0f}", fill=(230, 230, 230), font=font)
        d.text((4, 16), "raw", fill=(160, 160, 160), font=font)
        d.text((panels[0].shape[1] + 10, 16), "clean (noise-free)", fill=(160, 160, 160), font=font)
        a8 = np.asarray(img)
        w.append_data(a8[: a8.shape[0] // 2 * 2, : a8.shape[1] // 2 * 2])
    w.close()
    return f"ok {path} T={T} k={k} fps={fps} level={lvl}"


def _job(a):
    try:
        return render(*a)
    except Exception as e:  # noqa: BLE001
        return f"FAIL {a[0]} {type(e).__name__}: {e}"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seqs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max_px", type=int, default=640)
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--workers", type=int, default=9)
    a = ap.parse_args()
    with ProcessPoolExecutor(a.workers) as ex:
        for r in ex.map(_job, [(s, a.out, a.max_px, a.stride) for s in a.seqs]):
            print(r, flush=True)
