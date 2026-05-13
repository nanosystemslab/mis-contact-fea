#!/usr/bin/env python3
"""Load push + retention force curves, suppress single-point numerical
spikes, plot raw vs smoothed overlays per shape, and compute the
insertion/retention force ratio.

Outputs:
    force_smoothed_overlay.png   raw + smoothed for each shape (7x2 grid)
    force_smoothed_push.png      all smoothed push curves on one axis
    force_smoothed_retention.png all smoothed retention curves on one axis
    force_ratio.png              bar chart of |F_peak_push| / |F_peak_retention|
    force_peaks.csv              peak values and ratio per shape

Usage:
    python3 src/analyze_forces.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SHAPES = ["sphere", "oblate", "prolate", "cone", "cap", "capcone", "doublecone"]

# Colors matched to other plots in the project.
COLORS = {
    "sphere": "#1f77b4", "oblate": "#ff7f0e", "prolate": "#2ca02c",
    "cone": "#d62728", "cap": "#9467bd", "capcone": "#8c564b",
    "doublecone": "#e377c2",
}

PUSH_ROOT = Path("runs_all_shapes")
RETENTION_ROOT = Path("runs_retention")


def load_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = np.loadtxt(path, delimiter=",", skiprows=1)
    return rows[:, 0], rows[:, 1]


def hampel_filter(x: np.ndarray, half_window: int = 3, n_sigmas: float = 3.0) -> np.ndarray:
    """Replace points >n_sigmas * MAD from local median with that median.
    Pure numpy, no scipy. Robust to single- and multi-point spikes."""
    out = x.astype(float).copy()
    n = len(x)
    for i in range(n):
        lo, hi = max(0, i - half_window), min(n, i + half_window + 1)
        win = x[lo:hi]
        med = np.median(win)
        mad = 1.4826 * np.median(np.abs(win - med))
        if mad > 0 and abs(x[i] - med) > n_sigmas * mad:
            out[i] = med
    return out


def rolling_median(x: np.ndarray, half_window: int = 4) -> np.ndarray:
    """Robust to oscillating sign-flip chatter (median of alternating
    +/- values is near zero, washes out unphysical oscillations)."""
    n = len(x)
    out = np.empty_like(x, dtype=float)
    for i in range(n):
        lo, hi = max(0, i - half_window), min(n, i + half_window + 1)
        out[i] = np.median(x[lo:hi])
    return out


def moving_average(x: np.ndarray, half_window: int = 2) -> np.ndarray:
    """Centered moving average. Edges use shorter windows."""
    n = len(x)
    out = np.empty_like(x, dtype=float)
    for i in range(n):
        lo, hi = max(0, i - half_window), min(n, i + half_window + 1)
        out[i] = np.mean(x[lo:hi])
    return out


def detect_chatter_mask(force: np.ndarray, half_window: int = 5,
                        flips_threshold: int = 4) -> np.ndarray:
    """Return a boolean mask, True where the force curve has frequent
    sign-flips in a local window — a strong indicator of Newton-contact
    chatter rather than physical force. Peaks inside these regions are
    untrustworthy."""
    n = len(force)
    mask = np.zeros(n, dtype=bool)
    signs = np.sign(force)
    for i in range(n):
        lo, hi = max(0, i - half_window), min(n, i + half_window + 1)
        win = signs[lo:hi]
        flips = int(np.sum(np.abs(np.diff(win)) > 0))
        if flips >= flips_threshold:
            mask[i] = True
    return mask


def smooth_force(disp: np.ndarray, force: np.ndarray) -> np.ndarray:
    """Three-stage cleanup tuned for the contact-chatter we observed:
        1. rolling median (window=9) kills sign-flip oscillations
        2. hampel (window=7, t=2.5) catches any residual spikes
        3. short moving average for cosmetic smoothness
    The rolling median in step 1 is the key change — pure hampel can't
    handle multi-point chatter because its local median is also
    oscillating in chatter regions."""
    stage1 = rolling_median(force, half_window=4)        # 9-pt median
    stage2 = hampel_filter(stage1, half_window=3, n_sigmas=2.5)
    stage3 = moving_average(stage2, half_window=2)        # 5-pt mean
    return stage3


def safe_load(root: Path, shape: str) -> tuple[np.ndarray, np.ndarray] | None:
    p = root / shape / "results" / "force_displacement.csv"
    if not p.is_file():
        return None
    return load_csv(p)


def peak(disp: np.ndarray, force: np.ndarray) -> tuple[float, float]:
    """Return (|F|_peak, disp_at_peak) on the smoothed curve."""
    af = np.abs(force)
    i = int(np.argmax(af))
    return float(af[i]), float(disp[i])


def peak_excluding_chatter(disp: np.ndarray, force_smoothed: np.ndarray,
                           force_raw: np.ndarray) -> tuple[float, float, bool]:
    """Pick the peak from the smoothed curve but only in regions where
    the RAW curve isn't chattering. Returns (|F|_peak, disp_at_peak,
    in_chatter_region_flag) — flag True if the peak of the smoothed
    curve sits inside a chatter region (low confidence)."""
    chatter = detect_chatter_mask(force_raw, half_window=5, flips_threshold=4)
    af = np.abs(force_smoothed)
    masked = af.copy()
    masked[chatter] = 0.0
    if masked.max() == 0.0:
        # entire curve is chatter — fall back to the smoothed peak
        i = int(np.argmax(af))
        return float(af[i]), float(disp[i]), True
    i = int(np.argmax(masked))
    raw_i = int(np.argmax(af))
    return float(af[i]), float(disp[i]), bool(chatter[raw_i])


def build_overlay(out_path: Path) -> None:
    """Per-shape raw + smoothed plot for both push and retention."""
    n = len(SHAPES)
    fig, axes = plt.subplots(n, 2, figsize=(13, 2.4 * n), dpi=150, squeeze=False)
    for row, shape in enumerate(SHAPES):
        for col, (label, root) in enumerate([("push", PUSH_ROOT), ("retention", RETENTION_ROOT)]):
            ax = axes[row][col]
            d_f = safe_load(root, shape)
            if d_f is None:
                ax.set_title(f"{shape} — {label}: (no data)")
                continue
            d, f = d_f
            f_s = smooth_force(d, f)
            chatter = detect_chatter_mask(f, half_window=5, flips_threshold=4)
            # Shade chatter region
            if chatter.any():
                # Find contiguous chatter runs and shade them
                in_run = False
                for k, c in enumerate(chatter):
                    if c and not in_run:
                        run_start = k; in_run = True
                    elif not c and in_run:
                        ax.axvspan(d[run_start], d[k - 1], color="red", alpha=0.12, lw=0)
                        in_run = False
                if in_run:
                    ax.axvspan(d[run_start], d[-1], color="red", alpha=0.12, lw=0)
            ax.plot(d, np.abs(f), color="lightgray", linewidth=0.8, label="raw |F|")
            ax.plot(d, np.abs(f_s), color=COLORS[shape], linewidth=1.6, label="smoothed |F|")
            pk_clean, pk_clean_at, in_chatter = peak_excluding_chatter(d, f_s, f)
            ax.axvline(pk_clean_at, color="black", linestyle=":", linewidth=0.8)
            tag = "  (chatter!)" if in_chatter else ""
            ax.set_title(
                f"{shape} — {label}   clean peak {pk_clean:.0f} µN @ disp={pk_clean_at:.1f}{tag}"
            )
            ax.grid(alpha=0.2)
            ax.set_xlabel("displacement (µm)")
            ax.set_ylabel("|F_z| (µN/µm-depth)")
            ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def build_combined_plots(out_push: Path, out_retention: Path) -> None:
    """All shapes overlaid, smoothed |F| vs displacement."""
    for label, root, out in [("push", PUSH_ROOT, out_push), ("retention", RETENTION_ROOT, out_retention)]:
        fig, ax = plt.subplots(figsize=(11, 6), dpi=150)
        for shape in SHAPES:
            d_f = safe_load(root, shape)
            if d_f is None:
                continue
            d, f = d_f
            f_s = smooth_force(d, f)
            ax.plot(d, np.abs(f_s), color=COLORS[shape], label=shape, linewidth=1.8)
        ax.set_xlabel("displacement (µm)")
        ax.set_ylabel("|F_z| (µN per µm-depth)")
        ax.set_title(f"Smoothed {label} force vs displacement")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=10, loc="upper right")
        fig.tight_layout()
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved {out}")


def build_ratio_plot_and_csv(plot_path: Path, csv_path: Path) -> list[dict]:
    rows = []
    for shape in SHAPES:
        push = safe_load(PUSH_ROOT, shape)
        retn = safe_load(RETENTION_ROOT, shape)
        if push is None or retn is None:
            print(f"  WARN: {shape} missing push or retention")
            continue
        dp, fp = push
        dr, fr = retn
        fp_s = smooth_force(dp, fp)
        fr_s = smooth_force(dr, fr)
        # Primary peak: smoothed-curve max (rolling median handles chatter OK).
        push_peak_full, push_at_full = peak(dp, fp_s)
        retn_peak_full, retn_at_full = peak(dr, fr_s)
        # Diagnostic: peak excluding chatter region (more conservative).
        push_peak_clean, push_at_clean, push_chatter = peak_excluding_chatter(dp, fp_s, fp)
        retn_peak_clean, retn_at_clean, retn_chatter = peak_excluding_chatter(dr, fr_s, fr)
        ratio = retn_peak_full / push_peak_full if push_peak_full > 0 else float("nan")
        rows.append(dict(
            shape=shape,
            push_peak_uN_per_um_depth=round(push_peak_full, 3),
            push_at_disp_um=round(push_at_full, 3),
            push_in_chatter=push_chatter,
            push_peak_clean_uN=round(push_peak_clean, 3),
            push_at_clean_disp_um=round(push_at_clean, 3),
            retention_peak_uN_per_um_depth=round(retn_peak_full, 3),
            retention_at_disp_um=round(retn_at_full, 3),
            retention_in_chatter=retn_chatter,
            retention_peak_clean_uN=round(retn_peak_clean, 3),
            retention_at_clean_disp_um=round(retn_at_clean, 3),
            retention_insertion_ratio=round(ratio, 3),
        ))

    # CSV
    with csv_path.open("w", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Saved {csv_path}")

    # Bar chart of ratio (sorted high-to-low so best interlocks land left).
    sorted_rows = sorted(rows, key=lambda r: r["retention_insertion_ratio"], reverse=True)
    fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
    labels = [r["shape"] for r in sorted_rows]
    ratios = [r["retention_insertion_ratio"] for r in sorted_rows]
    colors = [COLORS[s] for s in labels]
    bars = ax.bar(labels, ratios, color=colors, edgecolor="black", linewidth=0.6)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8,
               label="ratio = 1 (retention = insertion)")
    for b, r in zip(bars, ratios):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                f"{r:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("peak retention / peak insertion")
    ax.set_title("Retention-to-insertion force ratio (larger = better interlock)")
    ax.grid(alpha=0.2, axis="y")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {plot_path}")

    return rows


def main():
    print("Building per-shape raw + smoothed overlays ...")
    build_overlay(Path("force_smoothed_overlay.png"))
    print("Building all-shapes smoothed comparison plots ...")
    build_combined_plots(Path("force_smoothed_push.png"),
                         Path("force_smoothed_retention.png"))
    print("Computing peaks and ratios ...")
    rows = build_ratio_plot_and_csv(Path("force_ratio.png"),
                                    Path("force_peaks.csv"))

    print()
    print("Peak summary — primary value is the smoothed-curve peak (rolling-median")
    print("handles chatter); CHATTER flag means the peak was in a chatter region.")
    print("R/I ratio: retention/insertion — larger = better interlock.")
    sorted_for_print = sorted(rows, key=lambda r: r["retention_insertion_ratio"], reverse=True)
    print(f"  {'shape':<11} {'push_peak':>11} {'@disp':>6} {'flag':>8} "
          f"{'ret_peak':>10} {'@disp':>6} {'flag':>8} {'R/I':>7}")
    for r in sorted_for_print:
        pq = "CHATTER" if r["push_in_chatter"] else "ok"
        rq = "CHATTER" if r["retention_in_chatter"] else "ok"
        print(f"  {r['shape']:<11} "
              f"{r['push_peak_uN_per_um_depth']:>11.2f} "
              f"{r['push_at_disp_um']:>6.1f} {pq:>8} "
              f"{r['retention_peak_uN_per_um_depth']:>10.2f} "
              f"{r['retention_at_disp_um']:>6.1f} {rq:>8} "
              f"{r['retention_insertion_ratio']:>7.2f}")
    print()
    print("Diagnostic — peak EXCLUDING chatter regions (lower bound):")
    print(f"  {'shape':<11} {'push_clean':>11} {'@disp':>6}   {'ret_clean':>10} {'@disp':>6}")
    for r in rows:
        print(f"  {r['shape']:<11} "
              f"{r['push_peak_clean_uN']:>11.2f} "
              f"{r['push_at_clean_disp_um']:>6.1f}   "
              f"{r['retention_peak_clean_uN']:>10.2f} "
              f"{r['retention_at_clean_disp_um']:>6.1f}")


if __name__ == "__main__":
    main()
