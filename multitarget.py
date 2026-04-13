from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

import singleframe as sf


def detect_multiple_targets(
    ra_power: np.ndarray,
    min_r: int = sf.MIN_RANGE_BIN,
    max_r: int = sf.MAX_RANGE_BIN,
    max_angle_deg: float = 60.0,
    max_peaks_per_range: int = 2,
    max_peaks: int = 4,
    db_below_global_max: float = 15.0,
    min_sep_angle_bins: int = 6,
    min_sep_range_bins: int = 4,
) -> list[tuple[float, float, float]]:
    """
    Multi-target detection with per-range angular peaks, then global 2D deduplication.

    For each range bin, up to ``max_peaks_per_range`` peaks are found along angle only. 
    That way several targets at the **same range** but different angles can all appear 
    without one 2D box wiping the other in range.

    Candidates are merged globally: sorted by power, keep a peak unless a stronger
    kept peak is within min_sep_angle_bins and min_sep_range_bins (same idea
    as 2D NMS: suppress duplicate/sidelobe hits that are close in both dimensions).

    Returns list of (range_m, angle_deg_reported, peak_linear_power), strongest first.
    """
    gated = np.asarray(ra_power, dtype=np.float64).copy()
    gated[:, :min_r] = -np.inf
    gated[:, max_r + 1 :] = -np.inf
    angle_mask = np.abs(sf.ANGLES_DEG) > max_angle_deg
    gated[angle_mask, :] = -np.inf

    valid = np.isfinite(gated)
    if not np.any(valid):
        return []

    pmax = float(np.max(gated[valid]))
    thresh = pmax * 10 ** (-db_below_global_max / 10.0)

    na, _nr = gated.shape
    da = max(1, min_sep_angle_bins // 2)

    # Stage 1: per range bin, greedy peaks along angle only 
    candidates: list[tuple[float, int, int]] = []
    for r_idx in range(min_r, max_r + 1):
        row = gated[:, r_idx].copy()
        for _ in range(max_peaks_per_range):
            p = float(np.nanmax(row))
            if not np.isfinite(p) or p < thresh:
                break
            a_idx = int(np.nanargmax(row))
            candidates.append((p, a_idx, r_idx))
            a0 = max(0, a_idx - da)
            a1 = min(na, a_idx + da + 1)
            row[a0:a1] = -np.inf

    if not candidates:
        return []

    # Strongest first for NMS
    candidates.sort(key=lambda t: -t[0])

    # Stage 2: global 2D NMS 
    kept: list[tuple[float, int, int]] = []
    for p, a_idx, r_idx in candidates:
        if len(kept) >= max_peaks:
            break
        too_close = False
        for _pk, ak, rk in kept:
            if (
                abs(a_idx - ak) < min_sep_angle_bins
                and abs(r_idx - rk) < min_sep_range_bins
            ):
                too_close = True
                break
        if not too_close:
            kept.append((p, a_idx, r_idx))

    out: list[tuple[float, float, float]] = []
    for p, a_idx, r_idx in kept:
        rng_m = float(r_idx * sf.RANGE_RES)
        ang_rep = float(-sf.ANGLES_DEG[a_idx])
        out.append((rng_m, ang_rep, p))

    return out


def main() -> None:
    print(
        f"Config: {sf.CONFIG_XML.name}  |  slope={sf.FREQ_SLOPE_MHZ_US} MHz/us  |  "
        f"f_s={sf.SAMPLE_RATE/1e6:.3f} MS/s  |  N_ADC={sf.NUM_ADC}  |  loops={sf.NUM_LOOPS}  |  "
        f"R_res={sf.RANGE_RES*100:.2f} cm/bin  |  R_max={sf.MAX_RANGE:.2f} m"
    )

    blob = sf.ADC_BIN.read_bytes()
    if len(blob) < sf.BYTES_PER_FRAME:
        raise SystemExit(
            f"{sf.ADC_BIN} is too small for one frame "
            f"(need {sf.BYTES_PER_FRAME} bytes, got {len(blob)})"
        )
    if len(blob) % sf.BYTES_PER_FRAME != 0:
        print(
            f"Warning: trailing {len(blob) % sf.BYTES_PER_FRAME} bytes ignored "
            f"(not a full frame)"
        )

    sf.print_diagnostics(blob)

    ra_db, _range_axis, ra_power, n_frames = sf.range_angle_map_multiframe(blob)
    dets = detect_multiple_targets(ra_power)

    print(f"Detections ({len(dets)} targets, reported angle convention):")
    for i, (rng_m, ang_deg, p) in enumerate(dets, start=1):
        print(f"  {i}: range = {rng_m:.3f} m,  angle = {ang_deg:+.2f} deg  (power={p:.3e})")

    ra_db_norm = ra_db - np.max(ra_db)
    display = ra_db_norm.T[:, ::-1]

    fig, ax = plt.subplots(figsize=(12, 7))
    im = ax.imshow(
        display,
        aspect="auto",
        origin="lower",
        cmap="inferno",
        extent=[
            float(sf.ANGLES_DEG[0]),
            float(sf.ANGLES_DEG[-1]),
            sf.DISPLAY_RANGE_MIN_M,
            sf.MAX_RANGE,
        ],
        interpolation="bilinear",
        vmin=-60,
        vmax=0.0,
    )
    ax.set_xlabel("Angle (deg)")
    ax.set_ylabel("Range (m)")

    if dets:
        xs = [d[1] for d in dets]
        ys = [d[0] for d in dets]
        ax.scatter(
            xs,
            ys,
            c="cyan",
            s=120,
            marker="+",
            linewidths=2.0,
            zorder=5,
            label=f"{len(dets)} targets",
        )
        for i, (rng_m, ang_deg, _p) in enumerate(dets, start=1):
            ax.annotate(
                str(i),
                (ang_deg, rng_m),
                textcoords="offset points",
                xytext=(6, 6),
                fontsize=10,
                color="white",
                fontweight="bold",
            )

    ax.legend(loc="upper right", fontsize=9)
    title_extra = (
        "  |  ".join(f"{i}: {d[0]:.2f}m @ {d[1]:+.0f}°" for i, d in enumerate(dets, start=1))
        if dets
        else "no detections"
    )
    ax.set_title(
        f"Multi-target range–angle  |  {sf.NUM_TX}x{sf.NUM_RX} MIMO  |  {n_frames} frames  |  "
        f"{sf.ADC_BIN.name}\n{title_extra}"
    )
    ax.set_ylim(sf.DISPLAY_RANGE_MIN_M, sf.DISPLAY_RANGE_MAX_M)
    plt.colorbar(
        im,
        ax=ax,
        fraction=0.03,
        pad=0.04,
        label="Relative power (dB re peak)",
    )
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
