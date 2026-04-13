"""
Live DCA1000 capture with snapshot-on-demand SAR aperture stepping.

Streams UDP from the DCA1000.
Each snapshot is labeled by along-track index in units of λ/2. UI matches sar.py:
live heatmap + keys b / ←/→ / c / u / f / q (plus d = diagnostics). Per
capture it averages the last N frames, runs singleframe / multitarget, stores
virtual-array slices for SAR, optional background subtraction, and f forms
the composite image.
"""

from __future__ import annotations

import datetime
import queue
import struct
import socket
import threading
import time
from collections import deque
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

import multitarget as mt
import singleframe as sf

# ── Network config ───────────────────────────────────────────
PC_IP = "192.168.33.30"
DATA_PORT = 4098

# ── Capture config ───────────────────────────────────────────
FRAMES_PER_SNAPSHOT = 8
RING_MAX_FRAMES = 128

# ── SAR imaging mode ─────────────────────────────────────────
# 'fft' : fast planar-wave FFT beamforming 
# 'bp'  : pixel backprojection 
SAR_MODE = "bp"

# ── Phase autofocus ──────────────────────────────────────────
# Set to the range (metres) of your static reference reflector.

AUTOFOCUS_REF_RANGE_M: float | None = 3.1

# ── Output dirs ──────────────────────────────────────────────
SNAPSHOT_SAVE_DIR: Path | None = None
SAR_OUTPUT_DIR: Path | None = None

# ── Angle grid ───────────────────────────────────────────────
N_SAR_ANGLES = 360


try:
    from numba import njit, prange
    _NUMBA = True
except ImportError:
    _NUMBA = False


# ── Backprojection core ──────────────────────────────────────

if _NUMBA:
    @njit(parallel=True, cache=True)
    def _bp_core(
        all_positions: np.ndarray,
        all_profiles_r: np.ndarray,
        all_profiles_i: np.ndarray,
        range_axis: np.ndarray,
        sar_angles_rad: np.ndarray,
        wavelength: float,
        range_res: float,
        num_range_bins: int,
    ) -> np.ndarray:
        n_ang = sar_angles_rad.shape[0]
        n_rng = num_range_bins
        n_elem = all_positions.shape[0]
        image_r = np.zeros((n_ang, n_rng), dtype=np.float64)
        image_i = np.zeros((n_ang, n_rng), dtype=np.float64)
        for ai in prange(n_ang):
            sin_t = np.sin(sar_angles_rad[ai])
            cos_t = np.cos(sar_angles_rad[ai])
            for ri in range(n_rng):
                R_n = range_axis[ri]
                if R_n < 0.05:
                    continue
                x_pix = R_n * sin_t
                y_pix = R_n * cos_t
                acc_r = 0.0
                acc_i = 0.0
                for ei in range(n_elem):
                    x_i = all_positions[ei]
                    dx = x_pix - x_i
                    R_slant = np.sqrt(dx * dx + y_pix * y_pix)
                    bin_f = R_slant / range_res
                    bin_lo = int(bin_f)
                    if bin_lo < 0 or bin_lo >= n_rng - 1:
                        continue
                    frac = bin_f - bin_lo
                    vr = ((1.0 - frac) * all_profiles_r[ei, bin_lo]
                          + frac * all_profiles_r[ei, bin_lo + 1])
                    vi = ((1.0 - frac) * all_profiles_i[ei, bin_lo]
                          + frac * all_profiles_i[ei, bin_lo + 1])
                    phase = 4.0 * np.pi / wavelength * R_slant
                    cp = np.cos(phase)
                    sp = np.sin(phase)
                    acc_r += vr * cp - vi * sp
                    acc_i += vr * sp + vi * cp
                image_r[ai, ri] = acc_r
                image_i[ai, ri] = acc_i
        return image_r * image_r + image_i * image_i

else:
    def _bp_core(
        all_positions, all_profiles_r, all_profiles_i,
        range_axis, sar_angles_rad, wavelength, range_res, num_range_bins,
    ):
        n_ang = len(sar_angles_rad)
        n_rng = num_range_bins
        image = np.zeros((n_ang, n_rng), dtype=np.complex128)
        all_profiles = all_profiles_r + 1j * all_profiles_i
        for ai, theta in enumerate(sar_angles_rad):
            sin_t = np.sin(theta)
            cos_t = np.cos(theta)
            x_pix = range_axis * sin_t
            y_pix = range_axis * cos_t
            for ei, x_i in enumerate(all_positions):
                dx = x_pix - x_i
                R_slant = np.sqrt(dx**2 + y_pix**2)
                bin_f = R_slant / range_res
                bin_lo = np.floor(bin_f).astype(np.int32)
                frac = bin_f - bin_lo
                valid = (bin_lo >= 0) & (bin_lo < n_rng - 1)
                bl = np.where(valid, bin_lo, 0)
                fr = np.where(valid, frac, 0.0)
                prof = all_profiles[ei]
                val = ((1.0 - fr) * prof[bl]
                       + fr * prof[np.minimum(bl + 1, n_rng - 1)])
                val = np.where(valid, val, 0.0)
                phase = 4.0 * np.pi / wavelength * R_slant
                image[ai] += val * np.exp(1j * phase)
        return np.abs(image)**2


def backproject_sar(
    sar_captures: list[np.ndarray],
    sar_offsets_m: list[float],
    sar_angles_rad: np.ndarray,
) -> np.ndarray:
    k = len(sar_captures)
    offsets = np.asarray(sar_offsets_m, dtype=np.float64)
    pos_local = sf.virtual_array_positions_m()
    n_virt = len(pos_local)
    range_axis = np.arange(sf.NUM_RANGE_BINS, dtype=np.float64) * sf.RANGE_RES

    all_positions = np.empty(k * n_virt, dtype=np.float64)
    all_profiles = np.empty((k * n_virt, sf.NUM_RANGE_BINS), dtype=np.complex64)
    for i in range(k):
        sl = slice(i * n_virt, (i + 1) * n_virt)
        all_positions[sl] = offsets[i] + pos_local
        all_profiles[sl, :] = sar_captures[i]

    sort_idx = np.argsort(all_positions)
    all_positions = all_positions[sort_idx]
    all_profiles = all_profiles[sort_idx, :]

    # Mirror aperture so angle axis matches physical convention
    all_positions = -all_positions

    win = np.hamming(len(all_positions)).astype(np.float32)
    all_profiles = all_profiles * win[:, np.newaxis]

    print(
        f"[BP] Backprojecting {len(all_positions)} phase centres × "
        f"{sf.NUM_RANGE_BINS} range bins × {len(sar_angles_rad)} angles"
        + (" [numba JIT]" if _NUMBA else " [numpy — may be slow]")
    )
    return _bp_core(
        all_positions,
        np.ascontiguousarray(all_profiles.real, dtype=np.float64),
        np.ascontiguousarray(all_profiles.imag, dtype=np.float64),
        range_axis,
        np.ascontiguousarray(sar_angles_rad, dtype=np.float64),
        float(sf.WAVELENGTH),
        float(sf.RANGE_RES),
        int(sf.NUM_RANGE_BINS),
    )


# ── Phase autofocus with unwrapped position refinement ───────

def autofocus_with_position_refinement(
    sar_captures: list[np.ndarray],
    sar_offsets_m: list[float],
) -> tuple[list[np.ndarray], list[float]]:
    ref_bin = int(round(AUTOFOCUS_REF_RANGE_M / sf.RANGE_RES))
    ref_bin = max(1, min(ref_bin, sf.NUM_RANGE_BINS - 2))
    print(f"[autofocus+pos] Reference bin {ref_bin} = "
          f"{ref_bin * sf.RANGE_RES:.3f} m")

    raw_phases = np.array([
        float(np.angle(cap[:, ref_bin].mean()))
        for cap in sar_captures
    ])

    nominal = np.array(sar_offsets_m)
    expected_slope = 4 * np.pi / sf.WAVELENGTH

    unwrapped = raw_phases.copy()
    for i in range(1, len(raw_phases)):
        expected_delta = expected_slope * (nominal[i] - nominal[0])
        raw_delta = raw_phases[i] - raw_phases[0]
        n_wraps = round((expected_delta - raw_delta) / (2 * np.pi))
        unwrapped[i] = raw_phases[i] + n_wraps * 2 * np.pi

    # Detect burst boundaries — large residual phase jumps that
    # exceed maximum position error (λ/4 = 90 deg)
    phase_jumps = np.abs(np.diff(unwrapped))
    burst_boundaries = np.where(phase_jumps > np.pi / 4)[0] + 1
    if len(burst_boundaries) > 0:
        print(f"[autofocus] Detected {len(burst_boundaries)} burst "
              f"boundaries at snapshots: {burst_boundaries.tolist()}")
        segments = np.split(np.arange(len(sar_captures)), burst_boundaries)
        ref_phase_per_segment = []
        for seg in segments:
            ref_phase_per_segment.append(unwrapped[seg[0]])
        for seg_idx, seg in enumerate(segments):
            offset = ref_phase_per_segment[seg_idx] - ref_phase_per_segment[0]
            unwrapped[seg] -= offset
    else:
        print("[autofocus] No burst boundaries detected — single session")

    focused = []
    print("[autofocus] Phase corrections:")
    for i, (cap, nominal_offset) in enumerate(
        zip(sar_captures, sar_offsets_m)
    ):
        delta = unwrapped[i] - unwrapped[0]
        correction = np.exp(-1j * delta).astype(np.complex64)
        focused.append(cap * correction)
        print(f"  snap {i:2d}: phase={np.rad2deg(raw_phases[i]):+7.1f}°  "
              f"unwrapped={np.rad2deg(unwrapped[i]):+8.1f}°  "
              f"correction={np.rad2deg(-delta):+7.1f}°")

    # Return nominal offsets — position refinement disabled
    return focused, list(sar_offsets_m)

# ── DC-corrected range profile ───────────────────────────────

def _virtual_range_profile_raw(blob: bytes, bg_reference: np.ndarray | None = None) -> np.ndarray:
    n_frames = len(blob) // sf.BYTES_PER_FRAME
    if n_frames == 0:
        raise ValueError("empty blob")
    win_range = np.blackman(sf.NUM_ADC).astype(np.float32)
    n_virt = len(sf.virtual_array_positions_m())
    acc = np.zeros((n_virt, sf.NUM_RANGE_BINS), dtype=np.complex128)
    for k in range(n_frames):
        iq = sf.raw_bytes_to_iq(
            blob[k * sf.BYTES_PER_FRAME : (k + 1) * sf.BYTES_PER_FRAME]
        )
        iq = iq - iq.mean(axis=2, keepdims=True)
        rfft = np.fft.fft(
            iq * win_range[np.newaxis, np.newaxis, :], axis=2
        )[:, :, : sf.NUM_RANGE_BINS]
        rfft[:, :, 0:2] = 0
        rfft_tx0 = rfft[0::2].mean(axis=0)
        rfft_tx1 = rfft[1::2].mean(axis=0)
        acc += np.vstack([rfft_tx0, rfft_tx1])
    return (acc / n_frames).astype(np.complex64)


def _sar_range_profile_from_blob(
    blob: bytes, bg_reference: np.ndarray | None = None
) -> np.ndarray:
    prof = _virtual_range_profile_raw(blob)
    if bg_reference is not None:
        prof = prof - bg_reference
    return prof


def _ra_power_from_profile(profile: np.ndarray) -> np.ndarray:
    pos_m = sf.virtual_array_positions_m()
    steering = np.exp(
        -1j * 2 * np.pi / sf.WAVELENGTH
        * np.sin(sf.ANGLES_RAD)[:, np.newaxis]
        * pos_m[np.newaxis, :]
    )
    return np.abs(steering @ profile) ** 2


# ── UDP / frame assembly ─────────────────────────────────────

def _udp_receiver(raw_queue: queue.Queue[bytes], stop: threading.Event) -> None:
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 2**26) 
    sock.bind((PC_IP, DATA_PORT))
    sock.settimeout(0.5)
    last_seq = -1
    print(f"UDP listening on {PC_IP}:{DATA_PORT} ...")
    while not stop.is_set():
        try:
            packet, _ = sock.recvfrom(4096)
            seq = struct.unpack("<I", packet[0:4])[0]
            data = packet[10:]
            if last_seq >= 0 and seq != last_seq + 1:
                gap = seq - last_seq - 1
                if gap > 0:
                    print(f"[warn] UDP gap: dropped ~{gap} packet(s)")
            last_seq = seq
            try:
                raw_queue.put_nowait(data)
            except queue.Full:
                pass
        except socket.timeout:
            continue
        except OSError:
            break
    sock.close()


def _frame_assembler(
    
    raw_queue: queue.Queue[bytes],
    frame_ring: deque,          
    ring_lock: threading.Lock,
    stop: threading.Event,
    bytes_per_frame: int,
) -> None:
   
    buffer = bytearray()
    frame_count = 0
    print("Frame assembler running.")
    while not stop.is_set():
        try:
            data = raw_queue.get(timeout=0.25)
        except queue.Empty:
            continue
        buffer.extend(data)

        # Drain queue to prevent overflow during slow backprojection
        while not raw_queue.empty():
            try:
                buffer.extend(raw_queue.get_nowait())
            except queue.Empty:
                break
        # In _frame_assembler, track the gap between frames:
        import time
        last_frame_time = time.monotonic()

        # After frame_count increment:
        now = time.monotonic()
        gap = now - last_frame_time
        if gap > 0.5:  
            print(f"[assembler] FRAME GAP: {gap:.2f}s after frame {frame_count}")
        last_frame_time = now

        while len(buffer) >= bytes_per_frame:
            raw = bytes(buffer[:bytes_per_frame])
            del buffer[:bytes_per_frame]
            with ring_lock:
                frame_ring.append((frame_count, raw)) 
            frame_count += 1
            if frame_count % 20 == 0:
                print(f"[assembler] {frame_count} frames  |  "
                      f"ring: {len(frame_ring)}  |  "
                      f"queue: {raw_queue.qsize()}")


def _wait_snapshot_blob_after(
    frame_ring: deque,
    ring_lock: threading.Lock,
    min_frame_id: int,
    deadline_s: float = 60.0,
) -> tuple[bytes, int] | None:
    deadline = time.monotonic() + deadline_s
    last_log = time.monotonic()
    while time.monotonic() < deadline:
        with ring_lock:
            all_ids = [fid for fid, _ in frame_ring]
            available = [(fid, raw) for fid, raw in frame_ring
                         if fid >= min_frame_id]
        
        now = time.monotonic()
        if now - last_log >= 2.0:
            min_id_in_ring = min(all_ids) if all_ids else -1
            max_id_in_ring = max(all_ids) if all_ids else -1
            print(f"[wait] need {FRAMES_PER_SNAPSHOT} frames with id>={min_frame_id}  |  "
                  f"ring has {len(all_ids)} frames  |  "
                  f"ids {min_id_in_ring}–{max_id_in_ring}  |  "
                  f"available: {len(available)}")
            last_log = now

        if len(available) >= FRAMES_PER_SNAPSHOT:
            selected = available[-FRAMES_PER_SNAPSHOT:]
            last_id = selected[-1][0]
            blob = b"".join(raw for _, raw in selected)
            return blob, last_id
        time.sleep(0.05)
    
    with ring_lock:
        all_ids = [fid for fid, _ in frame_ring]
    print(f"[wait] TIMEOUT — min_frame_id={min_frame_id}  |  "
          f"ring ids: {min(all_ids) if all_ids else 'empty'}–"
          f"{max(all_ids) if all_ids else 'empty'}")
    return None
# ── SAR image formation ──────────────────────────────────────

def finalize_sar_image(
    sar_captures: list[np.ndarray],
    sar_offsets_m: list[float],
) -> None:
    k = len(sar_captures)
    if k < 2:
        print(f"[SAR] Need >= 2 snapshots (have {k}).")
        return

    offsets = np.asarray(sar_offsets_m, dtype=np.float64)
    pos_local = sf.virtual_array_positions_m()
    n_virt = len(pos_local)
    positions_all = np.concatenate([offsets[i] + pos_local for i in range(k)])
    aperture_span = float(positions_all.max() - positions_all.min())
    theta_res_deg = (float(np.rad2deg(sf.WAVELENGTH / aperture_span))
                     if aperture_span > 0 else 0.0)

    print(
        f"[SAR] mode={SAR_MODE.upper()}  |  {k} positions  |  "
        f"{k * n_virt} virtual elements  |  "
        f"aperture ≈ {aperture_span * 1e3:.1f} mm  |  "
        f"θ_res ≈ {theta_res_deg:.1f}°"
    )

    # Phase autofocus with unwrapped position refinement
    sar_captures, sar_offsets_m = autofocus_with_position_refinement(
        sar_captures, list(sar_offsets_m)
    )

    # Angle grid — standard -90 to +90
    sar_angles_deg = np.linspace(-90.0, 90.0, N_SAR_ANGLES)
    sar_angles_rad = np.deg2rad(sar_angles_deg)
    range_axis = np.arange(sf.NUM_RANGE_BINS) * sf.RANGE_RES

    # ── Image formation ──────────────────────────────────────
    if SAR_MODE == "bp":
        t0 = time.monotonic()
        ra_power = backproject_sar(sar_captures, sar_offsets_m, sar_angles_rad)
        print(f"[BP] Complete in {time.monotonic() - t0:.1f} s")
        mode_label = "Backprojection"
    else:
        signal_ext = np.vstack(sar_captures)
        sort_idx = np.argsort(positions_all)
        positions_sorted = -positions_all[sort_idx]  
        signal_sorted = signal_ext[sort_idx, :]
        win = np.hamming(len(positions_sorted))
        signal_sorted = signal_sorted * win[:, np.newaxis]
        steering_ext = np.exp(
            -1j * 2 * np.pi / sf.WAVELENGTH
            * np.sin(sar_angles_rad)[:, np.newaxis]
            * positions_sorted[np.newaxis, :]
        )
        ra_power = np.abs(steering_ext @ signal_sorted) ** 2
        mode_label = "FFT beamforming"

    ra_db = 20 * np.log10(ra_power + 1e-12)
    ra_db_peak = float(np.nanmax(ra_db))
    ra_db_norm = ra_db - ra_db_peak

    # ── Cartesian reproject ──────────────────────────────────
    aa, rr = np.meshgrid(sar_angles_rad, range_axis, indexing="ij")
    x_flat = (rr * np.sin(aa)).ravel()
    y_flat = (rr * np.cos(aa)).ravel()
    z_flat = ra_db_norm.ravel()
    x_grid_1d = np.linspace(float(x_flat.min()), float(x_flat.max()), 600)
    y_grid_1d = np.linspace(max(float(y_flat.min()), 0.01),
                             float(y_flat.max()), 400)
    x_grid, y_grid = np.meshgrid(x_grid_1d, y_grid_1d)
    cart_image = griddata(
        (x_flat, y_flat), z_flat, (x_grid, y_grid), method="linear"
    )

    # ── Plot ─────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
    im1 = ax1.imshow(
        ra_db_norm.T, aspect="auto", origin="lower", cmap="inferno",
        extent=[sar_angles_deg[0], sar_angles_deg[-1], 0, sf.MAX_RANGE],
        interpolation="bilinear", vmin=-60, vmax=0,
    )
    ax1.set_xlabel("Angle (degrees)")
    ax1.set_ylabel("Range (m)")
    ax1.set_ylim(0.5, min(5.0, sf.MAX_RANGE))
    ax1.set_title(
        f"SAR range-angle  ({k} positions, {k * n_virt} virtual elements)\n"
        f"Mode: {mode_label}  |  Autofocus: ON  |  θ_res ≈ {theta_res_deg:.1f}°"
    )
    plt.colorbar(im1, ax=ax1, fraction=0.03, pad=0.04, label="dB re peak")

    im2 = ax2.imshow(
        cart_image, aspect="auto", origin="lower", cmap="inferno",
        extent=[x_grid_1d[0], x_grid_1d[-1], y_grid_1d[0], y_grid_1d[-1]],
        interpolation="bilinear", vmin=-60, vmax=0,
    )
    ax2.set_xlabel("Cross-range x (m)")
    ax2.set_ylabel("Down-range y (m)")
    ax2.set_ylim(0.5, min(5.0, sf.MAX_RANGE))
    ax2.set_title("SAR Cartesian image")
    plt.colorbar(im2, ax=ax2, fraction=0.03, pad=0.04, label="dB re peak")

    ref_str = (f"{AUTOFOCUS_REF_RANGE_M:.2f} m"
               if AUTOFOCUS_REF_RANGE_M is not None else "auto")
    fig.suptitle(
        f"Coherent SAR  |  {mode_label}  |  {k} positions  |  "
        f"λ = {sf.WAVELENGTH * 1e3:.2f} mm  |  "
        f"aperture ≈ {aperture_span * 1e3:.1f} mm  |  "
        f"θ_res ≈ {theta_res_deg:.1f}°  |  autofocus ref: {ref_str}",
        fontsize=12,
    )
    plt.tight_layout()

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base = SAR_OUTPUT_DIR if SAR_OUTPUT_DIR is not None else Path(".")
    base.mkdir(parents=True, exist_ok=True)
    png_path = base / f"sar_{SAR_MODE}_{ts}.png"
    npz_path = base / f"sar_{SAR_MODE}_{ts}.npz"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    print(f"[SAR] Figure saved → {png_path.resolve()}")
    np.savez_compressed(
        npz_path,
        captures=np.stack(sar_captures, axis=0),
        offsets_m=offsets,
        range_axis=range_axis,
        sar_angles_deg=sar_angles_deg,
        ra_db=ra_db_norm,
        ra_db_peak_offset_db=ra_db_peak,
        cart_image=cart_image,
        x_grid_1d=x_grid_1d,
        y_grid_1d=y_grid_1d,
        mode=SAR_MODE,
    )
    print(f"[SAR] Arrays saved → {npz_path.resolve()}")
    plt.show(block=False)


# ── Live plot ────────────────────────────────────────────────

def start_live_plot(
    frame_ring: deque,
    ring_lock: threading.Lock,
    stop: threading.Event,
    cmd_queue: queue.Queue[str | None],
    st: dict,
    st_lock: threading.Lock,
) -> None:
    fig, ax = plt.subplots(figsize=(11, 7))
    ny, nx = sf.NUM_RANGE_BINS, sf.NUM_ANGLE_BINS
    ra_img = ax.imshow(
        np.zeros((ny, nx)), aspect="auto", origin="lower", cmap="inferno",
        extent=[float(sf.ANGLES_DEG[0]), float(sf.ANGLES_DEG[-1]),
                sf.DISPLAY_RANGE_MIN_M, sf.MAX_RANGE],
        interpolation="bilinear",
    )
    scatter = ax.scatter(
        [], [], c="cyan", s=120, marker="+", linewidths=2.0, zorder=5
    )
    ax.set_xlabel("Angle (degrees)", fontsize=12)
    ax.set_ylabel("Range (m)", fontsize=12)
    ax.set_title(
        f"updatedsar — IWR + DCA1000  |  SAR mode: {SAR_MODE.upper()}",
        fontsize=13
    )
    ax.set_ylim(sf.DISPLAY_RANGE_MIN_M, sf.DISPLAY_RANGE_MAX_M)
    ax.grid(True, alpha=0.15, color="white")
    plt.colorbar(
        ra_img, ax=ax, fraction=0.03, pad=0.04,
        label="Rel. power (dB re peak)"
    )

    info = ax.text(
        0.01, 0.99, "", transform=ax.transAxes, color="white",
        fontsize=9, va="top",
        bbox=dict(facecolor="black", alpha=0.5, pad=3),
    )
    _ = ax.text(
        0.5, 0.01,
        "'b'=bg  |  ←/→=move λ/2+capture  |  'c'=capture  |  "
        "'r'=reset pos to 0  |  'u'=undo  |  'f'=finalize  |  'q'=quit",
        transform=ax.transAxes, color="yellow", fontsize=8,
        ha="center", va="bottom",
        bbox=dict(facecolor="black", alpha=0.6, pad=4),
    )

    display_queue: queue.Queue[dict] = queue.Queue()

    def on_key(event):
        if event.key is None:
            return
        mapping = {
            "q": None, "b": "bg", "right": "right", "left": "left",
            "c": "here", "u": "undo", "f": "finalize", "d": "diag", "r": "reset",
        }
        if event.key in mapping:
            cmd_queue.put(mapping[event.key])
            if event.key == "q":
                plt.close("all")

    fig.canvas.mpl_connect("key_press_event", on_key)

    def on_close(_event):
        if not stop.is_set():
            cmd_queue.put(None)
            stop.set()

    fig.canvas.mpl_connect("close_event", on_close)

    last_display: dict = {"ra_db": None, "dets": None, "lam": 0}

    def update(_):
        while True:
            try:
                msg = display_queue.get_nowait()
            except queue.Empty:
                break
            if msg["type"] == "snapshot":
                last_display.update(
                    ra_db=msg["ra_db"], dets=msg["dets"], lam=msg["lam"]
                )

        rdb = last_display["ra_db"]
        if rdb is not None:
            disp = (rdb - np.max(rdb)).T
            ra_img.set_data(disp)
            finite = disp[np.isfinite(disp)]
            if finite.size:
                lo = np.percentile(finite, 5)
                hi = np.percentile(finite, 99.5)
                ra_img.set_clim(lo, hi)

        dets = last_display["dets"]
        if dets:
            scatter.set_offsets(
                np.c_[[d[1] for d in dets], [d[0] for d in dets]]
            )
        else:
            scatter.set_offsets(np.empty((0, 2)))

        with st_lock:
            lam = st["along_track_lam"]
            ncap = len(st["sar_captures"])
            bg_ok = st["bg_reference"] is not None
            status = st.get("status", "Ready")

        info.set_text(
            f"SAR positions: {ncap}  |  x = {lam:+d} \u00d7 \u03bb/2  |  "
            f"BG ref: {'YES' if bg_ok else 'NO'}  |  {status}"
        )
        return ra_img, scatter, info

    def worker() -> None:
        while True:
            try:
                kind = cmd_queue.get(timeout=0.2)
            except queue.Empty:
                if stop.is_set():
                    break
                continue
            if kind is None:
                break
            try:
                with st_lock:
                    st["status"] = "Working..."

                # finalize
                if kind == "finalize":
                    with st_lock:
                        caps = list(st["sar_captures"])
                        offs = list(st["sar_offsets_m"])
                    if len(caps) < 2:
                        print(f"[SAR] Need >= 2 positions (have {len(caps)}).")
                    else:
                        finalize_sar_image(caps, offs)
                    with st_lock:
                        st["status"] = "Ready"
                    continue

                # undo
                if kind == "undo":
                    with st_lock:
                        if not st["sar_captures"]:
                            print("[SAR] Nothing to undo.")
                            st["status"] = "Ready"
                            continue
                        st["sar_captures"].pop()
                        st["sar_offsets_m"].pop()
                        st["snapshot_id"] = max(0, st["snapshot_id"] - 1)
                        st["along_track_lam"] = (
                            int(round(st["sar_offsets_m"][-1]
                                      / (sf.WAVELENGTH / 2)))
                            if st["sar_offsets_m"] else 0
                        )
                        lam = st["along_track_lam"]
                        nrem = len(st["sar_captures"])
                        st["status"] = "Ready"
                    print(f"[SAR] Undid last. {nrem} remain. "
                          f"x={lam:+d} × λ/2")
                    continue

                # Add after the undo handler:
                if kind == "reset":
                    with st_lock:
                        st["along_track_lam"] = 0
                        st["status"] = "Ready"
                    print("[SAR] Position counter reset to 0 — ready for right sweep")
                    continue

                # diagnostics
                if kind == "diag":
                    with st_lock:
                        last_id = st["last_frame_id"]
                    result = _wait_snapshot_blob_after(
                        frame_ring, ring_lock, last_id
                    )
                    if result is None:
                        print("[diag] Timeout waiting for fresh frames.")
                    else:
                        blob, _ = result
                        sf.print_diagnostics(blob)
                    with st_lock:
                        st["status"] = "Ready"
                    continue

                # move / capture / background
                with st_lock:
                    if kind == "right":
                        st["along_track_lam"] += 1
                    elif kind == "left":
                        st["along_track_lam"] -= 1
                    elif kind not in ("here", "bg"):
                        st["status"] = "Ready"
                        continue
                    lam = st["along_track_lam"]
                    pos_idx = len(st["sar_captures"]) + 1
                    last_id = st["last_frame_id"]

                label = {
                    "right": "RIGHT", "left": "LEFT",
                    "here": "HERE", "bg": "BACKGROUND"
                }.get(kind, kind.upper())
                print(
                    f"[SAR] [{label}] Waiting for {FRAMES_PER_SNAPSHOT} "
                    f"fresh frames (min_id={last_id})"
                    + (f" — position {pos_idx} (x={lam:+d} × λ/2)"
                       if kind != "bg" else "")
                    + "..."
                )

                result = _wait_snapshot_blob_after(
                    frame_ring, ring_lock, last_id
                )
                if result is None:
                    print("[SAR] Timeout waiting for fresh frames. "
                          "Check radar is still transmitting.")
                    with st_lock:
                        st["status"] = "Ready"
                    continue

                blob, new_last_id = result

                # Advance last_frame_id so next capture can't reuse these frames
                with st_lock:
                    st["last_frame_id"] = new_last_id + 1

                if kind == "bg":
                    raw = _virtual_range_profile_raw(blob)
                    with st_lock:
                        st["bg_reference"] = raw.copy()
                        st["status"] = "Ready"
                    print(f"[SAR] Background captured "
                          f"(frames up to id={new_last_id}).")
                    continue

                if SNAPSHOT_SAVE_DIR is not None:
                    SNAPSHOT_SAVE_DIR.mkdir(parents=True, exist_ok=True)
                    with st_lock:
                        sid = st["snapshot_id"]
                    out = (SNAPSHOT_SAVE_DIR
                           / f"aperture_{sid:04d}_lam{lam:+d}.bin")
                    out.write_bytes(blob)
                    print(f"Saved {out}")

                with st_lock:
                    bg_ref = st["bg_reference"]

                profile = _sar_range_profile_from_blob(
                    blob, bg_reference=bg_ref
                )
                ra_power = None
                ra_db = None
                dets = []

                offset_m = float(lam) * float(sf.WAVELENGTH) / 2.0
                with st_lock:
                    st["sar_captures"].append(profile)
                    st["sar_offsets_m"].append(offset_m)
                    st["snapshot_id"] += 1
                    ncap = len(st["sar_captures"])

                print(f"[SAR] Stored snapshot {ncap} @ {lam:+d} × λ/2 "
                      f"= {offset_m*1000:.2f} mm  "
                      f"(frame ids {last_id}–{new_last_id}).")
                for i, (rng_m, ang_deg, p) in enumerate(dets, start=1):
                    print(f"  det {i}: R={rng_m:.3f} m  "
                          f"\u2220={ang_deg:+.1f}\u00b0  P={p:.3e}")
                if not dets:
                    print("  (no detections above threshold)")

                if ra_db is not None:
                    display_queue.put({
                        "type": "snapshot", "ra_db": ra_db,
                        "dets": dets, "lam": lam
                    })
                else:
                    display_queue.put({
                        "type": "snapshot", "ra_db": None,
                        "dets": [], "lam": lam
                    })
                with st_lock:
                    st["status"] = "Ready"

            except Exception as e:
                print(f"[SAR] Error: {e}")
                import traceback
                traceback.print_exc()
                with st_lock:
                    st["status"] = "Ready"

    worker_thread = threading.Thread(target=worker, daemon=True)
    worker_thread.start()

    _ani = animation.FuncAnimation(
        fig, update, interval=100, blit=True, cache_frame_data=False,
    )
    plt.tight_layout()
    plt.show()

    stop.set()
    cmd_queue.put(None)
    worker_thread.join(timeout=5.0)


def main() -> None:
    print("IWR + DCA1000  updatedsar (virtual ULA + SAR buffer)")
    print(f"Config: {sf.CONFIG_XML.name}  |  \u03bb = {sf.WAVELENGTH * 1e3:.3f} mm")
    print(f"Frame size: {sf.BYTES_PER_FRAME / 1024:.1f} KB  |  "
          f"frames/capture: {FRAMES_PER_SNAPSHOT}")
    print(f"SAR mode: {SAR_MODE.upper()}"
          + (" [numba JIT]" if _NUMBA else
             " [numpy fallback — install numba for speed]"))
    print(f"Step size: λ/2 = {sf.WAVELENGTH / 2 * 1e3:.2f} mm per arrow key")
    ref_str = (f"{AUTOFOCUS_REF_RANGE_M:.2f} m"
               if AUTOFOCUS_REF_RANGE_M is not None else "auto (strongest bin)")
    print(f"Autofocus reference: {ref_str}")
    if AUTOFOCUS_REF_RANGE_M is not None:
        print(f"  *** Place a static metal reflector at "
              f"{AUTOFOCUS_REF_RANGE_M:.2f} m ***")
        print(f"  *** at a DIFFERENT range than your target ***")
    print()
    print("Keys (focus the plot window):")
    print("  'b'  = capture background reference")
    print("  ←/→  = move one λ/2 + capture")
    print("  'c'  = capture at current position")
    print("  'd'  = diagnostics")
    print("  'u'  = undo last capture")
    print("  'f'  = finalize SAR image")
    print("  'q'  = quit\n")

    raw_queue: queue.Queue[bytes] = queue.Queue(maxsize=4096)
    frame_ring: deque = deque(maxlen=RING_MAX_FRAMES)
    ring_lock = threading.Lock()
    stop = threading.Event()

    recv_thread = threading.Thread(
        target=_udp_receiver, args=(raw_queue, stop), daemon=True
    )
    asm_thread = threading.Thread(
        target=_frame_assembler,
        args=(raw_queue, frame_ring, ring_lock, stop, sf.BYTES_PER_FRAME),
        daemon=True,
    )
    recv_thread.start()
    asm_thread.start()

    print("Testing UDP reception for 5 seconds...")
    time.sleep(5)
    with ring_lock:
        n = len(frame_ring)
    print(f"Frames in ring after 5s: {n}")
    if n == 0:
        print("ERROR: No frames received. Check:")
        print("  1. Lua script running and showing 'Streaming...'")
        print("  2. PC firewall not blocking UDP port 4098")
        print("  3. PC IP is 192.168.33.30")
        print("  4. DCA1000 IP is 192.168.33.180")
    else:
        print(f"OK: receiving data ({n} frames buffered)")

    time.sleep(1)

    st: dict = {
        "along_track_lam": 0,
        "snapshot_id": 0,
        "sar_captures": [],
        "sar_offsets_m": [],
        "bg_reference": None,
        "status": "Ready",
        "last_frame_id": 0,     
    }
    st_lock = threading.Lock()
    cmd_queue: queue.Queue[str | None] = queue.Queue()

    try:
        start_live_plot(frame_ring, ring_lock, stop, cmd_queue, st, st_lock)
    finally:
        stop.set()
        recv_thread.join(timeout=2.0)
        asm_thread.join(timeout=2.0)
        print("Stopped.")


if __name__ == "__main__":
    main()