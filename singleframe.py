from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_SCRIPT_DIR = Path(__file__).resolve().parent
CONFIG_XML = _SCRIPT_DIR / "2txconfig.xml"

# Physical constants
SPEED_OF_LIGHT = 3e8

# Antenna geometry 
NUM_RX = 4
NUM_TX = 2
IQ = 2
BYTES_PER_SAMPLE = 2


def _params_under(root: ET.Element, section_tag: str) -> dict[str, str]:
    sec = root.find(section_tag)
    if sec is None:
        return {}
    return {p.get("name", ""): p.get("value", "") for p in sec.findall("param")}


def apply_config_from_xml(path: Path | None = None) -> None:
    global NUM_ADC, NUM_LOOPS, NUM_CHIRPS, BYTES_PER_FRAME
    global FREQ_SLOPE_MHZ_US, SAMPLE_RATE, START_FREQ, WAVELENGTH, ANT_SPACING, TX_SEP_M
    global CHIRP_SLOPE_HZ_S, RANGE_RES, NUM_RANGE_BINS, MAX_RANGE
    global MIN_RANGE_BIN, MAX_RANGE_BIN, NUM_ANGLE_BINS, ANGLES_DEG, ANGLES_RAD

    path = path or CONFIG_XML
    if not path.is_file():
        raise FileNotFoundError(f"Radar config not found: {path}")

    tree = ET.parse(path)
    root = tree.getroot()
    prof = _params_under(root, "apiname_profile_cfg")
    frame = _params_under(root, "apiname_frame_cfg")

    FREQ_SLOPE_MHZ_US = float(prof["freqSlopeConst"])
    NUM_ADC = int(prof["numAdcSamples"])
    SAMPLE_RATE = float(prof["digOutSampleRate"]) * 1e3
    START_FREQ = float(prof["startFreqConst"]) * 1e9

    NUM_LOOPS = int(frame["loopCount"])

    WAVELENGTH = SPEED_OF_LIGHT / START_FREQ
    ANT_SPACING = WAVELENGTH / 2.0
    TX_SEP_M = 4.0 * ANT_SPACING

    NUM_CHIRPS = NUM_TX * NUM_LOOPS
    BYTES_PER_FRAME = NUM_RX * NUM_CHIRPS * NUM_ADC * IQ * BYTES_PER_SAMPLE

    CHIRP_SLOPE_HZ_S = FREQ_SLOPE_MHZ_US * 1e12
    RANGE_RES = (SPEED_OF_LIGHT * SAMPLE_RATE) / (2.0 * CHIRP_SLOPE_HZ_S * float(NUM_ADC))

    NUM_RANGE_BINS = NUM_ADC // 2
    MAX_RANGE = RANGE_RES * NUM_RANGE_BINS

    MIN_RANGE_BIN = max(1, int(0.5 / RANGE_RES))
    MAX_RANGE_BIN = min(NUM_RANGE_BINS - 1, int(5.0 / RANGE_RES))

    NUM_ANGLE_BINS = 180
    ANGLES_DEG = np.linspace(-90, 90, NUM_ANGLE_BINS)
    ANGLES_RAD = np.deg2rad(ANGLES_DEG)


apply_config_from_xml()

DISPLAY_RANGE_MIN_M = 0.0
DISPLAY_RANGE_MAX_M = 6.0

ADC_BIN = Path(r"C:\ti\radar_data\adc_data_Raw_0.bin")

def raw_bytes_to_iq(raw: bytes) -> np.ndarray:
    samples = np.frombuffer(raw, dtype=np.int16).copy()
    raw4 = samples.reshape(-1, 4)
    lvds = np.empty(len(raw4) * 2, dtype=np.complex64)
    lvds[0::2] = raw4[:, 0].astype(np.float32) + 1j * raw4[:, 2].astype(np.float32)
    lvds[1::2] = raw4[:, 1].astype(np.float32) + 1j * raw4[:, 3].astype(np.float32)
    return lvds.reshape(NUM_CHIRPS, NUM_RX, NUM_ADC)


def virtual_array_positions_m() -> np.ndarray:
    rx_offsets = np.arange(NUM_RX, dtype=np.float64) * ANT_SPACING
    return np.concatenate([rx_offsets, TX_SEP_M + rx_offsets])


def range_angle_map_multiframe(blob: bytes) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    n_frames = len(blob) // BYTES_PER_FRAME
    if n_frames == 0:
        raise ValueError("No full frames found in file.")

    win_range = np.blackman(NUM_ADC).astype(np.float32)
    pos_m = virtual_array_positions_m()

    steering = np.exp(
        -1j * 2 * np.pi / WAVELENGTH
        * np.sin(ANGLES_RAD)[:, np.newaxis]
        * pos_m[np.newaxis, :]
    )  

    ra_power_sum = None

    for k in range(n_frames):
        iq = raw_bytes_to_iq(blob[k * BYTES_PER_FRAME:(k + 1) * BYTES_PER_FRAME])

        # Range FFT with Blackman window for low sidelobes
        rfft = np.fft.fft(
            iq * win_range[np.newaxis, np.newaxis, :], axis=2
        )[:, :, :NUM_RANGE_BINS] 

        rfft_tx0 = rfft[0::2].mean(axis=0)  
        rfft_tx1 = rfft[1::2].mean(axis=0) 

        range_virt = np.vstack([rfft_tx0, rfft_tx1])  

        ra_power = np.abs(steering @ range_virt) ** 2  

        ra_power_sum = ra_power if ra_power_sum is None else ra_power_sum + ra_power

    ra_power_avg = ra_power_sum / n_frames
    ra_db = 20 * np.log10(ra_power_avg + 1e-12)
    range_axis = np.arange(NUM_RANGE_BINS) * RANGE_RES
    return ra_db, range_axis, ra_power_avg, n_frames


def detect_strongest_target(
    ra_power: np.ndarray,
    min_r: int = MIN_RANGE_BIN,
    max_r: int = MAX_RANGE_BIN,
    max_angle_deg: float = 60.0,
) -> tuple[float, float, int, int]:
    gated = ra_power.copy()

    gated[:, :min_r] = -np.inf
    gated[:, max_r + 1:] = -np.inf

    # Angle gate: ignore beyond +/-60 deg to suppress grating lobes
    angle_mask = np.abs(ANGLES_DEG) > max_angle_deg
    gated[angle_mask, :] = -np.inf

    flat = np.argmax(gated)
    a_idx, r_idx = np.unravel_index(flat, gated.shape)
    range_m = float(r_idx * RANGE_RES)
    angle_deg = float(-ANGLES_DEG[a_idx])
    return range_m, angle_deg, int(a_idx), int(r_idx)


def print_diagnostics(blob: bytes) -> None:
    n_frames = len(blob) // BYTES_PER_FRAME
    iq = np.stack(
        [raw_bytes_to_iq(blob[k * BYTES_PER_FRAME:(k + 1) * BYTES_PER_FRAME])
         for k in range(n_frames)], axis=0
    ).mean(axis=0)

    chirp0 = iq[0, 0, :]
    chirp1 = iq[1, 0, :]
    chirp2 = iq[2, 0, :]
    print("Chirp 0 mean power:", np.abs(chirp0).mean())
    print("Chirp 1 mean power:", np.abs(chirp1).mean())
    print("Chirp 0 vs 1 correlation:", np.abs(np.dot(chirp0, np.conj(chirp1))) / (np.abs(chirp0).mean() * NUM_ADC))
    print("Chirp 0 vs 2 correlation:", np.abs(np.dot(chirp0, np.conj(chirp2))) / (np.abs(chirp0).mean() * NUM_ADC))

    win_range = np.blackman(NUM_ADC).astype(np.float32)
    rfft = np.fft.fft(iq * win_range[np.newaxis, np.newaxis, :], axis=2)[:, :, :NUM_RANGE_BINS]
    rfft_tx0 = rfft[0::2]
    rfft_tx1 = rfft[1::2]
    p0 = np.abs(rfft_tx0).mean()
    p1 = np.abs(rfft_tx1).mean()
    corr = np.abs(np.mean(rfft_tx0 * np.conj(rfft_tx1)))
    print(f"TX0 power: {p0:.2f}")
    print(f"TX1 power: {p1:.2f}")
    print(f"Normalized TX0/TX1 correlation: {corr / (p0 * p1):.4f}")


def main() -> None:
    print(
        f"Config: {CONFIG_XML.name}  |  slope={FREQ_SLOPE_MHZ_US} MHz/us  |  "
        f"f_s={SAMPLE_RATE/1e6:.3f} MS/s  |  N_ADC={NUM_ADC}  |  loops={NUM_LOOPS}  |  "
        f"R_res={RANGE_RES*100:.2f} cm/bin  |  R_max={MAX_RANGE:.2f} m"
    )

    blob = ADC_BIN.read_bytes()
    if len(blob) < BYTES_PER_FRAME:
        raise SystemExit(
            f"{ADC_BIN} is too small for one frame "
            f"(need {BYTES_PER_FRAME} bytes, got {len(blob)})"
        )
    if len(blob) % BYTES_PER_FRAME != 0:
        print(f"Warning: trailing {len(blob) % BYTES_PER_FRAME} bytes ignored (not a full frame)")

    print_diagnostics(blob)

    ra_db, _range_axis, ra_power, n_frames = range_angle_map_multiframe(blob)

    tgt_range_m, tgt_angle_deg, _a_idx, _r_idx = detect_strongest_target(ra_power)
    print(
        f"Detected target (strongest peak in range bins {MIN_RANGE_BIN} to {MAX_RANGE_BIN}, "
        f"angles +/-60 deg): range = {tgt_range_m:.3f} m,  angle = {tgt_angle_deg:+.2f} deg"
    )

    ra_db_norm = ra_db - np.max(ra_db)

    fig, ax = plt.subplots(figsize=(11, 7))
    
    display = ra_db_norm.T[:, ::-1]
    im = ax.imshow(
        display,
        aspect="auto",
        origin="lower",
        cmap="inferno",
        extent=[
            float(ANGLES_DEG[0]),
            float(ANGLES_DEG[-1]),
            DISPLAY_RANGE_MIN_M,
            MAX_RANGE,
        ],
        interpolation="bilinear",
        vmin=-60,
        vmax=0.0,
    )
    ax.set_xlabel("Angle (deg)")
    ax.set_ylabel("Range (m)")
    ax.plot(
        tgt_angle_deg,
        tgt_range_m,
        "+",
        color="cyan",
        markersize=14,
        markeredgewidth=2,
        label="detected peak",
    )
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title(
        f"Range-angle heatmap  |  {NUM_TX}x{NUM_RX} MIMO  |  {n_frames} frames (power avg)  |  "
        f"{ADC_BIN.name}\n"
        f"Detection: {tgt_range_m:.3f} m @ {tgt_angle_deg:+.1f} deg"
    )
    ax.set_ylim(DISPLAY_RANGE_MIN_M, DISPLAY_RANGE_MAX_M)
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