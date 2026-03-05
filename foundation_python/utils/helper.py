import numpy as np
import pandas as pd

import os
import matplotlib
import matplotlib.pyplot as plt
from numpy.lib.stride_tricks import sliding_window_view


def _ma7_2d(arr2d: np.ndarray) -> np.ndarray:
    df = pd.DataFrame(arr2d)
    out = df.rolling(window=7, min_periods=7).mean().to_numpy(dtype=np.float32)
    return out


def compute_2d_inci_peak(arr: np.ndarray, P_peak: dict) -> np.ndarray:
    """
    arr: (n_rows, n_cols) float array; each column is a time series.
    P_peak: {"peak_window": int (odd), "rel_thr": float}

    Returns:
      peaks: (n_rows - (half+6) - half, n_cols) uint8 array.
             Row i of 'peaks' corresponds to original row i + (half+6).
    """
    peak_window = int(P_peak["peak_window"])
    rel_thr     = float(P_peak["rel_thr"])
    assert peak_window % 2 == 1, "peak_window must be odd."
    half = (peak_window - 1) // 2

    n_rows, n_cols = arr.shape
    arr = arr.astype(np.float32, copy=False)

    ma = _ma7_2d(arr)  # (n_rows, n_cols) float32

    win = sliding_window_view(ma, window_shape=(2*half + 1), axis=0)
    # win shape: (n_rows - 2*half, n_cols, 2*half+1)

    if win.shape[0] <= 6:
        return np.zeros((0, n_cols), dtype=np.uint8)  # not enough data
    win = win[6:, :, :]  # (n_rows - 2*half - 6, n_cols, 2*half+1)

    center = win[:, :, half]                    # (n_valid, n_cols)
    nb_max = np.nanmax(win, axis=2)             # (n_valid, n_cols)

    past_min   = np.nanmin(win[:, :, :half], axis=2)      if half > 0 else np.full_like(center, np.nan, dtype=np.float32)
    future_min = np.nanmin(win[:, :, half+1:], axis=2)    if half > 0 else np.full_like(center, np.nan, dtype=np.float32)
    valid_sides = np.isfinite(past_min) & np.isfinite(future_min)

    is_local_max = (center >= nb_max) & np.isfinite(center) & (center > 0)

    prom_past   = center - past_min
    prom_future = center - future_min
    rel_ok = (prom_past >= (rel_thr * past_min)) & (prom_future >= (rel_thr * future_min))
    no_zero = (past_min > 0) & (future_min > 0)

    flags = (is_local_max & valid_sides & rel_ok & no_zero).astype(np.uint8)  # (n_valid, n_cols)

    return flags


def compute_windowed_inci_peak(
    inci_forecast: np.ndarray,            # (W, H=14, D) draws
    inci_observed_window: np.ndarray,     # (W, half+6, 1)
    P_peak: dict
) -> np.ndarray:
    """
    For each window w, build longer segments: [lead-in (half+6)] + [14 draws] + [next half draws]
    Then run peak detection and average 0/1 across draws to get probabilities.
    Returns: (W, 14, 1) float32 in [0,1].
    """
    peak_window = int(P_peak["peak_window"])
    half = (peak_window - 1) // 2

    W, H, D = inci_forecast.shape
    assert H == 14, "This function assumes horizon=14."
    assert inci_observed_window.shape == (W, half + 6, 1), "observed-window shape must be (W, half+6, 1)."

    out = np.zeros((W, H, 1), dtype=np.float32)

    for w in range(W):
        lead = inci_observed_window[w, :, 0]        # (half+6,)
        cur  = inci_forecast[w]                     # (14, D)
        nxt  = inci_forecast[w + 1, :half, :] if (w + 1 < W) else None  # (half, D) or None

        if nxt is not None:
            seg = np.concatenate([
                np.tile(lead[:, None], (1, D)),     # (half+6, D)
                cur,                                # (14, D)
                nxt                                 # (half, D)
            ], axis=0)                              # -> (L=2*half+20, D)
        else:
            right = np.zeros((half, D), dtype=np.float32)
            seg = np.concatenate([
                np.tile(lead[:, None], (1, D)),
                cur,
                right
            ], axis=0)

        flags = compute_2d_inci_peak(seg, P_peak).astype(np.float32)  # (14, D)
        probs = flags.mean(axis=1, dtype=np.float32)                   # (14,)

        if nxt is None:
            probs[-half:] = 0.0

        out[w, :, 0] = probs

    return out


def get_window_start_indices_from_dates(full_index: pd.DatetimeIndex, start_dates: np.ndarray) -> np.ndarray:
    """
    Map an array of Python date/datetime (start_dates) to integer row indices in full_index.
    Raises if any date is missing.
    full_index: 420 rows from 2020-01-19 to 2021-03-13 
    """
    # coerce to pandas Timestamps for robust matching
    sd = pd.to_datetime(start_dates)
    mapper = pd.Series(np.arange(len(full_index), dtype=np.int32), index=full_index)
    idx = mapper.reindex(sd)
    if idx.isna().any():
        missing = sd[idx.isna()]
        raise ValueError(f"Some start dates not found in index: {missing.tolist()}")
    return idx.to_numpy(dtype=np.int32)


def split_to_supple(
    inci_observed_long: np.ndarray,   # (n_rows, ) or (n_rows, 1) for one state
    start_indices: np.ndarray,        # (W,) window start row indices (t_i)
    P_peak: dict
) -> np.ndarray:
    """
    Returns observed lead-in for each window:
      shape = (W, half+6, 1), slices [t_i-(half+6) .. t_i-1]
    """
    peak_window = int(P_peak["peak_window"])
    half = (peak_window - 1) // 2

    y = np.asarray(inci_observed_long).reshape(-1, 1).astype(np.float32, copy=False)
    n_rows = y.shape[0]
    need = half + 6

    leads = []
    for t in start_indices:
        if t < need:
            raise ValueError(f"Window start {t} has insufficient past data (need {need}, have {t}).")
        leads.append(y[t - need: t, :])  # (need, 1)
    return np.stack(leads, axis=0)       # (W, need, 1)


def split_week_forecast(peak_window_forecast: np.ndarray,
                        start_indices: np.ndarray,
                        full_index: pd.DatetimeIndex,
                        window_len: int = 14,
                        stride_size: int = 7) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    peak_window_forecast: (W, 14, 1) float probs in [0,1] for each window day.
    Returns two DataFrames with date index:
      - one_week: forecast made 1 week ahead (closer window)
      - two_week: forecast made 2 weeks ahead (previous window’s second week)
    Note: two_week starts 7 days later; the first 7 days have no 2-week forecast.
    """
    W = peak_window_forecast.shape[0]
    H = window_len
    S = stride_size

    f = np.asarray(peak_window_forecast, dtype=np.float32).reshape(W, H)

    t0 = int(start_indices[0])
    t_last = int(start_indices[-1])
    end_pos = t_last + H

    # ---- 1-week ahead ----
    one = np.full((end_pos - t0,), np.nan, dtype=np.float32)
    # use first week [0..S-1] from each window
    for w, t in enumerate(start_indices):
        t = int(t)
        one[t - t0 : t - t0 + S] = f[w, :S]
    # plus the tail (second week) of the *last* window
    one[t_last - t0 + S : t_last - t0 + H] = f[-1, S:H]
    one_week = pd.DataFrame({'peak_prob_1w': one}, index=full_index[t0:end_pos])

    # ---- 2-week ahead ----
    # day d gets 2w forecast from previous window’s horizons [S..H-1]
    # therefore the series starts at start_indices[1]
    t1 = int(start_indices[1])  # first date with 2-week forecast
    two_len = end_pos - t1
    two = np.full((two_len,), np.nan, dtype=np.float32)

    for w in range(1, W):
        t = int(start_indices[w])
        # fill dates [t .. t+S-1] from previous window’s second week
        two[t - t1 : t - t1 + S] = f[w - 1, S:H]

    # no tail for 2-week: the last block already ends at end_pos
    two_week = pd.DataFrame({'peak_prob_2w': two}, index=full_index[t1:end_pos])

    return one_week, two_week


def sample_draws_from_quantiles(
    arr: np.ndarray, qt_levels: np.ndarray, n_draws: int = 2000, seed: int | None = None
) -> np.ndarray:
    """
    arr: (W=51, H=14, Q=19) quantiles; strictly increasing quantile levels required.
    qt_levels: shape (Q,)
    Returns:
      draws: (W, H, D) float32 with iid inverse-CDF sampling using linear interpolation.
      Tails (u<min_q, u>max_q) are clamped to end values.
    """
    W, H, Q = arr.shape
    qt_levels = np.asarray(qt_levels, dtype=np.float64)
    assert np.all(np.diff(qt_levels) > 0), "qt_levels must be strictly increasing and match arr's last axis."
    assert Q == qt_levels.size, "arr.shape[-1] must equal len(qt_levels)"

    rng = np.random.default_rng(seed)
    U = rng.random((W, H, n_draws), dtype=np.float64) # U: (W, H, D) in (0,1)

    # Extend quantile grid to [0,1] by clamping endpoints
    qx = np.concatenate(([0.0], qt_levels, [1.0]))

    # Flatten over (W,H) so we interpolate 714 short curves — very fast
    arr_flat = arr.reshape(W * H, Q).astype(np.float64, copy=False)
    U_flat   = U.reshape(W * H, n_draws)

    draws_flat = np.empty_like(U_flat, dtype=np.float32)

    for i in range(arr_flat.shape[0]):  # 714 iterations
        qy = arr_flat[i]
        # enforce monotone (defensive against tiny crossing)
        qy = np.maximum.accumulate(qy)
        qy_ext = np.concatenate(([qy[0]], qy, [qy[-1]]))
        draws_flat[i] = np.interp(U_flat[i], qx, qy_ext).astype(np.float32, copy=False)

    draws = draws_flat.reshape(W, H, n_draws)
    # Clip tiny negatives that sometimes appear around zero due to interpolation/tails
    np.maximum(draws, 0.0, out=draws)
    return draws.astype(np.float32, copy=False)


