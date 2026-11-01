"""Place-field computation, lifted directly from notebooks/02a-spike-stimulus-analysis.py.

Only generalized to accept arbitrary spike times / position arrays instead of
hard-coded globals, and to restrict input to a caller-supplied trial mask
(e.g., "Pre-Cooling" trials only) before binning. The occupancy-normalization
math is unchanged from the course notebook.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.ndimage import gaussian_filter


@dataclass
class RateMap:
    firing_rate_map: NDArray[np.floating]  # (n_x_bins, n_y_bins), Hz
    occupancy_time: NDArray[np.floating]  # seconds
    x_bins: NDArray[np.floating]
    y_bins: NDArray[np.floating]
    n_spikes_used: int


MIN_OCCUPANCY_SECONDS = 0.1  # Muller & Kubie (1987) convention: bins visited too
# briefly give unreliable rate estimates -- e.g. at this dataset's ~120 Hz position
# sampling rate, a single-sample bin (0.0083s) with one incidental spike reads as a
# 120 Hz "peak", which is pure sampling noise, not tuning. Below-threshold bins are
# treated as unvisited (excluded from rate map, spatial info, and peak rate).


def compute_rate_map(
    spike_times: NDArray[np.floating],
    x_position: NDArray[np.floating],
    y_position: NDArray[np.floating],
    position_timestamps: NDArray[np.floating],
    pos_bin_size: float = 2.0,
    min_occupancy_seconds: float = MIN_OCCUPANCY_SECONDS,
) -> RateMap:
    """Occupancy-normalized firing rate map. Same method as notebook 02a,
    plus a minimum-occupancy mask (not present in the notebook) to exclude
    sparsely-sampled bins from biasing peak rate / spatial information."""
    epoch_spike_times = spike_times[
        (spike_times >= position_timestamps.min())
        & (spike_times <= position_timestamps.max())
    ]

    position_ind_at_spike = np.searchsorted(position_timestamps, epoch_spike_times)
    position_ind_at_spike = np.clip(
        position_ind_at_spike, 0, len(position_timestamps) - 1
    )
    x_at_spike = x_position[position_ind_at_spike]
    y_at_spike = y_position[position_ind_at_spike]

    min_x, max_x = np.nanmin(x_position), np.nanmax(x_position)
    min_y, max_y = np.nanmin(y_position), np.nanmax(y_position)
    x_bins = np.arange(min_x, max_x + pos_bin_size, pos_bin_size)
    y_bins = np.arange(min_y, max_y + pos_bin_size, pos_bin_size)

    binned_occupancy, _, _ = np.histogram2d(x_position, y_position, bins=[x_bins, y_bins])
    binned_spikes, _, _ = np.histogram2d(x_at_spike, y_at_spike, bins=[x_bins, y_bins])

    sampling_rate = 1 / np.median(np.diff(position_timestamps))
    occupancy_time = binned_occupancy / sampling_rate

    with np.errstate(divide="ignore", invalid="ignore"):
        firing_rate_map = binned_spikes / occupancy_time
        firing_rate_map[~np.isfinite(firing_rate_map)] = 0.0

    under_sampled = (occupancy_time > 0) & (occupancy_time < min_occupancy_seconds)
    firing_rate_map[under_sampled] = 0.0
    occupancy_time = occupancy_time.copy()
    occupancy_time[under_sampled] = 0.0  # excluded from spatial-info occupancy weighting too

    return RateMap(
        firing_rate_map=firing_rate_map,
        occupancy_time=occupancy_time,
        x_bins=x_bins,
        y_bins=y_bins,
        n_spikes_used=len(x_at_spike),
    )


def spatial_information_bits_per_spike(rate_map: RateMap) -> float:
    """Skaggs et al. (1993) spatial information, in bits/spike.

    I = sum_i p_i * (r_i / r_bar) * log2(r_i / r_bar)
    where p_i is occupancy probability of bin i, r_i is the bin's firing rate,
    and r_bar is the overall mean firing rate (occupancy-weighted).
    """
    occ = rate_map.occupancy_time
    rate = rate_map.firing_rate_map
    visited = occ > 0
    if not np.any(visited):
        return np.nan

    p = occ[visited] / occ[visited].sum()
    r = rate[visited]
    r_bar = np.sum(p * r)
    if r_bar <= 0:
        return np.nan

    with np.errstate(divide="ignore", invalid="ignore"):
        term = p * r * np.log2(r / r_bar)
    term[~np.isfinite(term)] = 0.0  # bins with r_i = 0 contribute 0, not NaN
    return float(np.sum(term) / r_bar)


def field_count(
    rate_map: RateMap,
    threshold_frac: float = 0.2,
    min_bin_count: int = 3,
    smoothing_sigma_bins: float = 1.5,
) -> int:
    """Count spatially contiguous regions exceeding threshold_frac * peak rate.

    Simple 4-connected flood fill; matches the exercise 6 hint in notebook 02a
    (bins where firing rate exceeds 20% of the maximum). Smooths the map first
    (Muller & Kubie 1987 convention) — at 2 cm bins the raw map is speckly
    enough that unsmoothed thresholding fragments one field into many.
    """
    smoothed = gaussian_filter(rate_map.firing_rate_map, sigma=smoothing_sigma_bins)
    peak = np.nanmax(smoothed)
    if not np.isfinite(peak) or peak <= 0:
        return 0
    above = smoothed >= (threshold_frac * peak)

    visited = np.zeros_like(above, dtype=bool)
    n_fields = 0
    nx, ny = above.shape
    for i in range(nx):
        for j in range(ny):
            if above[i, j] and not visited[i, j]:
                stack = [(i, j)]
                size = 0
                while stack:
                    ci, cj = stack.pop()
                    if ci < 0 or ci >= nx or cj < 0 or cj >= ny:
                        continue
                    if visited[ci, cj] or not above[ci, cj]:
                        continue
                    visited[ci, cj] = True
                    size += 1
                    stack.extend([(ci + 1, cj), (ci - 1, cj), (ci, cj + 1), (ci, cj - 1)])
                if size >= min_bin_count:
                    n_fields += 1
    return n_fields
