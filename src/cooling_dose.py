"""Cooling-epoch extraction: per-session baseline-trial windows and cooling dose.

Uses the `trials` table's `cooling state` column (values observed so far:
"Pre-Cooling"; "Cooling"/"Post-Cooling" expected later in the same sessions)
plus the behavior module's `Temperature` series for a continuous dose measure.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from data_loading import SessionHandle


@dataclass
class CoolingSummary:
    baseline_start: float
    baseline_stop: float
    cooling_trial_seconds: float  # total time-in-trial spent in "Cooling" state
    mean_temperature_drop: float  # baseline mean temp minus cooling-epoch mean temp, degC (nan if unavailable)


def get_trials(session: SessionHandle):
    return session.nwbfile.trials.to_dataframe()


def baseline_window(trials_df) -> tuple[float, float]:
    """Time span covered by contiguous 'Pre-Cooling' trials at the start of the session."""
    pre = trials_df[trials_df["cooling state"] == "Pre-Cooling"]
    if len(pre) == 0:
        raise ValueError("No Pre-Cooling trials found in this session")
    return float(pre["start_time"].min()), float(pre["stop_time"].max())


COOLING_ON_LABELS = {"Cooling on", "Cooling off"}  # the probe-cooled epoch, on->off ramp


def cooling_trial_seconds(trials_df) -> float:
    cooling = trials_df[trials_df["cooling state"].isin(COOLING_ON_LABELS)]
    return float((cooling["stop_time"] - cooling["start_time"]).sum())


def temperature_drop(session: SessionHandle, trials_df) -> float:
    behavior_module = session.nwbfile.processing.get("behavior")
    if behavior_module is None or "Temperature" not in behavior_module.data_interfaces:
        return float("nan")
    temp_series = behavior_module.data_interfaces["Temperature"]
    temp_data = temp_series.data[:]
    if temp_series.timestamps is not None:
        temp_time = temp_series.timestamps[:]
    else:
        temp_time = np.arange(len(temp_data)) / temp_series.rate

    baseline_start, baseline_stop = baseline_window(trials_df)
    baseline_mask = (temp_time >= baseline_start) & (temp_time <= baseline_stop)

    cooling_rows = trials_df[trials_df["cooling state"].isin(COOLING_ON_LABELS)]
    if len(cooling_rows) == 0 or not np.any(baseline_mask):
        return float("nan")
    cooling_start, cooling_stop = float(cooling_rows["start_time"].min()), float(cooling_rows["stop_time"].max())
    cooling_mask = (temp_time >= cooling_start) & (temp_time <= cooling_stop)
    if not np.any(cooling_mask):
        return float("nan")

    return float(np.nanmean(temp_data[baseline_mask]) - np.nanmean(temp_data[cooling_mask]))


def summarize_cooling(session: SessionHandle) -> CoolingSummary:
    trials_df = get_trials(session)
    b_start, b_stop = baseline_window(trials_df)
    cool_secs = cooling_trial_seconds(trials_df)
    temp_drop = temperature_drop(session, trials_df)
    return CoolingSummary(
        baseline_start=b_start,
        baseline_stop=b_stop,
        cooling_trial_seconds=cool_secs,
        mean_temperature_drop=temp_drop,
    )


def mask_to_baseline(
    spike_times: NDArray[np.floating],
    x_position: NDArray[np.floating],
    y_position: NDArray[np.floating],
    position_timestamps: NDArray[np.floating],
    baseline_start: float,
    baseline_stop: float,
):
    """Restrict spikes and position samples to the Pre-Cooling window."""
    spike_mask = (spike_times >= baseline_start) & (spike_times <= baseline_stop)
    pos_mask = (position_timestamps >= baseline_start) & (position_timestamps <= baseline_stop)
    return (
        spike_times[spike_mask],
        x_position[pos_mask],
        y_position[pos_mask],
        position_timestamps[pos_mask],
    )
