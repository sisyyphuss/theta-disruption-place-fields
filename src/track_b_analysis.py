"""Track B: session-level place-field statistics vs. cumulative theta-disruption dose.

For one subject, streams every processed session, computes baseline (Pre-Cooling)
place fields for all good units using the notebook 1/2a-derived pipeline, reduces
each session to summary statistics, and regresses those against cumulative prior
cooling exposure and elapsed days since the first session. Significance is
assessed with a permutation null rather than a parametric p-value, since we
only have ~15 sessions per animal.

Results are cached to CSV per subject so re-runs don't re-stream from DANDI.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from data_loading import list_processed_sessions, open_session, get_good_units, get_position, session_datetime
from cooling_dose import summarize_cooling, get_trials, mask_to_baseline
from place_fields import compute_rate_map, spatial_information_bits_per_spike, field_count

CACHE_DIR = Path(__file__).parent / "cache"
CACHE_DIR.mkdir(exist_ok=True)


def process_subject(subject: str, min_baseline_spikes: int = 30) -> pd.DataFrame:
    cache_path = CACHE_DIR / f"{subject}_sessions.csv"
    if cache_path.exists():
        return pd.read_csv(cache_path, parse_dates=["session_date"])

    paths = list_processed_sessions(subject)
    rows = []
    for path in paths:
        try:
            session = open_session(path)
        except Exception as exc:  # noqa: BLE001 - report and skip malformed assets
            print(f"  [skip] {path}: failed to open ({exc})")
            continue

        try:
            units_df, good_idx = get_good_units(session)
            if units_df is None or len(good_idx) == 0:
                print(f"  [skip] {path}: no usable units")
                continue

            trials_df = get_trials(session)
            if "cooling state" not in trials_df.columns:
                print(f"  [skip] {path}: trials table has no 'cooling state' column "
                      f"(columns: {list(trials_df.columns)})")
                continue
            try:
                cooling = summarize_cooling(session)
            except ValueError as exc:
                print(f"  [skip] {path}: {exc}")
                continue
            x, y, pos_t = get_position(session)

            spatial_infos, n_fields_list, peak_rates = [], [], []
            for idx in good_idx:
                spike_times = session.nwbfile.units["spike_times"][idx]
                b_spikes, b_x, b_y, b_t = mask_to_baseline(
                    spike_times, x, y, pos_t, cooling.baseline_start, cooling.baseline_stop
                )
                if len(b_spikes) < min_baseline_spikes:
                    continue
                rm = compute_rate_map(b_spikes, b_x, b_y, b_t)
                si = spatial_information_bits_per_spike(rm)
                if not np.isfinite(si):
                    continue
                spatial_infos.append(si)
                n_fields_list.append(field_count(rm))
                peak_rates.append(float(np.nanmax(rm.firing_rate_map)))

            if len(spatial_infos) == 0:
                print(f"  [skip] {path}: no units passed spike-count threshold")
                continue

            rows.append(
                dict(
                    subject=subject,
                    path=path,
                    session_date=session_datetime(session).date().isoformat(),
                    session_hour=session_datetime(session).hour
                    + session_datetime(session).minute / 60,
                    n_units_used=len(spatial_infos),
                    n_trials=len(trials_df),
                    error_rate=float(trials_df["error"].mean()),
                    cooling_trial_seconds=cooling.cooling_trial_seconds,
                    mean_temperature_drop=cooling.mean_temperature_drop,
                    spatial_info_median=float(np.median(spatial_infos)),
                    field_count_median=float(np.median(n_fields_list)),
                    peak_rate_median=float(np.median(peak_rates)),
                )
            )
            print(f"  [ok] {path}: {len(spatial_infos)} units, "
                  f"SI_median={np.median(spatial_infos):.3f}")
        except Exception as exc:  # noqa: BLE001 - one bad session shouldn't kill the batch
            print(f"  [skip] {path}: unexpected error ({type(exc).__name__}: {exc})")
        finally:
            session.close()

    if len(rows) == 0:
        print(f"  no usable sessions at all for {subject}")
        df = pd.DataFrame(rows)
        df.to_csv(cache_path, index=False)
        return df

    df = pd.DataFrame(rows).sort_values("session_date").reset_index(drop=True)
    df["session_date"] = pd.to_datetime(df["session_date"])
    df["elapsed_days"] = (df["session_date"] - df["session_date"].min()).dt.days
    # Cumulative cooling dose accrued *before* this session (prior sessions only) --
    # this is the "chronic exposure so far" predictor, distinct from this session's own dose.
    df["cumulative_cooling_seconds_prior"] = df["cooling_trial_seconds"].cumsum() - df["cooling_trial_seconds"]

    df.to_csv(cache_path, index=False)
    return df


if __name__ == "__main__":
    subject = sys.argv[1] if len(sys.argv) > 1 else "MS21"
    print(f"Processing subject {subject}...")
    df = process_subject(subject)
    print(f"\n{len(df)} usable sessions")
    print(df[["session_date", "elapsed_days", "cooling_trial_seconds",
              "cumulative_cooling_seconds_prior", "spatial_info_median",
              "field_count_median", "n_units_used"]].to_string(index=False))
