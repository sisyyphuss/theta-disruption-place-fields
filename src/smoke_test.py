import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from data_loading import list_processed_sessions, open_session, get_good_units, get_position, session_datetime
from cooling_dose import summarize_cooling, get_trials, mask_to_baseline
from place_fields import compute_rate_map, spatial_information_bits_per_spike, field_count

paths = list_processed_sessions("MS21")
print(f"Found {len(paths)} sessions for MS21")

for path in paths[:4]:
    print("\n" + "=" * 70)
    print(path)
    session = open_session(path)
    print("session_start_time:", session_datetime(session))

    trials_df = get_trials(session)
    print("cooling state values:", trials_df["cooling state"].unique())

    cooling = summarize_cooling(session)
    print("cooling summary:", cooling)

    x, y, pos_t = get_position(session)
    units_df, good_idx = get_good_units(session)
    if units_df is None:
        print("  no units table in this session -- skipping")
        session.close()
        continue
    print(f"n good units: {len(good_idx)}")

    n_ok = 0
    for idx in good_idx[:5]:
        spike_times = session.nwbfile.units["spike_times"][idx]
        b_spikes, b_x, b_y, b_t = mask_to_baseline(
            spike_times, x, y, pos_t, cooling.baseline_start, cooling.baseline_stop
        )
        if len(b_t) < 10:
            continue
        rm = compute_rate_map(b_spikes, b_x, b_y, b_t)
        si = spatial_information_bits_per_spike(rm)
        nf = field_count(rm)
        print(f"  unit {idx}: n_spikes(baseline)={rm.n_spikes_used}, spatial_info={si:.3f} bits/spike, n_fields={nf}")
        n_ok += 1
    print(f"Successfully processed {n_ok}/5 sampled units")

    session.close()

print("\nSMOKE TEST PASSED")