"""Generate manuscript figures from the cached pooled session-level dataset."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pooled_regression import load_pooled

REPO_DIR = Path(__file__).resolve().parent.parent
FIG_DIR = REPO_DIR / "results" / "figures"

attempted = {"MS10": 9, "MS12": 8, "MS13": 5, "MS21": 15, "MS22": 13}


def generate_figures():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = load_pooled()
    subjects = sorted(df["subject"].cat.categories)
    colors = dict(zip(subjects, plt.cm.tab10.colors))

    # Figure 1: data yield per subject (attempted vs usable sessions)
    usable = df.groupby("subject", observed=True).size().reindex(attempted.keys(), fill_value=0)

    fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
    x = np.arange(len(attempted))
    ax.bar(x - 0.2, list(attempted.values()), width=0.4, label="Sessions attempted", color="lightgray")
    ax.bar(x + 0.2, usable.values, width=0.4, label="Sessions usable", color="steelblue")
    ax.set(xticks=x, xticklabels=list(attempted.keys()), ylabel="Number of sessions",
           title="Session yield per subject, Dandiset 000059")
    ax.legend()
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(FIG_DIR / "fig1_session_yield.png", dpi=150)
    plt.close(fig)

    # Figure 2: outcome measures vs elapsed days and cumulative cooling dose
    outcomes = [
        ("spatial_info_median", "Spatial information (bits/spike)"),
        ("field_count_median", "Field count"),
        ("peak_rate_median", "Peak rate (Hz)"),
    ]
    predictors = [
        ("elapsed_days", "Elapsed days since first session"),
        ("cumulative_cooling_seconds_prior", "Cumulative prior cooling exposure (s)"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(10, 10), layout="constrained")
    for row, (outcome, outcome_label) in enumerate(outcomes):
        for col, (predictor, predictor_label) in enumerate(predictors):
            ax = axes[row, col]
            for subject in subjects:
                sub_df = df[df["subject"] == subject].sort_values(predictor)
                ax.plot(sub_df[predictor], sub_df[outcome], "o-", color=colors[subject],
                         label=subject, markersize=6, linewidth=1)
            ax.set(xlabel=predictor_label if row == 2 else "",
                   ylabel=outcome_label if col == 0 else "")
            ax.spines[["top", "right"]].set_visible(False)
    axes[0, 1].legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8, title="Subject")
    fig.suptitle("Session-level place-field statistics vs. elapsed time and cumulative cooling exposure")
    fig.savefig(FIG_DIR / "fig2_trends.png", dpi=150)
    plt.close(fig)

    print("Figures written to", FIG_DIR)
    for f in sorted(FIG_DIR.iterdir()):
        print(" ", f.name)


if __name__ == "__main__":
    generate_figures()
