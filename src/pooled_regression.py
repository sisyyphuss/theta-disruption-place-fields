"""Pool all subjects' session-level place-field stats and test whether cumulative
prior theta-disruption exposure (cooling) predicts spatial coding quality,
controlling for subject identity and elapsed time, using a permutation null
rather than parametric p-values (n=15 sessions across 4 animals).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

REPO_DIR = Path(__file__).resolve().parent.parent
CACHE_DIR = REPO_DIR / "results" / "cache"
SUBJECTS = ["MS10", "MS12", "MS13", "MS21", "MS22"]
N_PERMUTATIONS = 5000
RNG_SEED = 0


def load_pooled() -> pd.DataFrame:
    frames = []
    for subject in SUBJECTS:
        path = CACHE_DIR / f"{subject}_sessions.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if len(df) > 0:
            frames.append(df)
    pooled = pd.concat(frames, ignore_index=True)
    pooled["subject"] = pooled["subject"].astype("category")
    return pooled


def fit_ols(
    df: pd.DataFrame,
    outcome: str,
    predictor: str = "cumulative_cooling_seconds_prior",
):
    """OLS with subject fixed effects: outcome ~ predictor + C(subject)."""
    model = smf.ols(f"{outcome} ~ {predictor} + C(subject)", data=df).fit()
    return model.params[predictor], model.pvalues[predictor], model


def permutation_test(
    df: pd.DataFrame,
    outcome: str,
    predictor: str = "cumulative_cooling_seconds_prior",
    rng: np.random.Generator = None,
) -> dict:
    if rng is None:
        rng = np.random.default_rng(RNG_SEED)
    observed_coef, parametric_p, observed_model = fit_ols(df, outcome, predictor)

    null_coefs = np.empty(N_PERMUTATIONS)
    for i in range(N_PERMUTATIONS):
        permuted = df.copy()
        # Shuffle predictor within each subject's own sessions
        permuted[predictor] = (
            permuted.groupby("subject", observed=True)[predictor]
            .transform(lambda s: rng.permutation(s.values))
        )
        try:
            coef, _, _ = fit_ols(permuted, outcome, predictor)
        except Exception:
            coef = np.nan
        null_coefs[i] = coef

    null_coefs = null_coefs[np.isfinite(null_coefs)]
    p_value = float(np.mean(np.abs(null_coefs) >= np.abs(observed_coef)))
    return dict(
        outcome=outcome,
        predictor=predictor,
        observed_coef=observed_coef,
        parametric_p=parametric_p,
        p_two_sided=p_value,
        n_permutations_used=len(null_coefs),
        null_mean=float(np.mean(null_coefs)),
        null_std=float(np.std(null_coefs)),
        model_summary=observed_model.summary().as_text(),
    )


if __name__ == "__main__":
    df = load_pooled()
    print(f"Pooled dataset: {len(df)} sessions across {df['subject'].nunique()} subjects")
    print(df.groupby("subject", observed=True).size().rename("n_sessions"))
    print()
    print(df[["subject", "date", "elapsed_days", "cumulative_cooling_seconds_prior",
              "spatial_info_median", "field_count_median", "peak_rate_median", "n_units_used"]].to_string(index=False))

    rng = np.random.default_rng(RNG_SEED)
    for outcome in ["spatial_info_median", "field_count_median", "peak_rate_median"]:
        print(f"\n{'=' * 70}\nOutcome: {outcome}")
        result = permutation_test(df, outcome, "cumulative_cooling_seconds_prior", rng)
        print(f"  observed coefficient: {result['observed_coef']:.6f}")
        print(f"  permutation p-value:  {result['p_two_sided']:.4f}")
        print(f"  parametric OLS p-val: {result['parametric_p']:.4f}")
