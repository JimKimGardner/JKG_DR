#!/usr/bin/env python3

from __future__ import annotations

import math
import re
import warnings as pywarnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from matplotlib.ticker import FormatStrFormatter
from patsy import build_design_matrices
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import braycurtis
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

INPUT_CANDIDATES = [
    "data/metaphlan_sgb.tsv",
    "data/metaphlan.tsv",
    "data/taxa.tsv",
    "metaphlan_sgb.tsv",
    "metaphlan.tsv",
    "data/processed/taxa/TAXA_matrix.tsv",
]
OUTPUT_DIR = Path("results")
TIMEPOINT_ORDER = ["w1", "w2", "w4", "w6", "w8", "w10", "w12", "w16", "w52"]
ALLOWED_WEEKS = {0, 1, 2, 4, 6, 8, 10, 12, 16, 52}
WEEK_MAP = {f"w{week}": week for week in sorted(ALLOWED_WEEKS)}
CORE_PREVALENCE_THRESHOLD = 0.70
DETECTION_THRESHOLD = 1e-4
GLOBAL_PREVALENCE_FILTER = 0.10
PSEUDOCOUNT = 1e-6
RECOVERY_THRESHOLD_BRAY = 0.10
DEFAULT_N_CLUSTERS = 3


@dataclass
class FitResult:
    metric: str
    compartment: str
    model_kind: str
    method: str
    reference_timepoint: str
    fitted: bool
    coefficients: pd.DataFrame
    estimated_means: pd.DataFrame
    pairwise: pd.DataFrame
    raw_result: object
    warnings: List[str]


def resolve_input_path() -> Path:
    for candidate in INPUT_CANDIDATES:
        path = Path(candidate)
        if path.exists():
            return path
    raise FileNotFoundError(f"No MetaPhlAn table found. Tried: {', '.join(INPUT_CANDIDATES)}")


def is_sample_like(value: str) -> bool:
    return bool(re.search(r"w0*\d+", str(value), flags=re.IGNORECASE))


def trim_delimiters(value: str) -> str:
    value = re.sub(r"^[\W_]+", "", value)
    value = re.sub(r"[\W_]+$", "", value)
    value = re.sub(r"\s+", " ", value)
    return value.strip()


def load_metaphlan_table(path: Path) -> Tuple[pd.DataFrame, str]:
    df = pd.read_csv(path, sep="\t")
    if df.shape[1] < 2:
        raise ValueError("Input table must contain at least one identifier column and one sample/data column.")

    id_values = df.iloc[:, 0].astype(str)
    data_df = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    col_sample_hits = sum(is_sample_like(col) for col in data_df.columns)
    row_sample_hits = sum(is_sample_like(val) for val in id_values)

    if col_sample_hits >= row_sample_hits:
      feature_mat = data_df.copy()
      feature_mat.index = id_values
      orientation = "features_in_rows_samples_in_columns"
    else:
      feature_mat = data_df.copy()
      feature_mat.index = id_values
      feature_mat = feature_mat.T
      orientation = "features_in_columns_samples_in_rows"

    sgb_mask = feature_mat.index.to_series().str.contains(r"(?:^|\|)t__SGB", regex=True)
    if sgb_mask.any():
        feature_mat = feature_mat.loc[sgb_mask]

    feature_mat = feature_mat.loc[(feature_mat != 0).any(axis=1)]
    if feature_mat.empty:
        raise ValueError("No non-zero SGB-level features remained after loading the MetaPhlAn table.")

    return feature_mat, orientation


def parse_sample_name(sample_id: str) -> Dict[str, object]:
    sample_id = str(sample_id)
    match = re.search(r"w0*\d+", sample_id, flags=re.IGNORECASE)
    if not match:
        return {
            "sample_id": sample_id,
            "patient_id": None,
            "timepoint": None,
            "weeks_post_return": np.nan,
            "parse_status": "unparsed",
            "exclusion_reason": "timepoint_not_recognized",
        }

    timepoint_raw = match.group(0).lower()
    digits = re.sub(r"^w0*", "", timepoint_raw)
    weeks = 0 if digits == "" else int(digits)

    patient_prefix = trim_delimiters(sample_id[: match.start()])
    patient_suffix = trim_delimiters(sample_id[match.end() :])
    patient_id = (patient_prefix or patient_suffix).lower()
    if not patient_id:
        return {
            "sample_id": sample_id,
            "patient_id": None,
            "timepoint": None,
            "weeks_post_return": np.nan,
            "parse_status": "unparsed",
            "exclusion_reason": "patient_or_timepoint_not_reliably_parsed",
        }

    timepoint = f"w{weeks}"
    if weeks not in ALLOWED_WEEKS:
        return {
            "sample_id": sample_id,
            "patient_id": patient_id,
            "timepoint": timepoint,
            "weeks_post_return": weeks,
            "parse_status": "unsupported_timepoint",
            "exclusion_reason": "timepoint_not_in_analysis_set",
        }

    return {
        "sample_id": sample_id,
        "patient_id": patient_id,
        "timepoint": timepoint,
        "weeks_post_return": weeks,
        "parse_status": "ok",
        "exclusion_reason": None,
    }


def build_sample_metadata(sample_ids: Sequence[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sample_meta = pd.DataFrame(parse_sample_name(sample_id) for sample_id in sample_ids)
    sample_meta["sample_order"] = np.arange(sample_meta.shape[0])

    ok_samples = sample_meta.loc[sample_meta["parse_status"] == "ok"].copy()
    ok_samples = ok_samples.sort_values(["patient_id", "weeks_post_return", "sample_order"])

    patient_summary = (
        ok_samples.groupby("patient_id")
        .agg(
            n_w0=("weeks_post_return", lambda x: int((x == 0).sum())),
            n_post=("weeks_post_return", lambda x: int((x > 0).sum())),
        )
        .reset_index()
    )
    baseline_lookup = (
        ok_samples.loc[ok_samples["weeks_post_return"] == 0]
        .sort_values(["patient_id", "sample_order"])
        .groupby("patient_id", as_index=False)
        .agg(
            baseline_sample_id=("sample_id", "first"),
            baseline_candidates=("sample_id", lambda x: "; ".join(x)),
        )
    )

    patient_summary = patient_summary.merge(baseline_lookup, on="patient_id", how="left")
    patient_summary["eligible"] = (patient_summary["n_w0"] > 0) & (patient_summary["n_post"] > 0)
    patient_summary["multiple_w0"] = patient_summary["n_w0"] > 1

    sample_meta = sample_meta.merge(
        patient_summary[["patient_id", "n_w0", "n_post", "baseline_sample_id", "baseline_candidates", "eligible", "multiple_w0"]],
        on="patient_id",
        how="left",
    )
    sample_meta["inclusion_status"] = np.select(
        [
            sample_meta["parse_status"] != "ok",
            sample_meta["n_w0"].fillna(0) == 0,
            sample_meta["n_post"].fillna(0) == 0,
            sample_meta["eligible"].fillna(False),
        ],
        [
            sample_meta["exclusion_reason"],
            "excluded_no_w0",
            "excluded_no_post",
            "included",
        ],
        default="excluded_other",
    )
    return sample_meta, patient_summary


def make_relative_abundance(feature_df: pd.DataFrame) -> pd.DataFrame:
    totals = feature_df.sum(axis=0)
    if (totals <= 0).any():
        bad = totals[totals <= 0].index.tolist()
        raise ValueError(f"Samples with zero total abundance encountered: {', '.join(bad)}")
    return feature_df.divide(totals, axis=1)


def apply_global_prevalence_filter(
    relative_df: pd.DataFrame,
    min_prevalence: float = GLOBAL_PREVALENCE_FILTER,
    detection_threshold: float = DETECTION_THRESHOLD,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    prevalence = (relative_df >= detection_threshold).mean(axis=1)
    prevalence_table = pd.DataFrame(
        {
            "feature_id": relative_df.index,
            "global_prevalence": prevalence.values,
            "n_present_samples": (relative_df >= detection_threshold).sum(axis=1).values,
        }
    )
    keep = prevalence >= min_prevalence
    filtered = relative_df.loc[keep].copy()
    filtered = filtered.loc[(filtered > 0).any(axis=1)]
    prevalence_table = prevalence_table.loc[prevalence_table["feature_id"].isin(filtered.index)].reset_index(drop=True)
    return filtered, prevalence_table


def define_core_accessory(
    relative_df: pd.DataFrame,
    sample_meta: pd.DataFrame,
    prevalence_threshold: float = CORE_PREVALENCE_THRESHOLD,
    detection_threshold: float = DETECTION_THRESHOLD,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    baseline_samples = (
        sample_meta.loc[(sample_meta["inclusion_status"] == "included") & (sample_meta["weeks_post_return"] == 0)]
        .sort_values(["patient_id", "sample_order"])
        .drop_duplicates("patient_id")
    )
    baseline_df = relative_df.loc[:, baseline_samples["sample_id"]]

    feature_table = pd.DataFrame(
        {
            "feature_id": relative_df.index,
            "prevalence_w0": (baseline_df >= detection_threshold).mean(axis=1).values,
            "mean_abundance_w0": baseline_df.mean(axis=1).values,
            "median_abundance_w0": baseline_df.median(axis=1).values,
        }
    )
    feature_table["is_core"] = feature_table["prevalence_w0"] >= prevalence_threshold
    feature_table["compartment"] = np.where(feature_table["is_core"], "core", "accessory")
    feature_table["core_definition_used"] = "cohort_core_from_w0_only"

    core_features = feature_table.loc[feature_table["is_core"], "feature_id"].tolist()
    accessory_features = feature_table.loc[~feature_table["is_core"], "feature_id"].tolist()
    return feature_table, core_features, accessory_features


def compute_clr(feature_df: pd.DataFrame, pseudocount: float = PSEUDOCOUNT) -> pd.DataFrame:
    adjusted = feature_df + pseudocount
    adjusted = adjusted.divide(adjusted.sum(axis=0), axis=1)
    logged = np.log(adjusted)
    return logged.subtract(logged.mean(axis=0), axis=1)


def compute_distance_long_for_compartment(
    relative_df: pd.DataFrame,
    sample_meta: pd.DataFrame,
    feature_ids: Sequence[str],
    compartment: str,
    pseudocount: float = PSEUDOCOUNT,
) -> pd.DataFrame:
    post_meta = (
        sample_meta.loc[(sample_meta["inclusion_status"] == "included") & (sample_meta["weeks_post_return"] > 0)]
        .sort_values(["patient_id", "weeks_post_return", "sample_order"])
        .copy()
    )
    if len(feature_ids) == 0:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "patient_id",
                "baseline_sample_id",
                "timepoint",
                "weeks_post_return",
                "compartment",
                "distance_type",
                "distance_to_baseline",
                "n_features",
            ]
        )

    subset_df = relative_df.loc[relative_df.index.intersection(feature_ids)].copy()
    subset_df = subset_df.loc[(subset_df > 0).any(axis=1)]
    if subset_df.empty:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "patient_id",
                "baseline_sample_id",
                "timepoint",
                "weeks_post_return",
                "compartment",
                "distance_type",
                "distance_to_baseline",
                "n_features",
            ]
        )

    subset_df = subset_df.divide(subset_df.sum(axis=0), axis=1).fillna(0.0)
    clr_df = compute_clr(subset_df, pseudocount=pseudocount)

    rows: List[Dict[str, object]] = []
    for row in post_meta.itertuples(index=False):
        post_vec = subset_df[row.sample_id].values
        baseline_vec = subset_df[row.baseline_sample_id].values
        bray = braycurtis(post_vec, baseline_vec)

        post_clr = clr_df[row.sample_id].values
        baseline_clr = clr_df[row.baseline_sample_id].values
        aitchison = float(np.linalg.norm(post_clr - baseline_clr))

        rows.extend(
            [
                {
                    "sample_id": row.sample_id,
                    "patient_id": row.patient_id,
                    "baseline_sample_id": row.baseline_sample_id,
                    "timepoint": row.timepoint,
                    "weeks_post_return": row.weeks_post_return,
                    "compartment": compartment,
                    "distance_type": "bray",
                    "distance_to_baseline": float(bray),
                    "n_features": int(subset_df.shape[0]),
                },
                {
                    "sample_id": row.sample_id,
                    "patient_id": row.patient_id,
                    "baseline_sample_id": row.baseline_sample_id,
                    "timepoint": row.timepoint,
                    "weeks_post_return": row.weeks_post_return,
                    "compartment": compartment,
                    "distance_type": "aitchison",
                    "distance_to_baseline": aitchison,
                    "n_features": int(subset_df.shape[0]),
                },
            ]
        )

    return pd.DataFrame(rows)


def make_distance_outputs(
    relative_df: pd.DataFrame,
    sample_meta: pd.DataFrame,
    core_features: Sequence[str],
    accessory_features: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    long_frames = [
        compute_distance_long_for_compartment(relative_df, sample_meta, relative_df.index.tolist(), "all"),
        compute_distance_long_for_compartment(relative_df, sample_meta, core_features, "core"),
        compute_distance_long_for_compartment(relative_df, sample_meta, accessory_features, "accessory"),
    ]
    distance_long = pd.concat(long_frames, ignore_index=True)
    distance_long["timepoint"] = pd.Categorical(distance_long["timepoint"], categories=TIMEPOINT_ORDER, ordered=True)

    distance_wide = (
        distance_long.pivot_table(
            index=["sample_id", "patient_id", "baseline_sample_id", "timepoint", "weeks_post_return"],
            columns=["compartment", "distance_type"],
            values="distance_to_baseline",
            observed=False,
        )
        .reset_index()
    )
    distance_wide.columns = [
        "_".join(col).strip("_") if isinstance(col, tuple) else col for col in distance_wide.columns.to_flat_index()
    ]
    return distance_long, distance_wide


def add_time_factor(data: pd.DataFrame) -> pd.DataFrame:
    present = [tp for tp in TIMEPOINT_ORDER if tp in set(data["timepoint"].astype(str))]
    data = data.copy()
    data["timepoint"] = pd.Categorical(data["timepoint"].astype(str), categories=present, ordered=True)
    data["time_factor"] = data["timepoint"]
    return data


def design_row_from_formula(result, row_df: pd.DataFrame) -> np.ndarray:
    design_info = result.model.data.design_info
    return np.asarray(build_design_matrices([design_info], row_df, return_type="dataframe")[0])[0]


def get_fixed_effects_and_cov(result) -> Tuple[pd.Series, pd.DataFrame]:
    if hasattr(result, "fe_params"):
        params = result.fe_params
        cov = result.cov_params().loc[params.index, params.index]
        return params, cov
    params = result.params
    cov = result.cov_params().loc[params.index, params.index]
    return params, cov


def fit_mixedlm_with_fallback(formula: str, data: pd.DataFrame, groups: pd.Series, vc_formula: Optional[Dict[str, str]] = None):
    warnings: List[str] = []
    try:
        model = smf.mixedlm(formula, data=data, groups=groups, vc_formula=vc_formula, re_formula="1")
        with pywarnings.catch_warnings(record=True) as captured:
            pywarnings.simplefilter("always")
            result = model.fit(reml=False, method="lbfgs", maxiter=200, disp=False)
        warnings.extend([str(w.message) for w in captured])
        return result, "MixedLM", warnings
    except Exception as exc:
        warnings.append(f"MixedLM fit failed ({exc}). Falling back to OLS with cluster-robust SE by patient.")
        with pywarnings.catch_warnings(record=True) as captured:
            pywarnings.simplefilter("always")
            ols = smf.ols(formula, data=data).fit(cov_type="cluster", cov_kwds={"groups": data["patient_id"]})
        warnings.extend([str(w.message) for w in captured])
        return ols, "OLS_clustered_by_patient", warnings


def compute_estimated_means_and_pairwise(
    result,
    data: pd.DataFrame,
    reference_timepoint: str,
    compartment: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    timepoints = [tp for tp in TIMEPOINT_ORDER if tp in set(data["timepoint"].astype(str))]
    params, cov = get_fixed_effects_and_cov(result)

    means_rows = []
    row_vectors: Dict[str, np.ndarray] = {}
    for timepoint in timepoints:
        row = {"timepoint": timepoint}
        if compartment is not None:
            row["compartment"] = compartment
        row_df = pd.DataFrame([row])
        x = design_row_from_formula(result, row_df)
        estimate = float(np.dot(x, params.values))
        se = float(np.sqrt(np.dot(x, np.dot(cov.values, x))))
        means_rows.append(
            {
                "timepoint": timepoint,
                "estimate": estimate,
                "std_error": se,
                "conf_low": estimate - 1.96 * se,
                "conf_high": estimate + 1.96 * se,
            }
        )
        row_vectors[timepoint] = x

    means_df = pd.DataFrame(means_rows)

    contrast_rows = []
    for idx, tp1 in enumerate(timepoints):
        for tp2 in timepoints[idx + 1 :]:
            diff = row_vectors[tp1] - row_vectors[tp2]
            estimate = float(np.dot(diff, params.values))
            se = float(np.sqrt(np.dot(diff, np.dot(cov.values, diff))))
            z_value = estimate / se if se > 0 else np.nan
            p_value = 2 * (1 - norm.cdf(abs(z_value))) if np.isfinite(z_value) else np.nan
            contrast_rows.append(
                {
                    "contrast": f"{tp1} - {tp2}",
                    "timepoint_1": tp1,
                    "timepoint_2": tp2,
                    "estimate": estimate,
                    "std_error": se,
                    "z_value": z_value,
                    "p_value": p_value,
                }
            )

    pairwise_df = pd.DataFrame(contrast_rows)
    if not pairwise_df.empty:
        valid_mask = pairwise_df["p_value"].notna()
        pairwise_df["p_value_bh"] = np.nan
        if valid_mask.any():
            pairwise_df.loc[valid_mask, "p_value_bh"] = multipletests(pairwise_df.loc[valid_mask, "p_value"], method="fdr_bh")[1]
    else:
        pairwise_df = pd.DataFrame(columns=["contrast", "timepoint_1", "timepoint_2", "estimate", "std_error", "z_value", "p_value", "p_value_bh"])

    return means_df, pairwise_df


def coefficients_to_frame(result, metric: str, compartment: str, model_kind: str, method: str) -> pd.DataFrame:
    if hasattr(result, "pvalues"):
        pvalues = result.pvalues
        params = result.params
        conf = result.conf_int()
        se = getattr(result, "bse", pd.Series(index=params.index, dtype=float))
    else:
        pvalues = pd.Series(dtype=float)
        params = pd.Series(dtype=float)
        conf = pd.DataFrame(columns=[0, 1])
        se = pd.Series(dtype=float)

    coef_df = pd.DataFrame(
        {
            "term": params.index,
            "estimate": params.values,
            "std_error": se.reindex(params.index).values,
            "conf_low": conf.reindex(params.index)[0].values,
            "conf_high": conf.reindex(params.index)[1].values,
            "p_value": pvalues.reindex(params.index).values,
            "metric": metric,
            "compartment": compartment,
            "model_kind": model_kind,
            "method": method,
        }
    )
    return coef_df


def fit_time_factor_model(data: pd.DataFrame, compartment: str, metric: str) -> FitResult:
    subset = data.loc[(data["compartment"] == compartment) & (data["distance_type"] == metric)].copy()
    subset = subset.dropna(subset=["distance_to_baseline", "patient_id", "timepoint"])
    subset = add_time_factor(subset)
    available = subset["timepoint"].cat.categories.tolist()
    reference = "w1" if "w1" in available else available[0]
    formula = f"distance_to_baseline ~ C(timepoint, Treatment(reference='{reference}'))"
    result, method, warnings = fit_mixedlm_with_fallback(formula, subset, subset["patient_id"])
    coef_df = coefficients_to_frame(result, metric, compartment, "time_factor_separate", method)
    means_df, pairwise_df = compute_estimated_means_and_pairwise(result, subset, reference)
    means_df["metric"] = metric
    means_df["compartment"] = compartment
    means_df["method"] = method
    pairwise_df["metric"] = metric
    pairwise_df["compartment"] = compartment
    pairwise_df["method"] = method
    return FitResult(metric, compartment, "time_factor_separate", method, reference, True, coef_df, means_df, pairwise_df, result, warnings)


def fit_joint_interaction_model(data: pd.DataFrame, metric: str) -> FitResult:
    subset = data.loc[(data["compartment"].isin(["core", "accessory"])) & (data["distance_type"] == metric)].copy()
    subset = subset.dropna(subset=["distance_to_baseline", "patient_id", "timepoint", "sample_id"])
    subset = add_time_factor(subset)
    available = subset["timepoint"].cat.categories.tolist()
    reference = "w1" if "w1" in available else available[0]
    subset["compartment"] = pd.Categorical(subset["compartment"], categories=["core", "accessory"], ordered=True)
    formula = (
        f"distance_to_baseline ~ C(timepoint, Treatment(reference='{reference}')) * "
        "C(compartment, Treatment(reference='core'))"
    )
    result, method, warnings = fit_mixedlm_with_fallback(
        formula,
        subset,
        subset["patient_id"],
        vc_formula={"sample_id": "0 + C(sample_id)"},
    )

    coef_df = coefficients_to_frame(result, metric, "core_vs_accessory", "time_factor_interaction", method)
    interaction_mask = coef_df["term"].str.contains("C\\(timepoint.*:C\\(compartment", regex=True)
    coef_df["p_value_bh_interaction"] = np.nan
    if interaction_mask.any():
        valid = coef_df.loc[interaction_mask, "p_value"].notna()
        if valid.any():
            coef_df.loc[interaction_mask[interaction_mask].index, "p_value_bh_interaction"] = multipletests(
                coef_df.loc[interaction_mask, "p_value"].fillna(1.0), method="fdr_bh"
            )[1]

    means_frames = []
    for compartment in ["core", "accessory"]:
        means_df, _ = compute_estimated_means_and_pairwise(result, subset, reference, compartment=compartment)
        means_df["metric"] = metric
        means_df["compartment"] = compartment
        means_df["method"] = method
        means_frames.append(means_df)
    estimated_means = pd.concat(means_frames, ignore_index=True)

    evidence = "no_clear_evidence"
    interaction_terms = coef_df.loc[interaction_mask].copy()
    if not interaction_terms.empty:
        if interaction_terms["p_value"].lt(0.05).any() or interaction_terms["p_value_bh_interaction"].lt(0.05).any():
            evidence = "interaction_detected"
    interaction_summary = interaction_terms.copy()
    interaction_summary["evidence"] = evidence

    return FitResult(metric, "core_vs_accessory", "time_factor_interaction", method, reference, True, coef_df, estimated_means, interaction_summary, result, warnings)


def confidence_summary(data: pd.DataFrame, group_cols: Sequence[str]) -> pd.DataFrame:
    summary = (
        data.groupby(list(group_cols), dropna=False)["distance_to_baseline"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "estimate", "std": "std_dev", "count": "n"})
    )
    summary["se"] = summary["std_dev"] / np.sqrt(summary["n"].clip(lower=1))
    summary["conf_low"] = summary["estimate"] - 1.96 * summary["se"]
    summary["conf_high"] = summary["estimate"] + 1.96 * summary["se"]
    return summary


def save_pdf_png(fig: plt.Figure, base_path: Path) -> List[Path]:
    pdf_path = base_path.with_suffix(".pdf")
    png_path = base_path.with_suffix(".png")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return [pdf_path, png_path]


def plot_distance_trajectories(distance_long: pd.DataFrame, metric: str, output_prefix: str) -> List[Path]:
    subset = distance_long.loc[distance_long["distance_type"] == metric].copy()
    compartments = ["all", "core", "accessory"]
    summary = confidence_summary(subset, ["compartment", "weeks_post_return"])
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharex=False)

    for ax, compartment in zip(axes, compartments):
        comp_df = subset.loc[subset["compartment"] == compartment].copy()
        comp_summary = summary.loc[summary["compartment"] == compartment].sort_values("weeks_post_return")
        for patient_id, patient_df in comp_df.groupby("patient_id"):
            patient_df = patient_df.sort_values("weeks_post_return")
            ax.plot(patient_df["weeks_post_return"], patient_df["distance_to_baseline"], alpha=0.45, linewidth=0.7)
            ax.scatter(patient_df["weeks_post_return"], patient_df["distance_to_baseline"], s=18, alpha=0.85)
        if not comp_summary.empty:
            ax.fill_between(
                comp_summary["weeks_post_return"].to_numpy(dtype=float),
                comp_summary["conf_low"].to_numpy(dtype=float),
                comp_summary["conf_high"].to_numpy(dtype=float),
                color="black",
                alpha=0.12,
            )
            ax.plot(comp_summary["weeks_post_return"], comp_summary["estimate"], color="black", linewidth=1.1)
            ax.scatter(comp_summary["weeks_post_return"], comp_summary["estimate"], color="black", s=24)
        ax.axhline(0, color="grey", linewidth=0.6)
        ax.set_title({"all": "All taxa", "core": "Core", "accessory": "Accessory"}[compartment], fontweight="bold")
        ax.set_xlabel("Weeks post return")
        ax.set_ylabel(f"{metric.capitalize()} distance to baseline")
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.grid(False)
        sns.despine(ax=ax)

    fig.suptitle(f"Distance trajectories for {metric}: all, core, and accessory", fontsize=14, fontweight="bold")
    fig.tight_layout()
    return save_pdf_png(fig, OUTPUT_DIR / output_prefix)


def plot_time_factor_model(
    distance_long: pd.DataFrame,
    model_results: Dict[Tuple[str, str], FitResult],
    metric: str,
    output_prefix: str,
) -> List[Path]:
    subset = distance_long.loc[(distance_long["distance_type"] == metric) & (distance_long["compartment"].isin(["core", "accessory"]))].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=False)

    for ax, compartment in zip(axes, ["core", "accessory"]):
        comp_df = subset.loc[subset["compartment"] == compartment].copy()
        comp_df["timepoint_str"] = comp_df["timepoint"].astype(str)
        x_map = {tp: idx for idx, tp in enumerate([tp for tp in TIMEPOINT_ORDER if tp in set(comp_df["timepoint_str"])])}
        for patient_id, patient_df in comp_df.groupby("patient_id"):
            patient_df = patient_df.sort_values("weeks_post_return")
            xs = [x_map[tp] for tp in patient_df["timepoint_str"]]
            ax.plot(xs, patient_df["distance_to_baseline"], alpha=0.35, linewidth=0.6)
            jitter = np.random.default_rng(42).normal(0, 0.04, size=len(xs))
            ax.scatter(np.array(xs) + jitter, patient_df["distance_to_baseline"], s=18, alpha=0.75)

        fit_res = model_results[(metric, compartment)]
        means = fit_res.estimated_means.copy()
        means["x"] = means["timepoint"].map(x_map)
        means = means.dropna(subset=["x"]).sort_values("x")
        ax.errorbar(
            means["x"],
            means["estimate"],
            yerr=[means["estimate"] - means["conf_low"], means["conf_high"] - means["estimate"]],
            color="black",
            linewidth=1.1,
            marker="o",
            markersize=5,
            capsize=3,
        )
        ax.axhline(0, color="grey", linewidth=0.6)
        ax.set_xticks(list(x_map.values()))
        ax.set_xticklabels(list(x_map.keys()))
        ax.set_title(f"{compartment.capitalize()} ({metric})", fontweight="bold")
        ax.set_xlabel("Timepoint")
        ax.set_ylabel("Distance to baseline")
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        sns.despine(ax=ax)

    fig.suptitle(f"Time-factor models for core vs accessory: {metric}", fontsize=14, fontweight="bold")
    fig.tight_layout()
    return save_pdf_png(fig, OUTPUT_DIR / output_prefix)


def compute_core_loss_accessory_gain(
    relative_df: pd.DataFrame,
    sample_meta: pd.DataFrame,
    core_features: Sequence[str],
    accessory_features: Sequence[str],
) -> pd.DataFrame:
    included = sample_meta.loc[sample_meta["inclusion_status"] == "included"].copy()
    post_rows = included.loc[included["weeks_post_return"] > 0].copy()
    rows = []

    for patient_id, patient_df in included.groupby("patient_id"):
        patient_df = patient_df.sort_values(["weeks_post_return", "sample_order"])
        baseline_row = patient_df.loc[patient_df["weeks_post_return"] == 0].iloc[0]
        baseline_vec = relative_df[baseline_row["sample_id"]]
        baseline_core_present = [f for f in core_features if baseline_vec.get(f, 0.0) >= DETECTION_THRESHOLD]
        baseline_accessory_absent = [f for f in accessory_features if baseline_vec.get(f, 0.0) < DETECTION_THRESHOLD]

        for row in patient_df.loc[patient_df["weeks_post_return"] > 0].itertuples(index=False):
            post_vec = relative_df[row.sample_id]
            core_loss = [f for f in baseline_core_present if post_vec.get(f, 0.0) < DETECTION_THRESHOLD]
            accessory_gain = [f for f in baseline_accessory_absent if post_vec.get(f, 0.0) >= DETECTION_THRESHOLD]
            rows.append(
                {
                    "sample_id": row.sample_id,
                    "patient_id": patient_id,
                    "baseline_sample_id": baseline_row["sample_id"],
                    "timepoint": row.timepoint,
                    "weeks_post_return": row.weeks_post_return,
                    "core_loss_count": len(core_loss),
                    "core_loss_fraction": len(core_loss) / len(baseline_core_present) if baseline_core_present else np.nan,
                    "accessory_gain_count": len(accessory_gain),
                    "accessory_gain_abundance": float(post_vec.loc[accessory_gain].sum()) if accessory_gain else 0.0,
                }
            )

    return pd.DataFrame(rows)


def plot_core_loss_accessory_gain(loss_gain: pd.DataFrame) -> List[Path]:
    long_df = loss_gain.melt(
        id_vars=["sample_id", "patient_id", "timepoint", "weeks_post_return"],
        value_vars=["core_loss_fraction", "accessory_gain_count"],
        var_name="metric",
        value_name="value",
    )
    long_df["metric"] = long_df["metric"].map(
        {"core_loss_fraction": "Core loss fraction", "accessory_gain_count": "Accessory gain count"}
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharex=False)
    for ax, metric in zip(axes, ["Core loss fraction", "Accessory gain count"]):
        metric_df = long_df.loc[long_df["metric"] == metric].copy()
        order = [tp for tp in TIMEPOINT_ORDER if tp in set(metric_df["timepoint"])]
        sns.violinplot(data=metric_df, x="timepoint", y="value", order=order, color="#d9e6f2", inner=None, ax=ax)
        sns.boxplot(data=metric_df, x="timepoint", y="value", order=order, width=0.18, color="white", fliersize=0, ax=ax)
        sns.stripplot(data=metric_df, x="timepoint", y="value", order=order, hue="patient_id", dodge=False, size=4, alpha=0.8, ax=ax)
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        ax.set_title(metric, fontweight="bold")
        ax.set_xlabel("Timepoint")
        ax.set_ylabel("")
        sns.despine(ax=ax)
    fig.suptitle("Core loss and accessory gain over time", fontsize=14, fontweight="bold")
    fig.tight_layout()
    return save_pdf_png(fig, OUTPUT_DIR / "fig_core_loss_accessory_gain_over_time")


def trapz_auc(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 2:
        return np.nan
    order = np.argsort(x)
    return float(np.trapezoid(y[order], x[order]))


def compute_recovery_time(weeks: np.ndarray, values: np.ndarray, threshold: float = RECOVERY_THRESHOLD_BRAY) -> float:
    mask = np.isfinite(weeks) & np.isfinite(values)
    weeks = weeks[mask]
    values = values[mask]
    if weeks.size == 0:
        return np.nan
    peak_idx = int(np.nanargmax(values))
    for idx in range(peak_idx, len(values)):
        if values[idx] <= threshold:
            return float(weeks[idx])
    return np.nan


def compute_trajectory_metrics(distance_long: pd.DataFrame, loss_gain: pd.DataFrame) -> pd.DataFrame:
    bray = distance_long.loc[distance_long["distance_type"] == "bray"].copy()
    rows = []
    for patient_id, patient_df in bray.groupby("patient_id"):
        patient_df = patient_df.sort_values("weeks_post_return")
        loss_df = loss_gain.loc[loss_gain["patient_id"] == patient_id].sort_values("weeks_post_return")
        row = {"patient_id": patient_id}
        for compartment in ["core", "accessory"]:
            comp = patient_df.loc[patient_df["compartment"] == compartment]
            values = comp["distance_to_baseline"].to_numpy(dtype=float)
            weeks = comp["weeks_post_return"].to_numpy(dtype=float)
            row[f"peak_{compartment}_bray"] = float(np.nanmax(values)) if values.size else np.nan
            row[f"auc_{compartment}_bray"] = trapz_auc(weeks, values)
            row[f"recovery_time_{compartment}_bray"] = compute_recovery_time(weeks, values)
            row[f"volatility_{compartment}_bray"] = float(np.nanmean(np.abs(np.diff(values)))) if values.size >= 2 else np.nan
        row["core_loss_fraction_peak"] = float(loss_df["core_loss_fraction"].max()) if not loss_df.empty else np.nan
        row["accessory_gain_count_peak"] = float(loss_df["accessory_gain_count"].max()) if not loss_df.empty else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def cluster_trajectories(metrics_df: pd.DataFrame, n_clusters: int = DEFAULT_N_CLUSTERS) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    metric_cols = [col for col in metrics_df.columns if col != "patient_id"]
    valid_cols = [col for col in metric_cols if metrics_df[col].notna().sum() >= 2 and metrics_df[col].std(skipna=True) > 0]
    if len(valid_cols) == 0 or metrics_df.shape[0] < 2:
        cluster_df = metrics_df[["patient_id"]].copy()
        cluster_df["cluster_id"] = 1
        cluster_df["cluster_label"] = "cluster_1"
        return metrics_df.copy(), cluster_df, "not_performed_insufficient_variation"

    matrix = metrics_df[valid_cols].copy()
    matrix = matrix.fillna(matrix.median())
    scaled = StandardScaler().fit_transform(matrix)
    linkage_mat = linkage(scaled, method="ward")
    k_used = min(max(1, n_clusters), metrics_df.shape[0])
    cluster_ids = fcluster(linkage_mat, t=k_used, criterion="maxclust")

    cluster_df = metrics_df[["patient_id"]].copy()
    cluster_df["cluster_id"] = cluster_ids
    cluster_df["cluster_label"] = cluster_df["cluster_id"].map(lambda x: f"cluster_{x}")

    z_metrics = pd.DataFrame(scaled, columns=valid_cols)
    z_metrics.insert(0, "patient_id", metrics_df["patient_id"].values)
    z_metrics = z_metrics.merge(cluster_df, on="patient_id", how="left")
    return z_metrics, cluster_df, "hierarchical_clustering_ward"


def plot_cluster_heatmap(z_metrics: pd.DataFrame) -> List[Path]:
    if z_metrics.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "Insufficient variation for clustering", ha="center", va="center")
        ax.axis("off")
        return save_pdf_png(fig, OUTPUT_DIR / "fig_core_accessory_cluster_heatmap")

    metric_cols = [col for col in z_metrics.columns if col not in {"patient_id", "cluster_id", "cluster_label"}]
    plot_df = z_metrics.sort_values(["cluster_id", "patient_id"]).set_index("patient_id")[metric_cols]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.heatmap(plot_df, cmap="vlag", center=0, linewidths=0.5, linecolor="white", ax=ax)
    ax.set_title("Cluster heatmap of core/accessory trajectory metrics", fontweight="bold")
    ax.set_xlabel("Metric")
    ax.set_ylabel("Patient")
    fig.tight_layout()
    return save_pdf_png(fig, OUTPUT_DIR / "fig_core_accessory_cluster_heatmap")


def write_summary(
    input_path: Path,
    orientation: str,
    sample_meta: pd.DataFrame,
    feature_table: pd.DataFrame,
    separate_results: Dict[Tuple[str, str], FitResult],
    joint_results: Dict[str, FitResult],
    output_paths: List[Path],
    reused_components: List[str],
    warnings: List[str],
) -> str:
    included = sample_meta.loc[sample_meta["inclusion_status"] == "included"].copy()
    timepoints = sorted(included.loc[included["weeks_post_return"] > 0, "timepoint"].astype(str).unique(), key=lambda x: WEEK_MAP[x])
    core_n = int(feature_table["is_core"].sum())
    accessory_n = int((~feature_table["is_core"]).sum())

    joint_evidence_lines = []
    for metric, fit in joint_results.items():
        interaction_terms = fit.coefficients.loc[fit.coefficients["term"].str.contains(":")]
        evidence = "no_clear_evidence"
        if not interaction_terms.empty and ((interaction_terms["p_value"] < 0.05).any() or (interaction_terms["p_value_bh_interaction"] < 0.05).any()):
            evidence = "evidence_for_different_core_accessory_trajectories"
        joint_evidence_lines.append(f"- {metric}: {evidence}")

    lines = [
        "Core/accessory split of MetaPhlAn SGB data in Python",
        f"Input file: {input_path.resolve()}",
        f"Detected orientation: {orientation}",
        f"Samples total: {sample_meta.shape[0]}",
        f"Included persons: {included['patient_id'].nunique()}",
        f"Included post samples: {(included['weeks_post_return'] > 0).sum()}",
        f"Core features: {core_n}",
        f"Accessory features: {accessory_n}",
        "Core definition: cohort-core from w0 only, prevalence >= 0.70, detection >= 1e-4",
        "Distance subsets treated as independent data spaces for Bray and CLR/Aitchison within each compartment.",
        f"Reused exploration components: {', '.join(reused_components)}",
        f"Separate core/accessory models fit: {all(f.fitted for f in separate_results.values())}",
        f"Joint interaction models fit: {all(f.fitted for f in joint_results.values())}",
        f"Timepoints included: {', '.join(timepoints)}",
        "Evidence for different core/accessory trajectories:",
        *joint_evidence_lines,
        "Warnings / limitations:",
    ]
    lines.extend([f"- {warning}" for warning in warnings] if warnings else ["- none"])
    lines.append("Generated files:")
    lines.extend([f"- {path.resolve()}" for path in output_paths])
    summary = "\n".join(lines)
    summary_path = OUTPUT_DIR / "core_accessory_time_factor_summary.txt"
    summary_path.write_text(summary)
    return summary


def main() -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    np.random.seed(42)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    warnings: List[str] = []
    reused_components = [
        "sample-name parsing and baseline logic from the existing distance-to-baseline workflow",
        "same Bray/Aitchison baseline-distance concept as in first exploration",
        "same trajectory/time-factor plot logic adapted to all/core/accessory",
    ]

    input_path = resolve_input_path()
    if input_path.name == "TAXA_matrix.tsv":
        warnings.append("Preferred MetaPhlAn table not found; using fallback data/processed/taxa/TAXA_matrix.tsv.")

    feature_df, orientation = load_metaphlan_table(input_path)
    sample_meta, patient_summary = build_sample_metadata(feature_df.columns.tolist())
    included_samples = sample_meta.loc[sample_meta["inclusion_status"] == "included", "sample_id"].tolist()
    feature_df = feature_df.loc[:, included_samples]
    relative_df = make_relative_abundance(feature_df)
    relative_df, prevalence_table = apply_global_prevalence_filter(relative_df)

    unsupported = sample_meta.loc[sample_meta["parse_status"] == "unsupported_timepoint", "sample_id"].tolist()
    if unsupported:
        warnings.append(f"Unsupported timepoints excluded: {', '.join(unsupported)}")
    if sample_meta["parse_status"].eq("unparsed").any():
        warnings.append(
            "Unparsed sample names excluded: " + ", ".join(sample_meta.loc[sample_meta["parse_status"] == "unparsed", "sample_id"].tolist())
        )
    multiple_w0 = patient_summary.loc[patient_summary["multiple_w0"], "patient_id"].tolist()
    if multiple_w0:
        warnings.append(f"Multiple w0 samples detected for {', '.join(multiple_w0)}; first baseline sample was used.")

    feature_table, core_features, accessory_features = define_core_accessory(relative_df, sample_meta)
    if len(accessory_features) > 5 * max(1, len(core_features)):
        warnings.append("Accessory compartment is much larger than the core compartment and may be relatively sparse.")

    feature_table.to_csv(OUTPUT_DIR / "core_accessory_feature_table.csv", index=False)
    prevalence_table.to_csv(OUTPUT_DIR / "core_accessory_global_prevalence.csv", index=False)

    distance_long, distance_wide = make_distance_outputs(relative_df, sample_meta, core_features, accessory_features)
    distance_long.to_csv(OUTPUT_DIR / "core_accessory_distance_long.csv", index=False)
    distance_wide.to_csv(OUTPUT_DIR / "core_accessory_distance_wide.csv", index=False)

    separate_results: Dict[Tuple[str, str], FitResult] = {}
    output_paths: List[Path] = [
        OUTPUT_DIR / "core_accessory_feature_table.csv",
        OUTPUT_DIR / "core_accessory_global_prevalence.csv",
        OUTPUT_DIR / "core_accessory_distance_long.csv",
        OUTPUT_DIR / "core_accessory_distance_wide.csv",
    ]

    for metric in ["bray", "aitchison"]:
        for compartment in ["core", "accessory"]:
            fit = fit_time_factor_model(distance_long, compartment, metric)
            separate_results[(metric, compartment)] = fit
            base_name = f"model_{compartment}_{metric}_time_factor"
            fit.coefficients.to_csv(OUTPUT_DIR / f"{base_name}.csv", index=False)
            fit.estimated_means.to_csv(OUTPUT_DIR / f"{base_name}_emmeans.csv", index=False)
            fit.pairwise.to_csv(OUTPUT_DIR / f"{base_name}_pairwise.csv", index=False)
            output_paths.extend(
                [
                    OUTPUT_DIR / f"{base_name}.csv",
                    OUTPUT_DIR / f"{base_name}_emmeans.csv",
                    OUTPUT_DIR / f"{base_name}_pairwise.csv",
                ]
            )
            warnings.extend(fit.warnings)

    joint_results: Dict[str, FitResult] = {}
    for metric in ["bray", "aitchison"]:
        fit = fit_joint_interaction_model(distance_long, metric)
        joint_results[metric] = fit
        fit.coefficients.to_csv(OUTPUT_DIR / f"model_joint_{metric}_time_factor_interaction.csv", index=False)
        fit.estimated_means.to_csv(OUTPUT_DIR / f"model_joint_{metric}_time_factor_interaction_emmeans.csv", index=False)
        output_paths.extend(
            [
                OUTPUT_DIR / f"model_joint_{metric}_time_factor_interaction.csv",
                OUTPUT_DIR / f"model_joint_{metric}_time_factor_interaction_emmeans.csv",
            ]
        )
        warnings.extend(fit.warnings)

    output_paths.extend(plot_distance_trajectories(distance_long, "bray", "fig_distance_trajectories_bray_all_core_accessory"))
    output_paths.extend(plot_distance_trajectories(distance_long, "aitchison", "fig_distance_trajectories_aitchison_all_core_accessory"))
    output_paths.extend(plot_time_factor_model(distance_long, separate_results, "bray", "fig_time_factor_core_accessory_bray"))
    output_paths.extend(plot_time_factor_model(distance_long, separate_results, "aitchison", "fig_time_factor_core_accessory_aitchison"))

    loss_gain = compute_core_loss_accessory_gain(relative_df, sample_meta, core_features, accessory_features)
    loss_gain.to_csv(OUTPUT_DIR / "core_loss_accessory_gain.csv", index=False)
    output_paths.append(OUTPUT_DIR / "core_loss_accessory_gain.csv")
    output_paths.extend(plot_core_loss_accessory_gain(loss_gain))

    trajectory_metrics = compute_trajectory_metrics(distance_long, loss_gain)
    trajectory_metrics.to_csv(OUTPUT_DIR / "core_accessory_trajectory_metrics.csv", index=False)
    output_paths.append(OUTPUT_DIR / "core_accessory_trajectory_metrics.csv")

    z_metrics, cluster_df, cluster_method = cluster_trajectories(trajectory_metrics)
    cluster_df.to_csv(OUTPUT_DIR / "core_accessory_trajectory_clusters.csv", index=False)
    output_paths.append(OUTPUT_DIR / "core_accessory_trajectory_clusters.csv")
    output_paths.extend(plot_cluster_heatmap(z_metrics))

    warnings = list(dict.fromkeys(warnings))
    summary_text = write_summary(
        input_path=input_path,
        orientation=orientation,
        sample_meta=sample_meta,
        feature_table=feature_table,
        separate_results=separate_results,
        joint_results=joint_results,
        output_paths=output_paths,
        reused_components=reused_components,
        warnings=warnings,
    )
    output_paths.append(OUTPUT_DIR / "core_accessory_time_factor_summary.txt")

    interaction_text = []
    for metric, fit in joint_results.items():
        interaction_terms = fit.coefficients.loc[fit.coefficients["term"].str.contains(":")]
        evidence = "no_clear_evidence"
        if not interaction_terms.empty and ((interaction_terms["p_value"] < 0.05).any() or (interaction_terms["p_value_bh_interaction"] < 0.05).any()):
            evidence = "evidence_for_different_trajectories"
        interaction_text.append(f"{metric}: {evidence}")

    print("\nCore/accessory time-factor analysis in Python completed.")
    print("Reused exploration:", "; ".join(reused_components))
    print(f"Core features: {len(core_features)}")
    print(f"Accessory features: {len(accessory_features)}")
    print("Models successfully fitted:")
    for (metric, compartment), fit in separate_results.items():
        print(f" - separate {compartment} {metric}: {fit.method}")
    for metric, fit in joint_results.items():
        print(f" - joint interaction {metric}: {fit.method}")
    print("Joint interaction evidence:")
    for line in interaction_text:
        print(f" - {line}")
    print("Generated files:")
    for path in output_paths:
        print(f" - {path.resolve()}")
    if warnings:
        print("Warnings / limitations:")
        for warning in warnings:
            print(f" - {warning}")
    print("\nSummary:")
    print(summary_text)


if __name__ == "__main__":
    main()
