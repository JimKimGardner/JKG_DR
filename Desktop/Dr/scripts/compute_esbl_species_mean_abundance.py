#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ELIGIBLE_SAMPLES_PATH = PROJECT_ROOT / "results" / "esbl_non_esbl_daa" / "esbl_eligible_samples.tsv"
RELATIVE_ABUNDANCE_PATH = PROJECT_ROOT / "data" / "processed" / "taxa" / "derived" / "TAXA_relative_abundance.tsv"
FEATURE_METADATA_PATH = PROJECT_ROOT / "data" / "processed" / "taxa" / "TAXA_feature_metadata.tsv"
OUTPUT_PATH = PROJECT_ROOT / "results" / "esbl_non_esbl_daa" / "esbl_non_esbl_species_mean_abundance.tsv"


def main() -> None:
    eligible_samples = pd.read_csv(ELIGIBLE_SAMPLES_PATH, sep="\t")
    relative_abundance = pd.read_csv(RELATIVE_ABUNDANCE_PATH, sep="\t")
    feature_metadata = pd.read_csv(FEATURE_METADATA_PATH, sep="\t")

    species_meta = (
        feature_metadata.loc[feature_metadata["terminal_rank"] == "SGB", [
            "feature_id",
            "species_label",
            "sgb_label",
            "display_label",
            "prevalence",
            "n_present_samples",
        ]]
        .copy()
    )

    selected_sample_ids = eligible_samples["sample_id"].tolist()
    group_map = eligible_samples.set_index("sample_id")["esbl_status"].to_dict()

    species_df = (
        relative_abundance.loc[
            relative_abundance["feature_id"].isin(species_meta["feature_id"]),
            ["feature_id", *selected_sample_ids],
        ]
        .copy()
        .set_index("feature_id")
    )

    transposed = species_df.transpose()
    transposed.index.name = "sample_id"
    transposed = transposed.reset_index()
    transposed["esbl_status"] = transposed["sample_id"].map(group_map)

    long_df = transposed.melt(
        id_vars=["sample_id", "esbl_status"],
        var_name="feature_id",
        value_name="relative_abundance",
    )

    summary = (
        long_df.groupby(["feature_id", "esbl_status"], as_index=False)
        .agg(
            n_samples=("sample_id", "size"),
            n_present=("relative_abundance", lambda x: int((x > 0).sum())),
            mean_relative_abundance=("relative_abundance", "mean"),
            median_relative_abundance=("relative_abundance", "median"),
        )
    )

    test_rows = []
    for feature_id, feature_df in long_df.groupby("feature_id", sort=False):
        esbl_values = feature_df.loc[
            feature_df["esbl_status"] == "esbl_positive", "relative_abundance"
        ].to_numpy()
        no_esbl_values = feature_df.loc[
            feature_df["esbl_status"] == "no_esbl", "relative_abundance"
        ].to_numpy()

        welch_t_pvalue = float("nan")
        mann_whitney_pvalue = float("nan")

        if len(esbl_values) >= 2 and len(no_esbl_values) >= 2:
            welch_t_pvalue = stats.ttest_ind(
                esbl_values,
                no_esbl_values,
                equal_var=False,
                nan_policy="omit",
            ).pvalue

        if len(esbl_values) >= 1 and len(no_esbl_values) >= 1:
            mann_whitney_pvalue = stats.mannwhitneyu(
                esbl_values,
                no_esbl_values,
                alternative="two-sided",
            ).pvalue

        test_rows.append(
            {
                "feature_id": feature_id,
                "welch_t_pvalue": welch_t_pvalue,
                "mann_whitney_pvalue": mann_whitney_pvalue,
            }
        )

    test_df = pd.DataFrame(test_rows)

    valid_welch = test_df["welch_t_pvalue"].notna()
    test_df["welch_t_qvalue_bh"] = float("nan")
    if valid_welch.any():
        test_df.loc[valid_welch, "welch_t_qvalue_bh"] = multipletests(
            test_df.loc[valid_welch, "welch_t_pvalue"],
            method="fdr_bh",
        )[1]

    valid_mw = test_df["mann_whitney_pvalue"].notna()
    test_df["mann_whitney_qvalue_bh"] = float("nan")
    if valid_mw.any():
        test_df.loc[valid_mw, "mann_whitney_qvalue_bh"] = multipletests(
            test_df.loc[valid_mw, "mann_whitney_pvalue"],
            method="fdr_bh",
        )[1]

    wide = (
        summary.pivot(index="feature_id", columns="esbl_status")
        .sort_index(axis=1)
    )
    wide.columns = [f"{metric}_{group}" for metric, group in wide.columns]
    wide = wide.reset_index()

    result = (
        species_meta.merge(wide, on="feature_id", how="left")
        .merge(test_df, on="feature_id", how="left")
    )

    result["feature_short"] = result["species_label"].fillna("").str.strip()
    missing_species = result["feature_short"] == ""
    result.loc[missing_species, "feature_short"] = result.loc[missing_species, "display_label"].fillna("")
    still_missing = result["feature_short"] == ""
    result.loc[still_missing, "feature_short"] = result.loc[still_missing, "feature_id"]
    result["feature_short"] = result["feature_short"] + " / " + result["sgb_label"].fillna(result["feature_id"])

    result["difference_esbl_positive_minus_no_esbl"] = (
        result["mean_relative_abundance_esbl_positive"] - result["mean_relative_abundance_no_esbl"]
    )
    result["difference_no_esbl_minus_esbl_positive"] = (
        result["mean_relative_abundance_no_esbl"] - result["mean_relative_abundance_esbl_positive"]
    )
    result["absolute_difference"] = (
        result["difference_esbl_positive_minus_no_esbl"].abs()
    )

    result = result.sort_values(
        by=["absolute_difference", "difference_esbl_positive_minus_no_esbl"],
        ascending=[False, False],
        na_position="last",
    )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT_PATH, sep="\t", index=False)

    print(f"Wrote {OUTPUT_PATH}")
    print(f"Rows: {len(result)}")


if __name__ == "__main__":
    main()
