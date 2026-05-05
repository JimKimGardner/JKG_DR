from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
from matplotlib.colors import ListedColormap


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "results" / "esbl_status"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STATUS_TSV = OUT_DIR / "esbl_status_from_screenshot.tsv"
MATCHED_TSV = OUT_DIR / "esbl_status_metaphlan_matched.tsv"
HEATMAP_PNG = OUT_DIR / "esbl_status_from_screenshot_heatmap.png"

# Cell geometry was calibrated on the provided screenshot.
X_BOUNDS = [108, 177, 246, 315, 384, 453, 522, 591, 660, 729, 798, 871]
X_CENTERS = [(X_BOUNDS[i] + X_BOUNDS[i + 1]) // 2 for i in range(len(X_BOUNDS) - 1)]

ROW_LABELS = [
    "P40", "P39", "P38", "P37", "P36", "P35", "P34", "P32", "P31", "P30",
    "P29", "P28", "P27", "P26", "P25", "P24", "P22", "P21", "P20", "P19",
    "P18", "P17", "P16", "P15", "P14", "P13", "P12", "P11", "P10", "P9",
    "P8", "P7", "P6", "P5", "P4", "P3", "P2", "P1",
]
Y_START = 51
Y_STEP = 23.45
Y_CENTERS = [round(Y_START + i * Y_STEP) for i in range(len(ROW_LABELS))]

SCREENSHOT_TIMEPOINTS = ["w-1", "w0", "w2", "w4", "w6", "w8", "w10", "w12", "w16", "w20", "w52"]
SEQUENCING_TIMEPOINT_MAP = {
    "w-1": "w0",
    "w0": "w1",
    "w2": "w2",
    "w4": "w4",
    "w6": "w6",
    "w8": "w8",
    "w10": "w10",
    "w12": "w12",
    "w16": "w16",
    "w20": "w20",
    "w52": "w52",
}
SEQUENCING_WEEK_MAP = {
    "w-1": 0,
    "w0": 1,
    "w2": 2,
    "w4": 4,
    "w6": 6,
    "w8": 8,
    "w10": 10,
    "w12": 12,
    "w16": 16,
    "w20": 20,
    "w52": 52,
}

RGB_TO_STATUS = {
    (240, 201, 176): "esbl_positive",
    (200, 229, 243): "no_esbl",
    (166, 166, 166): "na",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract ESBL status calls from the calibrated screenshot."
    )
    parser.add_argument(
        "--image",
        dest="image_path",
        default=os.environ.get(
            "ESBL_SCREENSHOT_PATH",
            "/Users/jimmyg/Library/Mobile Documents/com~apple~CloudDocs/Images/Screenshot 2026-05-04 at 11.17.13.png",
        ),
        help="Path to the screenshot image used for ESBL status extraction.",
    )
    return parser.parse_args()


def classify_rgb(rgb: tuple[int, int, int]) -> str:
    best_label = None
    best_distance = None
    for ref_rgb, label in RGB_TO_STATUS.items():
        distance = sum((channel - ref_channel) ** 2 for channel, ref_channel in zip(rgb, ref_rgb))
        if best_distance is None or distance < best_distance:
            best_distance = distance
            best_label = label
    return best_label or "unknown"


def participant_to_id(label: str) -> str:
    return f"p{int(label[1:]):03d}"


def build_status_table(image_path: Path) -> pd.DataFrame:
    image = Image.open(image_path).convert("RGB")
    records: list[dict[str, str | int | tuple[int, int, int]]] = []

    for row_label, y in zip(ROW_LABELS, Y_CENTERS):
        for screenshot_timepoint, x in zip(SCREENSHOT_TIMEPOINTS, X_CENTERS):
            rgb = image.getpixel((x, y))
            records.append(
                {
                    "participant_label": row_label,
                    "participant_id": participant_to_id(row_label),
                    "screenshot_timepoint": screenshot_timepoint,
                    "sequencing_timepoint": SEQUENCING_TIMEPOINT_MAP[screenshot_timepoint],
                    "sequencing_week": SEQUENCING_WEEK_MAP[screenshot_timepoint],
                    "cell_rgb": ",".join(str(channel) for channel in rgb),
                    "esbl_status": classify_rgb(rgb),
                }
            )

    return pd.DataFrame.from_records(records)


def build_metaphlan_matched_table(status_df: pd.DataFrame) -> pd.DataFrame:
    sample_meta_path = PROJECT_ROOT / "data" / "processed" / "taxa" / "TAXA_sample_metadata.tsv"
    sample_meta = pd.read_csv(sample_meta_path, sep="\t", dtype=str)
    sample_meta["participant_id"] = sample_meta["participant_id"].astype(str)
    sample_meta["timepoint_code"] = sample_meta["timepoint_code"].astype(str).str.lower()
    sample_meta["week"] = sample_meta["week"].astype(int)

    joined = status_df.merge(
        sample_meta[["sample_id", "participant_id", "timepoint_code", "week"]],
        left_on=["participant_id", "sequencing_week"],
        right_on=["participant_id", "week"],
        how="inner",
    )

    joined = joined[
        [
            "sample_id",
            "participant_label",
            "participant_id",
            "screenshot_timepoint",
            "sequencing_timepoint",
            "sequencing_week",
            "week",
            "esbl_status",
        ]
    ].sort_values(["participant_id", "week", "sample_id"])

    return joined


def save_heatmap(status_df: pd.DataFrame) -> None:
    plot_df = status_df.pivot(index="participant_label", columns="screenshot_timepoint", values="esbl_status")
    plot_df = plot_df.loc[ROW_LABELS, SCREENSHOT_TIMEPOINTS]

    value_map = {"esbl_positive": 0, "no_esbl": 1, "na": 2}
    heatmap_values = plot_df.apply(lambda column: column.map(value_map)).to_numpy()

    cmap = ListedColormap(["#f0c9b0", "#c8e5f3", "#a6a6a6"])

    fig, ax = plt.subplots(figsize=(10.8, 12.2))
    ax.imshow(heatmap_values, cmap=cmap, aspect="auto", interpolation="nearest")

    ax.set_xticks(range(len(SCREENSHOT_TIMEPOINTS)))
    ax.set_xticklabels(SCREENSHOT_TIMEPOINTS, fontsize=10)
    ax.set_yticks(range(len(ROW_LABELS)))
    ax.set_yticklabels(ROW_LABELS, fontsize=8)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, length=0)

    ax.set_xticks([x - 0.5 for x in range(1, len(SCREENSHOT_TIMEPOINTS))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(ROW_LABELS))], minor=True)
    ax.grid(which="minor", color="#7b8794", linestyle="-", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.set_title("ESBL E. coli status extracted from screenshot", fontsize=15, loc="left", pad=24)

    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor="#f0c9b0", edgecolor="#111827", label="ESBL E. coli detected"),
        Patch(facecolor="#c8e5f3", edgecolor="#111827", label="No ESBL E. coli"),
        Patch(facecolor="#a6a6a6", edgecolor="#111827", label="NA"),
    ]
    ax.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.07), ncol=3, frameon=False)

    fig.tight_layout()
    fig.savefig(HEATMAP_PNG, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    image_path = Path(args.image_path).expanduser().resolve()
    if not image_path.exists():
        raise FileNotFoundError(f"Screenshot image not found: {image_path}")

    status_df = build_status_table(image_path)
    status_df.to_csv(STATUS_TSV, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    matched_df = build_metaphlan_matched_table(status_df)
    matched_df.to_csv(MATCHED_TSV, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    save_heatmap(status_df)

    print(f"Created {STATUS_TSV}")
    print(f"Created {MATCHED_TSV}")
    print(f"Created {HEATMAP_PNG}")
    print(f"Image used: {image_path}")


if __name__ == "__main__":
    main()
