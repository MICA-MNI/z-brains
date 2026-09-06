#!/usr/bin/env python3
"""Plot pooled AUC benchmark values by processing-feature category."""

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import pandas as pd


BENCHMARK_DIR = Path(__file__).resolve().parent / "processing_method_benchmark"
INPUT_CSV = BENCHMARK_DIR / "choice_consensus.csv"
OUTPUT_PNG = BENCHMARK_DIR / "processing_feature_auc_by_category.png"
OUTPUT_SVG = BENCHMARK_DIR / "processing_feature_auc_by_category.svg"

TOP_N_PER_MODALITY = 4

FAMILY_ORDER = [
    "depth_contrast",
    "intensity_depth_self_normalization",
    "gw_transition_blur",
    "blur_quantification",
    "gradient_flattening_self_normalization",
    "multisurface_intensity",
    "intensity_self_normalization",
]

FAMILY_LABELS = {
    "depth_contrast": "Depth Contrast",
    "intensity_depth_self_normalization": "Intensity-Depth Self-Norm",
    "gw_transition_blur": "Gray-White Transition",
    "blur_quantification": "Blur / Profile Distance",
    "gradient_flattening_self_normalization": "Gradient Flattening Self-Norm",
    "multisurface_intensity": "Multisurface Intensity",
    "intensity_self_normalization": "Intensity Self-Norm",
}

MODALITY_COLORS = {
    "FLAIR": "#2f7f95",
    "T1w": "#b85c38",
}


def wrapped_label(modality: str, strategy: str, width: int = 34) -> str:
    label = f"{modality}: {strategy}"
    return "\n".join(textwrap.wrap(label, width=width, break_long_words=False))


def top_rows_by_family(df: pd.DataFrame, family: str) -> pd.DataFrame:
    rows = []
    family_df = df[df["family"].eq(family)]
    for modality in ["FLAIR", "T1w"]:
        modality_df = family_df[family_df["modality"].eq(modality)]
        top = modality_df.sort_values(
            ["mean_primary_score", "top3_count", "top1_count"],
            ascending=[False, False, False],
        ).head(TOP_N_PER_MODALITY)
        rows.append(top)
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    df = pd.read_csv(INPUT_CSV)
    auc_df = df[
        df["axis"].eq("feature_processing")
        & df["metric_name"].eq("pooled_auc")
        & df["family"].isin(FAMILY_ORDER)
    ].copy()

    fig, axes = plt.subplots(4, 2, figsize=(18, 22), constrained_layout=True)
    axes = axes.ravel()

    for ax, family in zip(axes, FAMILY_ORDER):
        panel = top_rows_by_family(auc_df, family)
        panel = panel.sort_values(
            ["modality", "mean_primary_score"],
            ascending=[True, True],
        ).reset_index(drop=True)

        colors = [MODALITY_COLORS.get(modality, "#777777") for modality in panel["modality"]]
        xerr = [
            panel["mean_primary_score"] - panel["min_primary_score"],
            panel["max_primary_score"] - panel["mean_primary_score"],
        ]
        labels = [
            wrapped_label(row.modality, row.strategy)
            for row in panel.itertuples(index=False)
        ]

        bars = ax.barh(
            range(len(panel)),
            panel["mean_primary_score"],
            xerr=xerr,
            color=colors,
            alpha=0.88,
            error_kw={"elinewidth": 1, "capthick": 1, "capsize": 3, "ecolor": "#333333"},
        )

        ax.set_yticks(range(len(panel)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(0.45, 0.9)
        ax.set_xlabel("Mean pooled AUC across benchmark splits")
        ax.set_title(FAMILY_LABELS.get(family, family), fontsize=13, fontweight="bold")
        ax.grid(axis="x", color="#d0d0d0", linewidth=0.6, alpha=0.7)
        ax.set_axisbelow(True)

        for bar, score in zip(bars, panel["mean_primary_score"]):
            ax.text(
                min(score + 0.007, 0.885),
                bar.get_y() + bar.get_height() / 2,
                f"{score:.3f}",
                va="center",
                ha="left",
                fontsize=8,
                color="#111111",
            )

    extra_ax = axes[len(FAMILY_ORDER)]
    extra_ax.axis("off")
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=color, alpha=0.88)
        for color in MODALITY_COLORS.values()
    ]
    extra_ax.legend(
        legend_handles,
        list(MODALITY_COLORS.keys()),
        loc="center",
        title="Modality",
        frameon=False,
        fontsize=11,
        title_fontsize=12,
    )
    extra_ax.text(
        0.5,
        0.30,
        "Bars show mean pooled AUC across MICs Overall,\n"
        "NOEL Overall, and pooled benchmark splits.\n"
        "Error bars show min-max AUC across those splits.",
        ha="center",
        va="center",
        fontsize=10,
        color="#333333",
    )

    fig.suptitle(
        "Processing Feature Benchmark: Top Pooled AUC Values by Feature Category",
        fontsize=18,
        fontweight="bold",
    )
    fig.savefig(OUTPUT_PNG, dpi=300)
    fig.savefig(OUTPUT_SVG)
    print(f"Wrote {OUTPUT_PNG}")
    print(f"Wrote {OUTPUT_SVG}")


if __name__ == "__main__":
    main()
