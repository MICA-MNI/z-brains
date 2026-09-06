#!/usr/bin/env python3
"""Plot disease-specific AUC by processing category rather than benchmark family."""

from __future__ import annotations

from math import ceil
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_DIR = REPO_ROOT / "results" / "processing_method_benchmark"
SCREENING_CSV = BENCHMARK_DIR / "screening_long.csv"

TOP_N_PER_MODALITY = 4
TOP_N_WSCORE_PER_FEATURE = 4
XLIM = (0.35, 0.90)

DATASET_COLORS = {
    "MICs": "#3f77b5",
    "NOEL": "#c26b3f",
}

CATEGORIES = [
    {
        "key": "intensity_feature_definition",
        "title": "Intensity Feature Processing",
        "description": "T1w/FLAIR features from depth contrasts, depth normalization, and multisurface summaries.",
        "families": {"depth_contrast", "intensity_depth_self_normalization", "multisurface_intensity"},
    },
    {
        "key": "blur_feature_definition",
        "title": "Blur Features From Intensities",
        "description": "Blur/profile and gray-white transition features calculated from intensity profiles.",
        "families": {"blur_quantification", "gw_transition_blur"},
    },
    {
        "key": "blur_preprocessing",
        "title": "Blur Feature Pre-Processing",
        "description": "Pre-processing/normalization of intensity-derived blur maps.",
        "families": {"gradient_flattening_self_normalization"},
    },
    {
        "key": "intensity_preprocessing",
        "title": "Intensity Pre-Processing",
        "description": "Subject-level pre-W-score centering, scaling, and rank/robust transforms.",
        "families": {"intensity_self_normalization"},
        "strategy_prefixes": ("pre_",),
    },
    {
        "key": "intensity_postprocessing",
        "title": "Intensity Post-Processing",
        "description": "Post-W-score spatial centering, scaling, and robust transforms.",
        "families": {"intensity_self_normalization"},
        "strategy_prefixes": ("post_",),
    },
    {
        "key": "covariate_partial_normalization",
        "title": "Covariate / Partial Normalization",
        "description": "Covariate-adjusted and partial subject-level normalization variants.",
        "families": {"intensity_self_normalization"},
        "strategy_prefixes": ("covariate_", "partial_", "raw_wscore"),
    },
]

FAMILY_SHORT = {
    "depth_contrast": "depth",
    "intensity_depth_self_normalization": "depth-norm",
    "multisurface_intensity": "multisurface",
    "blur_quantification": "profile",
    "gw_transition_blur": "transition",
    "gradient_flattening_self_normalization": "gradient-flat",
    "intensity_self_normalization": "self-norm",
    "tle_adc_fa_wscore": "w-score",
}


def wrap_label(text: str, width: int = 40) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def category_rows(screening: pd.DataFrame, disease: str, category: dict[str, object]) -> pd.DataFrame:
    rows = screening[
        screening["axis"].eq("feature_processing")
        & screening["metric_name"].eq("mean_subject_auc_by_disease")
        & screening["disease"].eq(disease)
        & screening["dataset"].isin(["MICs", "NOEL"])
        & screening["family"].isin(category["families"])
    ].copy()
    prefixes = category.get("strategy_prefixes")
    if prefixes:
        rows = rows[
            rows["strategy"].astype(str).apply(
                lambda strategy: any(strategy.startswith(prefix) for prefix in prefixes)
            )
        ].copy()
    elif category["families"] == {"intensity_self_normalization"}:
        rows = rows.iloc[0:0].copy()

    if rows.empty:
        return rows

    rows["source_family"] = rows["family"].map(FAMILY_SHORT).fillna(rows["family"])
    selected = []
    for modality in ["T1w", "FLAIR"]:
        modality_rows = rows[rows["modality"].eq(modality)]
        if modality_rows.empty:
            continue
        summary = (
            modality_rows.groupby(["family", "source_family", "modality", "strategy"], dropna=False)
            .agg(mean_auc=("primary_score", "mean"), available_datasets=("dataset", "nunique"))
            .reset_index()
            .sort_values(["mean_auc", "available_datasets"], ascending=[False, False])
            .head(TOP_N_PER_MODALITY)
        )
        selected.append(
            modality_rows.merge(summary[["family", "modality", "strategy"]])
        )

    if not selected:
        return pd.DataFrame(columns=rows.columns)
    return pd.concat(selected, ignore_index=True)


def adc_fa_wscore_csv(disease: str) -> Path:
    return REPO_ROOT / "results" / f"wscore_candidate_benchmark_{disease.lower()}_adc_fa" / "feature_auc.csv"


def adc_fa_rows(disease: str) -> pd.DataFrame:
    path = adc_fa_wscore_csv(disease)
    if not path.exists():
        return pd.DataFrame()
    rows = pd.read_csv(path)
    rows = rows[rows["dataset"].isin(["MICs", "NOEL"]) & rows["feature"].isin(["ADC", "FA"])].copy()
    rows["family"] = "tle_adc_fa_wscore"
    rows["source_family"] = "w-score"
    rows["modality"] = rows["feature"]
    rows["strategy"] = rows["method"]
    rows["primary_score"] = rows["auc"]
    rows["disease"] = disease
    if rows.empty:
        return rows

    selected = []
    for feature in ["ADC", "FA"]:
        feature_rows = rows[rows["feature"].eq(feature)]
        summary = (
            feature_rows.groupby(["family", "source_family", "modality", "strategy"], dropna=False)
            .agg(mean_auc=("primary_score", "mean"), available_datasets=("dataset", "nunique"))
            .reset_index()
            .sort_values(["mean_auc", "available_datasets"], ascending=[False, False])
            .head(TOP_N_WSCORE_PER_FEATURE)
        )
        selected.append(
            feature_rows.merge(summary[["family", "modality", "strategy"]])
        )
    return pd.concat(selected, ignore_index=True)


def panel_table(rows: pd.DataFrame) -> pd.DataFrame:
    if rows.empty:
        return pd.DataFrame()
    summary = (
        rows.groupby(["source_family", "modality", "strategy"], dropna=False)
        .agg(mean_auc=("primary_score", "mean"), available_datasets=("dataset", "nunique"))
        .reset_index()
    )
    pivot = rows.pivot_table(
        index=["source_family", "modality", "strategy"],
        columns="dataset",
        values="primary_score",
        aggfunc="mean",
    ).reset_index()
    table = summary.merge(pivot, on=["source_family", "modality", "strategy"], how="left")
    table["sort_modality"] = (
        table["modality"].map({"T1w": 0, "FLAIR": 1, "ADC": 0, "FA": 1}).fillna(2)
    )
    table = table.sort_values(
        ["sort_modality", "mean_auc"],
        ascending=[True, True],
    ).reset_index(drop=True)
    return table


def plot_panel(ax, rows: pd.DataFrame, title: str, description: str) -> None:
    table = panel_table(rows)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlim(*XLIM)
    ax.set_xlabel("Disease-specific AUC")
    ax.grid(axis="x", color="#d0d0d0", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    if table.empty:
        ax.text(0.5, 0.5, "No rows available", transform=ax.transAxes, ha="center", va="center")
        ax.set_yticks([])
        return

    y_positions = range(len(table))
    labels = [
        wrap_label(f"{row.modality}: {row.source_family} / {row.strategy}")
        for row in table.itertuples(index=False)
    ]
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_ylim(-0.8, len(table) - 0.2)
    ax.text(
        0.01,
        0.98,
        wrap_label(description, width=62),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        color="#444444",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 2.0},
    )

    for y, row in zip(y_positions, table.itertuples(index=False)):
        mics = getattr(row, "MICs", float("nan"))
        noel = getattr(row, "NOEL", float("nan"))
        if pd.notna(mics) and pd.notna(noel):
            ax.plot([mics, noel], [y, y], color="#8d8d8d", linewidth=1.2, zorder=1)
            delta = noel - mics
            ax.text(
                min(max(mics, noel) + 0.008, XLIM[1] - 0.005),
                y,
                f"Δ {delta:+.3f}",
                va="center",
                ha="left",
                fontsize=7,
                color="#333333",
            )
        elif pd.notna(mics):
            ax.text(
                min(mics + 0.008, XLIM[1] - 0.005),
                y,
                "MICs only",
                va="center",
                ha="left",
                fontsize=7,
                color="#555555",
            )
        elif pd.notna(noel):
            ax.text(
                min(noel + 0.008, XLIM[1] - 0.005),
                y,
                "NOEL only",
                va="center",
                ha="left",
                fontsize=7,
                color="#555555",
            )

        for dataset, value in (("MICs", mics), ("NOEL", noel)):
            if pd.notna(value):
                text_y_offset = 0.19 if dataset == "MICs" else -0.19
                text_va = "bottom" if dataset == "MICs" else "top"
                ax.scatter(
                    [value],
                    [y],
                    s=34,
                    color=DATASET_COLORS[dataset],
                    edgecolor="white",
                    linewidth=0.7,
                    zorder=2,
                )
                ax.text(
                    value,
                    y + text_y_offset,
                    f"{value:.3f}",
                    ha="center",
                    va=text_va,
                    fontsize=6.5,
                    color="#111111",
                )


def make_figure(disease: str) -> None:
    screening = pd.read_csv(SCREENING_CSV)
    panels = list(CATEGORIES)
    if adc_fa_wscore_csv(disease).exists():
        panels = panels + [
            {
                "key": "tle_adc_fa_wscore",
                "title": f"{disease} ADC / FA W-Score Post-Processing",
                "description": f"{disease}-specific ADC and FA candidate W-score transforms where maps are available.",
            }
        ]

    ncols = 2
    nrows = ceil(len(panels) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(18, 5.2 * nrows),
        constrained_layout=True,
    )
    axes = axes.ravel()

    plotted_rows = []
    for ax, panel in zip(axes, panels):
        if panel["key"] == "tle_adc_fa_wscore":
            rows = adc_fa_rows(disease)
        else:
            rows = category_rows(screening, disease, panel)
        if not rows.empty:
            rows = rows.copy()
            rows["processing_category"] = panel["key"]
            plotted_rows.append(rows)
        plot_panel(ax, rows, panel["title"], panel["description"])

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="white",
            markersize=8,
        )
        for color in DATASET_COLORS.values()
    ]
    for ax in axes[len(panels):]:
        ax.axis("off")
        ax.legend(
            handles,
            list(DATASET_COLORS.keys()),
            title="Dataset",
            loc="center",
            frameon=False,
            fontsize=11,
            title_fontsize=12,
        )
        ax.text(
            0.5,
            0.30,
            "Lines connect matched processing choices.\n"
            "Delta is NOEL AUC minus MICs AUC.\n"
            "Single markers indicate a missing dataset.",
            ha="center",
            va="center",
            fontsize=10,
            color="#333333",
        )

    if len(panels) == len(axes):
        axes[0].legend(
            handles,
            list(DATASET_COLORS.keys()),
            title="Dataset",
            loc="lower right",
            frameon=False,
            fontsize=8,
            title_fontsize=9,
        )

    fig.suptitle(
        f"{disease} Feature-Engineering Benchmark by Processing Category",
        fontsize=18,
        fontweight="bold",
    )

    output_base = BENCHMARK_DIR / f"processing_feature_auc_by_processing_category_{disease.lower()}_mics_vs_noel"
    fig.savefig(output_base.with_suffix(".png"), dpi=300)
    fig.savefig(output_base.with_suffix(".svg"))
    print(f"Wrote {output_base.with_suffix('.png')}")
    print(f"Wrote {output_base.with_suffix('.svg')}")

    if plotted_rows:
        plotted = pd.concat(plotted_rows, ignore_index=True)
        cols = [
            "processing_category",
            "dataset",
            "disease",
            "source_family",
            "family",
            "modality",
            "strategy",
            "primary_score",
            "metric_name",
        ]
        available = [col for col in cols if col in plotted.columns]
        plotted[available].to_csv(
            BENCHMARK_DIR / f"processing_feature_auc_by_processing_category_{disease.lower()}_plotted_rows.csv",
            index=False,
        )


def main() -> None:
    make_figure("FCD")
    make_figure("TLE")


if __name__ == "__main__":
    main()
