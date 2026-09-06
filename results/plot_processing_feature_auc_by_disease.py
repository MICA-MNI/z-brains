#!/usr/bin/env python3
"""Plot disease-specific AUC values with MICs-vs-NOEL differences."""

from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_DIR = REPO_ROOT / "results" / "processing_method_benchmark"
SCREENING_CSV = BENCHMARK_DIR / "screening_long.csv"
TLE_ADC_FA_CSV = REPO_ROOT / "results" / "wscore_candidate_benchmark_tle_adc_fa" / "feature_auc.csv"

TOP_N_PROCESSING_PER_MODALITY = 4
TOP_N_WSCORE_PER_FEATURE = 4
XLIM = (0.35, 0.90)

DATASET_COLORS = {
    "MICs": "#3f77b5",
    "NOEL": "#c26b3f",
}

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
    "tle_adc_fa_wscore": "TLE ADC / FA W-Score Candidates",
}


def wrap_label(text: str, width: int = 37) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def select_processing_rows(screening: pd.DataFrame, disease: str, family: str) -> pd.DataFrame:
    rows = screening[
        screening["axis"].eq("feature_processing")
        & screening["metric_name"].eq("mean_subject_auc_by_disease")
        & screening["disease"].eq(disease)
        & screening["family"].eq(family)
        & screening["dataset"].isin(["MICs", "NOEL"])
    ].copy()
    if rows.empty:
        return rows

    selected = []
    for modality in ["T1w", "FLAIR"]:
        modality_rows = rows[rows["modality"].eq(modality)]
        if modality_rows.empty:
            continue
        summary = (
            modality_rows.groupby(["family", "modality", "strategy"], dropna=False)
            .agg(mean_auc=("primary_score", "mean"), available_datasets=("dataset", "nunique"))
            .reset_index()
            .sort_values(["mean_auc", "available_datasets"], ascending=[False, False])
            .head(TOP_N_PROCESSING_PER_MODALITY)
        )
        selected.append(
            modality_rows.merge(summary[["family", "modality", "strategy"]])
        )

    if not selected:
        return pd.DataFrame(columns=rows.columns)
    return pd.concat(selected, ignore_index=True)


def select_tle_adc_fa_rows() -> pd.DataFrame:
    rows = pd.read_csv(TLE_ADC_FA_CSV)
    rows = rows[rows["dataset"].isin(["MICs", "NOEL"]) & rows["feature"].isin(["ADC", "FA"])].copy()
    rows["family"] = "tle_adc_fa_wscore"
    rows["modality"] = rows["feature"]
    rows["strategy"] = rows["method"]
    rows["primary_score"] = rows["auc"]
    rows["metric_name"] = "feature_auc_by_tle"
    rows["disease"] = "TLE"

    selected = []
    for feature in ["ADC", "FA"]:
        feature_rows = rows[rows["feature"].eq(feature)]
        summary = (
            feature_rows.groupby(["family", "modality", "strategy"], dropna=False)
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
        rows.groupby(["modality", "strategy"], dropna=False)
        .agg(mean_auc=("primary_score", "mean"), available_datasets=("dataset", "nunique"))
        .reset_index()
    )
    pivot = rows.pivot_table(
        index=["modality", "strategy"],
        columns="dataset",
        values="primary_score",
        aggfunc="mean",
    ).reset_index()
    table = summary.merge(pivot, on=["modality", "strategy"], how="left")
    table["sort_modality"] = table["modality"].map({"T1w": 0, "FLAIR": 1, "ADC": 0, "FA": 1}).fillna(2)
    table = table.sort_values(
        ["sort_modality", "mean_auc"],
        ascending=[True, True],
    ).reset_index(drop=True)
    return table


def plot_panel(ax, rows: pd.DataFrame, title: str) -> None:
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
        wrap_label(f"{row.modality}: {row.strategy}")
        for row in table.itertuples(index=False)
    ]
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_ylim(-0.6, len(table) - 0.4)

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
    families = FAMILY_ORDER.copy()
    if disease == "TLE":
        families.append("tle_adc_fa_wscore")

    fig, axes = plt.subplots(4, 2, figsize=(18, 22), constrained_layout=True)
    axes = axes.ravel()

    for ax, family in zip(axes, families):
        if family == "tle_adc_fa_wscore":
            rows = select_tle_adc_fa_rows()
        else:
            rows = select_processing_rows(screening, disease, family)
        plot_panel(ax, rows, FAMILY_LABELS.get(family, family))

    for ax in axes[len(families):]:
        ax.axis("off")
        handles = [
            plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="white", markersize=8)
            for color in DATASET_COLORS.values()
        ]
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
            "Lines connect matched methods.\n"
            "Delta is NOEL AUC minus MICs AUC.\n"
            "Single markers indicate a missing dataset.",
            ha="center",
            va="center",
            fontsize=10,
            color="#333333",
        )

    if disease == "TLE":
        handles = [
            plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="white", markersize=8)
            for color in DATASET_COLORS.values()
        ]
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
        f"{disease} Processing Feature Benchmark: MICs vs NOEL AUC Differences",
        fontsize=18,
        fontweight="bold",
    )

    output_base = BENCHMARK_DIR / f"processing_feature_auc_by_category_{disease.lower()}_mics_vs_noel"
    fig.savefig(output_base.with_suffix(".png"), dpi=300)
    fig.savefig(output_base.with_suffix(".svg"))
    print(f"Wrote {output_base.with_suffix('.png')}")
    print(f"Wrote {output_base.with_suffix('.svg')}")


def main() -> None:
    make_figure("FCD")
    make_figure("TLE")


if __name__ == "__main__":
    main()
