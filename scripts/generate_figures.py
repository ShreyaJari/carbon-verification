"""
generate_figures.py

Satellite MRV Toolkit for Carbon Removal Verification
Visualization: generates the core figures from Stage 2/3/4 output CSVs.

Produces four PNGs in outputs/figures/:
    1. agb_by_block.png       -- AGB density per block, Block III flagged
    2. ndvi_trend_by_block.png -- Sentinel-2 NDVI trend, 2019-2024
    3. xco2_regional_trend.png -- OCO-2/3 sampled regional XCO2 points
    4. comparison_vs_registry.png -- independent estimate vs. registry claim

Reads from data/processed/ (must have already run Stages 2-4 to
produce the input CSVs). Does not recompute anything -- purely
visualizes what those stages already produced.

Requirements: pip install matplotlib
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "outputs" / "figures"

BLOCK_III_ID = "III"


def read_csv_rows(path):
    with open(path, newline="") as csv_file:
        return list(csv.DictReader(csv_file))


def to_float_or_none(value):
    if value in ("", "None", None):
        return None
    return float(value)


# ---------------------------------------------------------------------------
# Figure 1: AGB density per block
# ---------------------------------------------------------------------------
def plot_agb_by_block():
    rows = read_csv_rows(DATA_DIR / "agb_mgha_pooled_by_block.csv")

    block_ids = [r["block_id"] for r in rows]
    means = [to_float_or_none(r["agbd_mgha_mean"]) or 0 for r in rows]
    stderrs = [to_float_or_none(r["agbd_mgha_stderr"]) or 0 for r in rows]
    colors = [
        "#c0392b" if r["block_id"] == BLOCK_III_ID else "#2e7d5b"
        for r in rows
    ]
    n_shots = [r["gedi_valid_pixel_count"] for r in rows]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(block_ids, means, yerr=stderrs, capsize=6, color=colors)

    for bar, stderr, n in zip(bars, stderrs, n_shots):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + stderr + max(means) * 0.03,
            f"n={n}",
            ha="center", fontsize=9,
        )

    ax.set_xlabel("KFP block")
    ax.set_ylabel("Aboveground biomass density (Mg/ha)")
    ax.set_title(
        "GEDI L4A pooled AGB density by KFP block (2019-2024)\n"
        "Red = Block III, flagged as a likely geometry-contamination "
        "outlier (see README)"
    )
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "agb_by_block.png", dpi=150)
    plt.close(fig)
    print("Wrote agb_by_block.png")


# ---------------------------------------------------------------------------
# Figure 2: Sentinel-2 NDVI trend per block
# ---------------------------------------------------------------------------
def plot_ndvi_trend():
    rows = read_csv_rows(DATA_DIR / "sentinel2_ndvi_by_block_year.csv")

    by_block = {}
    for row in rows:
        block_id = row["block_id"]
        year = int(row["year"])
        ndvi = to_float_or_none(row["sentinel2_ndvi_mean"])
        by_block.setdefault(block_id, []).append((year, ndvi))

    fig, ax = plt.subplots(figsize=(8, 5))
    for block_id, points in sorted(by_block.items()):
        points.sort(key=lambda p: p[0])
        years = [p[0] for p in points]
        ndvi_values = [p[1] for p in points]
        linestyle = "--" if block_id == BLOCK_III_ID else "-"
        ax.plot(years, ndvi_values, marker="o", linestyle=linestyle,
                 label=f"Block {block_id}")

    ax.set_xlabel("Year")
    ax.set_ylabel("Median annual NDVI")
    ax.set_title(
        "Sentinel-2 NDVI trend per KFP block (2019-2024)\n"
        "Context only -- not a biomass estimate"
    )
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "ndvi_trend_by_block.png", dpi=150)
    plt.close(fig)
    print("Wrote ndvi_trend_by_block.png")


# ---------------------------------------------------------------------------
# Figure 3: Regional XCO2 trend (OCO-2/3)
# ---------------------------------------------------------------------------
def plot_xco2_trend():
    rows = read_csv_rows(DATA_DIR / "xco2_ppm_regional.csv")

    mission_colors = {"OCO-2": "#2166ac", "OCO-3": "#b2182b"}
    fig, ax = plt.subplots(figsize=(8, 5))

    for mission, color in mission_colors.items():
        mission_rows = [r for r in rows if r["mission"] == mission]
        mission_rows.sort(key=lambda r: r["sample_date"])
        dates = [r["sample_date"] for r in mission_rows]
        means = [to_float_or_none(r["xco2_ppm_mean"]) for r in mission_rows]
        stderrs = [
            to_float_or_none(r["xco2_ppm_stderr"]) or 0
            for r in mission_rows
        ]
        n_soundings = [r["n_soundings"] for r in mission_rows]

        if not dates:
            continue

        ax.errorbar(
            dates, means, yerr=stderrs, marker="o", linestyle="none",
            color=color, label=mission, capsize=4,
        )
        for x, y, n in zip(dates, means, n_soundings):
            ax.annotate(f"n={n}", (x, y), textcoords="offset points",
                        xytext=(0, 8), fontsize=7, ha="center")

    ax.set_xlabel("Sample date")
    ax.set_ylabel("Regional mean XCO2 (ppm)")
    ax.set_title(
        "Regional XCO2, ~110 km box around KFP (sparse sampling)\n"
        "NOT a project-scale signal -- see README for why"
    )
    ax.tick_params(axis="x", rotation=45)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "xco2_regional_trend.png", dpi=150)
    plt.close(fig)
    print("Wrote xco2_regional_trend.png")


# ---------------------------------------------------------------------------
# Figure 4: Independent estimate vs. registry claim (the headline chart)
# ---------------------------------------------------------------------------
def plot_comparison_vs_registry():
    rows = read_csv_rows(DATA_DIR / "comparison_table.csv")

    labels = []
    values = []
    errors = []
    colors = []

    for row in rows:
        scenario_key = row["scenario"]
        if scenario_key == "excluding_block_III":
            label = "Excluding\nBlock III"
        elif scenario_key == "including_all_blocks":
            label = "Including\nAll Blocks"
        else:
            label = scenario_key.replace("_", " ").title()
        labels.append(label)
        values.append(float(row["independent_co2e_estimate_tco2e"]))
        errors.append(
            float(row["independent_co2e_stderr_tco2e"])
            if row["independent_co2e_stderr_tco2e"] not in ("", "None")
            else 0
        )
        colors.append("#2e7d5b")

    registry_value = float(rows[0]["registry_claimed_cumulative_tco2e"])
    labels.append("Registry Claimed\n(issued CERs)")
    values.append(registry_value)
    errors.append(0)
    colors.append("#555555")

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(labels, values, yerr=errors, capsize=6, color=colors)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(values) * 0.02,
            f"{value:,.0f}",
            ha="center", fontsize=9,
        )

    ax.set_ylabel("Cumulative CO2e (tCO2e)")
    ax.set_title(
        "Independent satellite-based estimate vs. KFP registry-claimed\n"
        "cumulative removals through 2020"
    )
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "comparison_vs_registry.png", dpi=150)
    plt.close(fig)
    print("Wrote comparison_vs_registry.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    missing = [
        name for name in [
            "agb_mgha_pooled_by_block.csv",
            "sentinel2_ndvi_by_block_year.csv",
            "xco2_ppm_regional.csv",
            "comparison_table.csv",
        ]
        if not (DATA_DIR / name).exists()
    ]
    if missing:
        print(f"Missing required input file(s): {missing}. "
              f"Run Stages 2-4 first.")
        return

    plot_agb_by_block()
    plot_ndvi_trend()
    plot_xco2_trend()
    plot_comparison_vs_registry()

    print(f"\nAll figures written to {FIGURES_DIR}")


if __name__ == "__main__":
    main()