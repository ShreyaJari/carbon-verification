"""
stage4_comparison_and_uncertainty.py

Satellite MRV Toolkit for Carbon Removal Verification
Stage 4: Comparison -- independent satellite-based estimate vs.
registry-claimed sequestration, with uncertainty bounds.

Case study: Kachung Forest Project (KFP), CDM Project 4653, Dokolo
District, Uganda. See stage1/2/3 scripts for full background.

Purpose
--------
Convert Stage 2's independent GEDI-based AGB estimate into a total
CO2e figure comparable to KFP's registry-claimed issued CERs, using
documented, standard (not project-specific) conversion factors, and
report the comparison alongside every assumption that makes it
approximate rather than precise.

This script is pure computation on Stage 2's output CSV plus the
registry figures sourced from the CDM PDD -- no APIs, no new
credentials, nothing that can fail for environment reasons like
Stages 1-3 did.

Read before trusting the output: two things make this an approximate
comparison, not a precise reconciliation
------------------------------------------------------------------
1. STOCK vs. NET REMOVALS: GEDI gives a biomass STOCK during the
   2019-2024 measurement window. The registry's issued CERs are a NET
   CUMULATIVE REMOVAL figure (project scenario minus a counterfactual
   baseline scenario, ex-post verified for 2006-2020). Treating them
   as comparable assumes the pre-project baseline biomass was
   approximately zero (the PDD describes the land as degraded
   grassland/shrubland, which supports this, but it is still an
   assumption, not a measured baseline). There's also a real
   1-4 year gap between the registry claim's 2020 cutoff and GEDI's
   2019-2024 measurement window, during which further growth would
   have occurred.
2. CONVERSION FACTORS ARE IPCC DEFAULTS, NOT KFP'S OWN VALUES: the
   root-to-shoot ratio and any species-specific biomass expansion
   factors the PDD's own AR-AM0004 methodology actually used are not
   reproduced here (they weren't captured from the PDD excerpts
   available at Stage 1). Using IPCC 2006 defaults instead introduces
   a real, documented source of divergence from the registry's own
   numbers, separate from any true measurement discrepancy.

Requirements: none beyond the Python standard library.
"""

import csv
import math
from pathlib import Path

# ---------------------------------------------------------------------------
# 1. Registry-claimed figures (source: CDM Project 4653 monitoring and
#    verification reports, as summarized in Stage 1's docstring)
# ---------------------------------------------------------------------------
REGISTRY_CERS_PERIOD_1 = 30_492   # tCO2e, monitoring period 2006-2012
REGISTRY_CERS_PERIOD_2 = 314_672  # tCO2e, monitoring period 2013-2020
REGISTRY_CUMULATIVE_CERS_THROUGH_2020 = (
    REGISTRY_CERS_PERIOD_1 + REGISTRY_CERS_PERIOD_2
)

# ---------------------------------------------------------------------------
# 2. Conversion factors -- IPCC 2006 Guidelines defaults, NOT KFP's own
#    project-specific values (see module docstring, point 2)
# ---------------------------------------------------------------------------
# Root-to-shoot ratio (R): IPCC 2006 GL for National GHG Inventories,
# AFOLU Vol. 4, Table 4.4, tropical forest default. KFP's actual PDD
# methodology (AR-AM0004) may use a different, species-specific value
# for Pinus caribaea / Eucalyptus -- flagged explicitly, not assumed
# equivalent.
ROOT_TO_SHOOT_RATIO = 0.24

# Carbon fraction of dry biomass (CF): IPCC 2006 GL default, applies
# broadly across forest types, comparatively low-uncertainty input.
CARBON_FRACTION = 0.47

# CO2 to C molecular weight ratio (44.01 / 12.01)
CO2_TO_C_RATIO = 44.01 / 12.01

KFP_TOTAL_PROJECT_AREA_HA = 2_099  # from CDM PDD, Section G1.3

INPUT_CSV = (
    Path(__file__).resolve().parent.parent
    / "data" / "processed" / "agb_mgha_pooled_by_block.csv"
)
OUTPUT_CSV = (
    Path(__file__).resolve().parent.parent
    / "data" / "processed" / "comparison_table.csv"
)

BLOCK_III_ID = "III"


# ---------------------------------------------------------------------------
# 3. Load Stage 2 output and compute area-weighted project estimates
# ---------------------------------------------------------------------------
def load_block_agb_rows():
    with open(INPUT_CSV, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        rows = []
        for row in reader:
            rows.append({
                "block_id": row["block_id"],
                "block_area_ha": float(row["block_area_ha"]),
                "agbd_mgha_mean": (
                    float(row["agbd_mgha_mean"])
                    if row["agbd_mgha_mean"] not in ("", "None")
                    else None
                ),
                "agbd_mgha_stderr": (
                    float(row["agbd_mgha_stderr"])
                    if row["agbd_mgha_stderr"] not in ("", "None")
                    else None
                ),
            })
        return rows


def area_weighted_agb_with_uncertainty(rows, exclude_block_id=None):
    """Area-weighted mean AGB density (Mg/ha) across blocks, with
    combined standard error assuming block-level estimates are
    independent (a simplifying assumption -- adjacent blocks likely
    have some spatial autocorrelation in practice, which would make
    the true combined uncertainty somewhat larger than this)."""
    usable = [
        r for r in rows
        if r["agbd_mgha_mean"] is not None and r["block_id"] != exclude_block_id
    ]
    if not usable:
        return None, None, 0.0

    total_area = sum(r["block_area_ha"] for r in usable)
    weighted_mean = sum(
        r["agbd_mgha_mean"] * r["block_area_ha"] for r in usable
    ) / total_area

    # Combined variance: sum of (weight * stderr)^2 across blocks with
    # a defined stderr; blocks with stderr=None (n<=1) are excluded
    # from the uncertainty combination but still contribute to the
    # mean above -- flagged in printed output, not silently ignored.
    variance_terms = [
        (r["block_area_ha"] / total_area * r["agbd_mgha_stderr"]) ** 2
        for r in usable
        if r["agbd_mgha_stderr"] is not None
    ]
    combined_stderr = math.sqrt(sum(variance_terms)) if variance_terms else None

    return weighted_mean, combined_stderr, total_area


# ---------------------------------------------------------------------------
# 4. Convert AGB density -> total CO2e for the full project area
# ---------------------------------------------------------------------------
def agb_density_to_total_co2e(agbd_mgha_mean, agbd_mgha_stderr, area_ha):
    """Apply the AGB -> AGB+BGB -> carbon -> CO2e conversion chain and
    scale to the full project area, propagating relative uncertainty
    through each multiplicative step (valid since every step here is
    a simple multiplication, so relative uncertainties add in
    quadrature)."""
    agb_total_mg = agbd_mgha_mean * area_ha
    bgb_total_mg = agb_total_mg * ROOT_TO_SHOOT_RATIO
    total_biomass_mg = agb_total_mg + bgb_total_mg
    carbon_stock_mg = total_biomass_mg * CARBON_FRACTION
    co2e_total_mg = carbon_stock_mg * CO2_TO_C_RATIO

    if agbd_mgha_stderr is not None and agbd_mgha_mean:
        relative_stderr = agbd_mgha_stderr / agbd_mgha_mean
        co2e_stderr = co2e_total_mg * relative_stderr
    else:
        co2e_stderr = None

    return co2e_total_mg, co2e_stderr


# ---------------------------------------------------------------------------
# 5. Main pipeline
# ---------------------------------------------------------------------------
def run_stage4():
    if not INPUT_CSV.exists():
        print(f"Stage 2 output not found at {INPUT_CSV}. Run Stage 2 first.")
        return

    rows = load_block_agb_rows()

    scenarios = {
        "excluding_block_III": area_weighted_agb_with_uncertainty(
            rows, exclude_block_id=BLOCK_III_ID
        ),
        "including_all_blocks": area_weighted_agb_with_uncertainty(rows),
    }

    print("=" * 70)
    print("KFP Stage 4: independent estimate vs. registry-claimed CERs")
    print("=" * 70)
    print(f"\nRegistry-claimed cumulative issued CERs through 2020: "
          f"{REGISTRY_CUMULATIVE_CERS_THROUGH_2020:,} tCO2e")
    print(f"  (period 1, 2006-2012: {REGISTRY_CERS_PERIOD_1:,} tCO2e)")
    print(f"  (period 2, 2013-2020: {REGISTRY_CERS_PERIOD_2:,} tCO2e)")

    output_rows = []

    for scenario_name, (mean_agbd, stderr_agbd, area_covered_ha) in scenarios.items():
        if mean_agbd is None:
            print(f"\n{scenario_name}: no usable AGB data, skipping.")
            continue

        # Scale the area-weighted density to the FULL CDM-eligible
        # project area (2,099 ha), not just the area actually covered
        # by blocks with valid GEDI data -- this assumes spatial
        # homogeneity across the project, an explicit approximation
        # flagged here, not hidden.
        co2e_total, co2e_stderr = agb_density_to_total_co2e(
            mean_agbd, stderr_agbd, KFP_TOTAL_PROJECT_AREA_HA
        )

        difference = co2e_total - REGISTRY_CUMULATIVE_CERS_THROUGH_2020
        pct_difference = (
            100 * difference / REGISTRY_CUMULATIVE_CERS_THROUGH_2020
        )

        print(f"\n--- Scenario: {scenario_name} ---")
        print(f"  Area-weighted AGB density: {mean_agbd:.1f} +/- "
              f"{stderr_agbd if stderr_agbd is None else round(stderr_agbd, 1)} Mg/ha "
              f"(density measured over {area_covered_ha:.0f} ha of blocks, "
              f"scaled here to full {KFP_TOTAL_PROJECT_AREA_HA} ha project area)")
        print(f"  Independent CO2e stock estimate (2019-2024 GEDI window): "
              f"{co2e_total:,.0f} +/- "
              f"{co2e_stderr if co2e_stderr is None else round(co2e_stderr):,} tCO2e")
        print(f"  Registry-claimed cumulative removals (through 2020): "
              f"{REGISTRY_CUMULATIVE_CERS_THROUGH_2020:,} tCO2e")
        print(f"  Difference (independent - registry): {difference:,.0f} tCO2e "
              f"({pct_difference:+.1f}%)")

        output_rows.append({
            "scenario": scenario_name,
            "agbd_mgha_mean": round(mean_agbd, 2),
            "agbd_mgha_stderr": (
                round(stderr_agbd, 2) if stderr_agbd is not None else ""
            ),
            "area_covered_by_gedi_ha": round(area_covered_ha, 1),
            "area_scaled_to_ha": KFP_TOTAL_PROJECT_AREA_HA,
            "independent_co2e_estimate_tco2e": round(co2e_total),
            "independent_co2e_stderr_tco2e": (
                round(co2e_stderr) if co2e_stderr is not None else ""
            ),
            "registry_claimed_cumulative_tco2e": (
                REGISTRY_CUMULATIVE_CERS_THROUGH_2020
            ),
            "difference_tco2e": round(difference),
            "pct_difference": round(pct_difference, 1),
        })

    with open(OUTPUT_CSV, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=output_rows[0].keys())
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"\nWrote comparison table to {OUTPUT_CSV}")

    print("\n" + "=" * 70)
    print("How to interpret this comparison -- read before using it")
    print("=" * 70)
    print("1. This compares a GEDI-measured STOCK (2019-2024) against a")
    print("   registry NET CUMULATIVE REMOVALS figure (2006-2020). They")
    print("   are related but not identical quantities -- see module")
    print("   docstring point 1. Treat any agreement as supportive, not")
    print("   confirmatory, and any disagreement as partly attributable")
    print("   to this mismatch, not solely to a measurement error by")
    print("   either party.")
    print("2. The root-to-shoot ratio and carbon fraction used are IPCC")
    print("   defaults, not KFP's own project-specific values -- see")
    print("   module docstring point 2. This is a genuine, quantifiable")
    print("   source of divergence separate from the biomass measurement")
    print("   itself.")
    print("3. Compare the 'excluding_block_III' and 'including_all_blocks'")
    print("   scenarios -- if they bracket the registry figure similarly,")
    print("   that's a reasonably strong result. If registry sits far")
    print("   outside both, that's the headline finding to write up,")
    print("   not something to explain away.")


if __name__ == "__main__":
    run_stage4()