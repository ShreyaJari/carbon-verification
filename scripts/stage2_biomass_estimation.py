"""
stage2_biomass_estimation.py

Satellite MRV Toolkit for Carbon Removal Verification
Stage 2: Biomass stage (bottom-up independent estimate)

Case study: Kachung Forest Project (KFP), CDM Project 4653, Dokolo
District, Uganda. See stage1_data_access_verification.py for project
background and registry-claimed sequestration figures.

Revision history
------------------
v1: Queried LARSE/GEDI/GEDI04_A_002_INDEX as if it held shot-level
    AGBD values. That collection is only an INDEX (table_id,
    time_start, time_end pointers to per-orbit tables) -- filtering it
    on agbd/quality properties silently returned zero matches. Fixed
    by switching to the LARSE/GEDI/GEDI04_A_002_MONTHLY raster
    collection, which rasterizes shot AGBD to footprint locations.

v2 (this version): The v1 fix ran successfully but revealed GEDI shot
    density per block PER YEAR is too sparse to support an annual
    breakdown -- most block-years had 0-5 valid pixels, several with
    n=1 (a meaningless stderr=0.0). This version instead pools shots
    across the full available GEDI mission span per block, for a
    single defensible mean with a real sample size. It also flags
    Block III specifically: a v1 run showed an anomalous 457 Mg/ha
    reading there (n=5), plausibly because the block's circular
    area-matched approximation clips into a remnant government
    Gmelina arborea plantation that the PDD (Section G1.2) explicitly
    describes as excluded from the CDM-eligible area for predating the
    project. This version reports Block III's pooled estimate
    alongside a project-wide total computed BOTH with and without
    Block III included, rather than silently folding a possibly
    contaminated reading into one number.

Two documented approximations carried over from v1 (read before
trusting the output)
------------------------------------------------------------------
1. PROJECT BOUNDARY: block geometries are circles of matching area
   centered on the PDD's block centroids, NOT the true surveyed
   polygon (the PDD only gives centroids + areas, not vertices). This
   is very likely the mechanism behind the Block III anomaly above.
2. GEDI L4A BIOMASS MODEL: GEDI L4A's predictive AGBD model depends on
   a global plant-functional-type (PFT) map, which may not correctly
   classify a pine/eucalyptus plantation established on former
   savanna/shrubland. Not yet surfaced per-pixel in this output.

Requirements: same environment as Stage 1 (earthengine-api,
authenticated with project="carbon-verification-toolkit"). No new
packages needed.

Run this locally, paste the output back before Stage 3.
"""

import csv
import math
import sys
from pathlib import Path

import ee

# ---------------------------------------------------------------------------
# 1. Project setup
# ---------------------------------------------------------------------------
EE_PROJECT_ID = "carbon-verification-toolkit"

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
AGB_OUTPUT_CSV = OUTPUT_DIR / "agb_mgha_pooled_by_block.csv"
NDVI_OUTPUT_CSV = OUTPUT_DIR / "sentinel2_ndvi_by_block_year.csv"

# GEDI L4A quality thresholds (GEDI L4A user guide recommended values).
GEDI_MIN_SENSITIVITY = 0.95
GEDI_REQUIRE_QUALITY_FLAG = 1
GEDI_REQUIRE_DEGRADE_FLAG = 0

# Pool across the full available GEDI L4A mission span rather than
# breaking out by year (see revision history above for why). Catalog
# entry for GEDI04_A_002_MONTHLY lists availability through Nov 2024.
GEDI_POOL_START_DATE = "2019-01-01"
GEDI_POOL_END_DATE = "2024-11-30"

# Sentinel-2 NDVI trend is still broken out annually -- that data was
# dense throughout in v1, no sparsity issue there.
NDVI_ANALYSIS_YEARS = list(range(2019, 2025))

BLOCK_III_ID = "III"
BLOCK_III_FLAG_NOTE = (
    "Anomalously high AGBD in prior per-year runs (457 Mg/ha at n=5); "
    "circular area-matched geometry for this block plausibly clips into "
    "a remnant, non-CDM-eligible Gmelina arborea plantation described in "
    "PDD Section G1.2. Treat this block's estimate with caution; "
    "project-wide totals are reported both with and without it."
)

# ---------------------------------------------------------------------------
# 2. KFP block definitions
# ---------------------------------------------------------------------------
# Source: CCBA/CDM PDD, Table G.1.3.1. UTM Zone 36N. Centroids only,
# not full polygon vertices -- see module docstring, approximation #1.
KFP_BLOCKS = [
    {"block_id": "I", "area_ha": 275.2, "utm_easting": 491742, "utm_northing": 222989},
    {"block_id": "II", "area_ha": 815.5, "utm_easting": 493280, "utm_northing": 220742},
    {"block_id": "III", "area_ha": 121.9, "utm_easting": 493490, "utm_northing": 223484},
    {"block_id": "IV", "area_ha": 228.1, "utm_easting": 495325, "utm_northing": 222157},
    {"block_id": "V", "area_ha": 658.2, "utm_easting": 497792, "utm_northing": 222592},
]
KFP_UTM_EPSG = "EPSG:32636"  # UTM Zone 36N, WGS84

TOTAL_BLOCK_AREA_HA = sum(b["area_ha"] for b in KFP_BLOCKS)
print(f"Sum of block areas: {TOTAL_BLOCK_AREA_HA:.1f} ha "
      f"(PDD states 2,099 ha total CDM-eligible area -- should match "
      f"closely as a sanity check).")


def build_block_geometry(block):
    """Circular area-matched approximation of a KFP block. See module
    docstring, approximation #1, for the documented limitation."""
    area_m2 = block["area_ha"] * 10_000
    radius_m = math.sqrt(area_m2 / math.pi)
    centroid = ee.Geometry.Point(
        [block["utm_easting"], block["utm_northing"]], KFP_UTM_EPSG
    )
    return centroid.buffer(radius_m)


# ---------------------------------------------------------------------------
# 3. GEDI L4A extraction -- pooled across the full mission span
# ---------------------------------------------------------------------------
def get_gedi_agb_pooled(block_geometry):
    """Quality-masked mean AGBD (Mg/ha) over a block's geometry, pooled
    across the full GEDI L4A mission span rather than broken out by
    year. Returns mean, stderr, and valid-pixel count (a proxy for
    quality-passing GEDI shot count within the block)."""
    monthly_images = (
        ee.ImageCollection("LARSE/GEDI/GEDI04_A_002_MONTHLY")
        .filterBounds(block_geometry)
        .filterDate(GEDI_POOL_START_DATE, GEDI_POOL_END_DATE)
    )

    if monthly_images.size().getInfo() == 0:
        return {"shot_count": 0, "agbd_mgha_mean": None, "agbd_mgha_stderr": None}

    def quality_mask(image):
        return (
            image.updateMask(
                image.select("l4_quality_flag").eq(GEDI_REQUIRE_QUALITY_FLAG)
            )
            .updateMask(
                image.select("degrade_flag").eq(GEDI_REQUIRE_DEGRADE_FLAG)
            )
            .updateMask(
                image.select("sensitivity").gte(GEDI_MIN_SENSITIVITY)
            )
        )

    masked_agbd = monthly_images.map(quality_mask).select("agbd").mosaic()

    stats = masked_agbd.reduceRegion(
        reducer=ee.Reducer.mean()
        .combine(ee.Reducer.stdDev(), sharedInputs=True)
        .combine(ee.Reducer.count(), sharedInputs=True),
        geometry=block_geometry,
        scale=25,  # native GEDI footprint diameter
        maxPixels=1e9,
    ).getInfo()

    valid_pixel_count = stats.get("agbd_count") or 0
    mean_agbd = stats.get("agbd_mean")
    std_agbd = stats.get("agbd_stdDev")
    stderr_agbd = (
        std_agbd / math.sqrt(valid_pixel_count)
        if std_agbd is not None and valid_pixel_count > 0
        else None
    )

    return {
        "shot_count": valid_pixel_count,
        "agbd_mgha_mean": mean_agbd,
        "agbd_mgha_stderr": stderr_agbd,
    }


# ---------------------------------------------------------------------------
# 4. Sentinel-2 NDVI context (unchanged approach, annual, dense data)
# ---------------------------------------------------------------------------
def get_sentinel2_ndvi_for_block_year(block_geometry, year):
    """Median annual NDVI over a block for canopy-development context.
    NOT a biomass estimate -- corroborating trend evidence only."""
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(block_geometry)
        .filterDate(f"{year}-01-01", f"{year}-12-31")
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
    )
    if s2.size().getInfo() == 0:
        return None

    def add_ndvi(image):
        return image.addBands(
            image.normalizedDifference(["B8", "B4"]).rename("ndvi")
        )

    median_ndvi = s2.map(add_ndvi).select("ndvi").median()
    ndvi_stats = median_ndvi.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=block_geometry,
        scale=20,
        maxPixels=1e9,
    ).getInfo()
    return ndvi_stats.get("ndvi")


# ---------------------------------------------------------------------------
# 5. Main pipeline
# ---------------------------------------------------------------------------
def run_stage2():
    try:
        ee.Initialize(project=EE_PROJECT_ID)
    except Exception as init_error:
        print(f"ee.Initialize(project='{EE_PROJECT_ID}') failed: "
              f"{init_error}")
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    agb_rows = []
    ndvi_rows = []
    block_geometries = {}

    print("\n" + "-" * 70)
    print("Pooled GEDI L4A biomass estimate per block "
          f"({GEDI_POOL_START_DATE} to {GEDI_POOL_END_DATE})")
    print("-" * 70)
    for block in KFP_BLOCKS:
        block_geometry = build_block_geometry(block)
        block_geometries[block["block_id"]] = block_geometry

        gedi_result = get_gedi_agb_pooled(block_geometry)
        flag_note = BLOCK_III_FLAG_NOTE if block["block_id"] == BLOCK_III_ID else ""

        print(
            f"Block {block['block_id']} ({block['area_ha']} ha): "
            f"n={gedi_result['shot_count']}, "
            f"agbd_mgha_mean={gedi_result['agbd_mgha_mean']}, "
            f"agbd_mgha_stderr={gedi_result['agbd_mgha_stderr']}"
            + (f"  [FLAGGED: {flag_note}]" if flag_note else "")
        )

        agb_rows.append({
            "block_id": block["block_id"],
            "block_area_ha": block["area_ha"],
            "gedi_valid_pixel_count": gedi_result["shot_count"],
            "agbd_mgha_mean": gedi_result["agbd_mgha_mean"],
            "agbd_mgha_stderr": gedi_result["agbd_mgha_stderr"],
            "flag_note": flag_note,
        })

    # --- Project-wide area-weighted mean, with and without Block III ---
    def area_weighted_mean(rows, exclude_block_id=None):
        usable = [
            r for r in rows
            if r["agbd_mgha_mean"] is not None
            and r["block_id"] != exclude_block_id
        ]
        if not usable:
            return None, 0.0
        total_area = sum(r["block_area_ha"] for r in usable)
        weighted_sum = sum(
            r["agbd_mgha_mean"] * r["block_area_ha"] for r in usable
        )
        return weighted_sum / total_area, total_area

    mean_with_all, area_with_all = area_weighted_mean(agb_rows)
    mean_excl_iii, area_excl_iii = area_weighted_mean(
        agb_rows, exclude_block_id=BLOCK_III_ID
    )

    print("\n" + "-" * 70)
    print("Project-wide area-weighted AGB estimate")
    print("-" * 70)
    print(f"Including all 5 blocks ({area_with_all:.1f} ha covered): "
          f"{mean_with_all}")
    print(f"Excluding Block III ({area_excl_iii:.1f} ha covered): "
          f"{mean_excl_iii}")
    print("Compare these two -- if they differ substantially, that's")
    print("direct evidence the Block III geometry issue matters and")
    print("needs fixing (tighter boundary) before Stage 4, not just")
    print("footnoting.")

    # --- Sentinel-2 NDVI trend, annual, per block (context only) ---
    print("\n" + "-" * 70)
    print("Sentinel-2 NDVI annual trend per block (context, not biomass)")
    print("-" * 70)
    for block in KFP_BLOCKS:
        block_geometry = block_geometries[block["block_id"]]
        for year in NDVI_ANALYSIS_YEARS:
            ndvi_mean = get_sentinel2_ndvi_for_block_year(block_geometry, year)
            print(f"  Block {block['block_id']}, {year}: ndvi_mean={ndvi_mean}")
            ndvi_rows.append({
                "block_id": block["block_id"],
                "year": year,
                "sentinel2_ndvi_mean": ndvi_mean,
            })

    # --- Write outputs ---
    with open(AGB_OUTPUT_CSV, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=agb_rows[0].keys())
        writer.writeheader()
        writer.writerows(agb_rows)
    print(f"\nWrote {len(agb_rows)} pooled block AGB rows to {AGB_OUTPUT_CSV}")

    with open(NDVI_OUTPUT_CSV, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=ndvi_rows[0].keys())
        writer.writeheader()
        writer.writerows(ndvi_rows)
    print(f"Wrote {len(ndvi_rows)} block-year NDVI rows to {NDVI_OUTPUT_CSV}")

    print("\nReminders before trusting this output:")
    print("  - Block geometries remain circular area-matched")
    print("    approximations, not the true surveyed boundary.")
    print("  - Even pooled, check each block's n before treating its")
    print("    mean as reliable -- pooling raises n but doesn't")
    print("    guarantee it's now large enough for every block.")
    print("  - GEDI's per-pixel PFT model choice is still not surfaced")
    print("    here; flag if any block's AGBD looks implausible for a")
    print("    young pine/eucalyptus plantation.")


if __name__ == "__main__":
    run_stage2()