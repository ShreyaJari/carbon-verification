"""
stage3_atmospheric_crosscheck.py

Satellite MRV Toolkit for Carbon Removal Verification
Stage 3: Atmospheric stage (top-down, honestly scoped regional
cross-check)

Case study: Kachung Forest Project (KFP), CDM Project 4653, Dokolo
District, Uganda. See stage1_data_access_verification.py and
stage2_biomass_estimation.py for project background and the biomass
stage's results.

Purpose and scope -- read this before interpreting any output
------------------------------------------------------------------
KFP is ~2,099 ha. OCO-2/OCO-3 footprints are ~1.29 km x 2.25 km with
sparse orbital tracks -- a project this size will essentially never be
resolved as a distinct atmospheric CO2 drawdown signal against
background variability and instrument noise. This script does NOT
attempt project-scale detection. It computes a REGIONAL (roughly
+/-1 degree, ~110 km box around KFP) mean XCO2 time series as a coarse
sanity check only: does regional atmospheric CO2 look broadly
consistent with expected trends, with no attempt to attribute any
year-to-year wiggle in this regional mean to KFP specifically. Treat
any local minimum or maximum in the printed series as noise unless
proven otherwise -- this script does not perform that proof.

Sampling approach
--------------------
OCO-2/OCO-3 "Lite FP" granules are DAILY GLOBAL files (not small
regional extracts) -- downloading every available day across a
2014-2024 span would mean thousands of files, impractical and
unnecessary for a coarse regional check. This script instead samples
a fixed number of representative days per year (spread across
March/July/November to avoid seasonal bias), searches for the nearest
available granule to each target date, and downloads only those. This
is a deliberate scope decision, not a hidden shortcut -- see the
module docstring above for why project-scale exhaustive sampling
wouldn't be meaningful anyway.

Requirements (install locally, in addition to Stage 1/2's packages)
------------------------------------------------------------------
    pip install xarray h5netcdf

Run this locally (same earthaccess login as Stage 1), paste the
output back before Stage 4.
"""

import csv
import sys
from datetime import date, timedelta
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# 1. Region and sampling definition
# ---------------------------------------------------------------------------
# Same KFP bounding box as Stage 1, padded to a regional box -- see
# that script's module docstring for the source of the KFP coordinates.
KFP_LAT_MIN_DEG = 1.98222
KFP_LAT_MAX_DEG = 2.04222
KFP_LON_MIN_DEG = 32.91528
KFP_LON_MAX_DEG = 32.99528

REGIONAL_PAD_DEG = 1.0
REGIONAL_LON_MIN = KFP_LON_MIN_DEG - REGIONAL_PAD_DEG
REGIONAL_LON_MAX = KFP_LON_MAX_DEG + REGIONAL_PAD_DEG
REGIONAL_LAT_MIN = KFP_LAT_MIN_DEG - REGIONAL_PAD_DEG
REGIONAL_LAT_MAX = KFP_LAT_MAX_DEG + REGIONAL_PAD_DEG

# Representative sample dates: 3 per year (spread across seasons),
# for each year both missions could plausibly have data.
SAMPLE_MONTHS_DAYS = [(3, 15), (7, 15), (11, 15)]
OCO2_SAMPLE_YEARS = list(range(2015, 2025))  # OCO-2 launched July 2014
OCO3_SAMPLE_YEARS = list(range(2019, 2025))  # OCO-3 installed May 2019

XCO2_GOOD_QUALITY_FLAG = 0  # per OCO-2/3 Lite FP product convention

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw" / "oco2_oco3"
OUTPUT_CSV = OUTPUT_DIR / "xco2_ppm_regional.csv"


def build_sample_dates(years):
    """Return a list of date objects: one per (year, month, day) in the
    sampling grid, for the given years."""
    sample_dates = []
    for year in years:
        for month, day in SAMPLE_MONTHS_DAYS:
            try:
                sample_dates.append(date(year, month, day))
            except ValueError:
                continue
    return sample_dates


# ---------------------------------------------------------------------------
# 2. Search and download sampled granules
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 2. Search and download sampled granules
# ---------------------------------------------------------------------------
# NOTE on a corrected approach: OCO-2/3 Lite FP files are DAILY GLOBAL
# granules -- a file's overall bounding box touches nearly every region
# on Earth regardless of that day's actual orbit track, since OCO's
# swath is only ~10 km wide with a 16-day repeat cycle. Searching by
# bounding box and taking the single nearest-date match (as an earlier
# version of this script did) reliably "finds" a granule that then
# turns out to have ZERO soundings actually near KFP -- confirmed by
# inspecting a downloaded file directly (global lat/lon extent, 0
# soundings in the regional box). This version instead searches a wide
# window covering multiple repeat cycles, downloads and checks several
# candidate days in sequence, and keeps the first one with real
# in-region coverage rather than trusting the first search match.
SEARCH_WINDOW_HALF_DAYS = 30  # ~2 full OCO-2 repeat cycles either side
MAX_CANDIDATES_PER_PERIOD = 8


def find_and_download_granule_near_date(short_name, target_date, download_dir):
    """Search a wide window around target_date, download candidate
    granules one at a time, and return the local path of the first one
    that actually has soundings in the regional box (checked via
    extract_regional_xco2 by the caller) -- or None if none of the
    candidates checked have real coverage."""
    import earthaccess

    window_start = (target_date - timedelta(days=SEARCH_WINDOW_HALF_DAYS)).isoformat()
    window_end = (target_date + timedelta(days=SEARCH_WINDOW_HALF_DAYS)).isoformat()

    results = earthaccess.search_data(
        short_name=short_name,
        bounding_box=(
            REGIONAL_LON_MIN,
            REGIONAL_LAT_MIN,
            REGIONAL_LON_MAX,
            REGIONAL_LAT_MAX,
        ),
        temporal=(window_start, window_end),
        count=MAX_CANDIDATES_PER_PERIOD,
    )
    if not results:
        return []

    download_dir.mkdir(parents=True, exist_ok=True)
    downloaded_paths = earthaccess.download(results, str(download_dir))
    return [p for p in downloaded_paths if p]


# ---------------------------------------------------------------------------
# 3. Extract quality-filtered XCO2 soundings from a granule
# ---------------------------------------------------------------------------
def extract_regional_xco2(granule_path):
    """Open an OCO-2/3 Lite FP granule, filter soundings to the
    regional box and good quality flag, return (mean, stderr, n) for
    xco2 in ppm, or (None, None, 0) if nothing passes."""
    import xarray as xr

    with xr.open_dataset(granule_path) as ds:
        lat = ds["latitude"].values
        lon = ds["longitude"].values
        xco2 = ds["xco2"].values
        quality_flag = ds["xco2_quality_flag"].values

    in_region = (
        (lat >= REGIONAL_LAT_MIN)
        & (lat <= REGIONAL_LAT_MAX)
        & (lon >= REGIONAL_LON_MIN)
        & (lon <= REGIONAL_LON_MAX)
    )
    good_quality = quality_flag == XCO2_GOOD_QUALITY_FLAG
    keep = in_region & good_quality

    n_soundings = int(np.sum(keep))
    if n_soundings == 0:
        return None, None, 0

    values = xco2[keep]
    mean_xco2 = float(np.mean(values))
    stderr_xco2 = (
        float(np.std(values, ddof=1) / np.sqrt(n_soundings))
        if n_soundings > 1
        else None
    )
    return mean_xco2, stderr_xco2, n_soundings


# ---------------------------------------------------------------------------
# 4. Main pipeline
# ---------------------------------------------------------------------------
def run_stage3():
    try:
        import earthaccess
    except ImportError:
        print("earthaccess is not installed. Run: pip install earthaccess")
        sys.exit(1)
    try:
        import xarray  # noqa: F401
    except ImportError:
        print("xarray is not installed. Run: pip install xarray h5netcdf")
        sys.exit(1)

    try:
        earthaccess.login()
    except Exception as login_error:
        print(f"earthaccess.login() failed: {login_error}")
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    missions = [
        ("OCO2_L2_Lite_FP", "OCO-2", build_sample_dates(OCO2_SAMPLE_YEARS)),
        ("OCO3_L2_Lite_FP", "OCO-3", build_sample_dates(OCO3_SAMPLE_YEARS)),
    ]

    for short_name, mission_label, sample_dates in missions:
        print(f"\n{'-' * 70}")
        print(f"{mission_label}: sampling {len(sample_dates)} target dates")
        print(f"{'-' * 70}")

        for target_date in sample_dates:
            candidate_paths = find_and_download_granule_near_date(
                short_name, target_date, RAW_DIR / mission_label.lower()
            )
            if not candidate_paths:
                print(f"  {target_date}: no candidate granules found in "
                      f"+/-{SEARCH_WINDOW_HALF_DAYS} day window")
                continue

            found_coverage = False
            for candidate_path in candidate_paths:
                mean_xco2, stderr_xco2, n_soundings = extract_regional_xco2(
                    candidate_path
                )
                if n_soundings and n_soundings > 0:
                    print(
                        f"  {target_date}: used {Path(candidate_path).name}, "
                        f"n={n_soundings}, xco2_ppm_mean={mean_xco2}, "
                        f"xco2_ppm_stderr={stderr_xco2}"
                    )
                    rows.append({
                        "mission": mission_label,
                        "sample_date": target_date.isoformat(),
                        "granule_used": Path(candidate_path).name,
                        "n_soundings": n_soundings,
                        "xco2_ppm_mean": mean_xco2,
                        "xco2_ppm_stderr": stderr_xco2,
                    })
                    found_coverage = True
                    break

            if not found_coverage:
                print(f"  {target_date}: checked "
                      f"{len(candidate_paths)} candidate granules, none had "
                      f"soundings actually within the regional box "
                      f"(expected -- OCO's narrow swath doesn't cross "
                      f"every region every cycle)")

    if not rows:
        print("\nNo usable soundings found across any sample date. "
              "Nothing written.")
        return

    with open(OUTPUT_CSV, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} regional XCO2 samples to {OUTPUT_CSV}")
    print("\nHow to read this output:")
    print("  - This is a REGIONAL (~110 km box) mean, not a KFP-specific")
    print("    signal. Do not interpret any dip or rise as attributable")
    print("    to KFP's sequestration -- the project is far too small")
    print("    relative to OCO's footprint and noise floor for that.")
    print("  - A gently rising xco2_ppm_mean across years, roughly")
    print("    matching the global ~2-3 ppm/yr background trend, is the")
    print("    expected and unremarkable result. That is what Stage 4")
    print("    should describe this data as showing (or not showing).")
    print("  - Low n_soundings dates reflect real data gaps (cloud")
    print("    cover, orbit geometry) -- expected, not an error.")


if __name__ == "__main__":
    run_stage3()