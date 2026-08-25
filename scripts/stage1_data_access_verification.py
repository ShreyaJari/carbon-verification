"""
stage1_data_access_verification.py

Satellite MRV Toolkit for Carbon Removal Verification
Stage 1: Case-study selection and data access verification

Case study: Kachung Forest Project (KFP)
    - CDM Project 4653, host party Uganda
    - Registered 4 Apr 2011, methodology AR-AM0004 v4
    - Location: Kachung Central Forest Reserve, Dokolo District, Uganda
    - CDM-eligible project area: 2,099 ha, bounded 1 58'56"N-2 02'32"N,
      32 54'55"E-32 59'43"E (from the CCBA/CDM Project Design Document,
      Section G1.3)
    - Species: ~90% Pinus caribaea hondurensis, remainder Eucalyptus spp.
      and Maesopsis eminii
    - Registry-claimed sequestration (for later comparison stage):
        issued CERs 2006-2012 monitoring period: 30,492 tCO2e
        issued CERs 2013-2020 monitoring period: 314,672 tCO2e

Purpose of this script
-----------------------
Before building any inversion or comparison logic, confirm that the three
satellite data sources this project depends on actually have usable
coverage over the KFP project area and monitoring period:
    1. GEDI L4A/L4B (aboveground biomass, via Google Earth Engine)
    2. Sentinel-2 surface reflectance (canopy cover / NDVI, via GEE)
    3. OCO-2/OCO-3 XCO2 column retrievals (via NASA Earthdata / earthaccess)

This script only counts/reports data availability. It does not compute
any biomass or carbon estimates yet -- that is Stage 2.

Requirements (install locally, not run in this environment)
--------------------------------------------------------------
    pip install earthengine-api earthaccess

Before running:
    1. Authenticate Earth Engine once: `earthengine authenticate`
       (or run ee.Authenticate() interactively the first time)
    2. Have a free NASA Earthdata Login account for the earthaccess checks
       (https://urs.earthdata.nasa.gov/) -- earthaccess.login() will
       prompt for credentials or read them from a .netrc file.

Run this script locally in VS Code and paste the full terminal output
back so we can decide, from real results, whether to proceed to Stage 2
or adjust the case study / date range.
"""

import sys
from datetime import date

# ---------------------------------------------------------------------------
# 1. Project area definition
# ---------------------------------------------------------------------------
# Bounding box taken directly from the CCBA/CDM PDD, Section G1.3:
# "The project boundary area of land is 2,099 ha confined within 3,500 ha
#  of Reserve land, located between 1 58'56" N to 2 02'32" N and
#  32 54'55" E to 32 59'43" E."
KFP_LAT_MIN_DEG = 1.98222   # 1 deg 58' 56" N
KFP_LAT_MAX_DEG = 2.04222   # 2 deg 02' 32" N
KFP_LON_MIN_DEG = 32.91528  # 32 deg 54' 55" E
KFP_LON_MAX_DEG = 32.99528  # 32 deg 59' 43" E

# Monitoring periods per the registered CERs, used to scope the date
# ranges we check data availability for.
MONITORING_PERIOD_1_START = date(2006, 10, 1)
MONITORING_PERIOD_1_END = date(2012, 11, 22)
MONITORING_PERIOD_2_START = date(2013, 1, 1)
MONITORING_PERIOD_2_END = date(2020, 12, 31)

# GEDI's ISS-orbit inclination limits coverage to +/-51.6 deg latitude.
# KFP sits at ~2 deg N, so this is not a coverage constraint here --
# recorded explicitly so the assumption is visible, not assumed silently.
GEDI_MAX_ABS_LATITUDE_DEG = 51.6

# Google Cloud project registered for Earth Engine access (required as of
# EE's project-based auth model -- ee.Initialize() with no project arg
# fails even with valid credentials if none is specified).
EE_PROJECT_ID = "carbon-verification-toolkit"


def check_gedi_latitude_coverage():
    """Sanity check that KFP falls inside the GEDI ISS-orbit latitude band."""
    in_range = abs(KFP_LAT_MAX_DEG) <= GEDI_MAX_ABS_LATITUDE_DEG
    print(f"GEDI latitude coverage check: KFP max latitude "
          f"{KFP_LAT_MAX_DEG:.4f} deg N, GEDI limit "
          f"+/-{GEDI_MAX_ABS_LATITUDE_DEG} deg -> "
          f"{'OK' if in_range else 'OUT OF RANGE'}")
    return in_range


# ---------------------------------------------------------------------------
# 2. Google Earth Engine checks: GEDI L4A, GEDI L4B, Sentinel-2
# ---------------------------------------------------------------------------
def run_earth_engine_checks():
    """Query GEE for GEDI L4A/L4B footprint counts and Sentinel-2 scene
    counts intersecting the KFP bounding box, across both monitoring
    periods where the collection's temporal coverage allows."""
    try:
        import ee
    except ImportError:
        print("earthengine-api is not installed locally. "
              "Run: pip install earthengine-api")
        return

    try:
        ee.Initialize(project=EE_PROJECT_ID)
    except Exception as init_error:
        print(f"ee.Initialize(project='{EE_PROJECT_ID}') failed: "
              f"{init_error}")
        print("Run `earthengine authenticate` in your terminal first, "
              "then retry. Also confirm the Earth Engine API is enabled "
              "for this project in Google Cloud Console and that the "
              "project is registered at https://code.earthengine."
              "google.com/register.")
        return

    aoi = ee.Geometry.Rectangle(
        [KFP_LON_MIN_DEG, KFP_LAT_MIN_DEG, KFP_LON_MAX_DEG, KFP_LAT_MAX_DEG]
    )

    # --- GEDI L4A (footprint-level aboveground biomass density) ---
    # Mission operated Apr 2019 - present (with a funding-driven gap;
    # verify current end date in the printed collection metadata below).
    gedi_l4a = ee.ImageCollection("LARSE/GEDI/GEDI04_A_002_MONTHLY")
    gedi_l4a_aoi = gedi_l4a.filterBounds(aoi)
    gedi_l4a_count = gedi_l4a_aoi.size().getInfo()
    print(f"\nGEDI L4A monthly composite images intersecting KFP AOI: "
          f"{gedi_l4a_count}")
    if gedi_l4a_count > 0:
        first_date = ee.Date(
            gedi_l4a_aoi.first().get("system:time_start")
        ).format("YYYY-MM-dd").getInfo()
        print(f"  Example available month: {first_date}")

    # --- GEDI L4B (single global 1km gridded mean AGBD raster) ---
    # Unlike L4A, L4B is published as one static Image asset, not a
    # time-indexed collection -- check for valid (non-masked) pixels
    # over the AOI instead of counting collection members.
    gedi_l4b = ee.Image("LARSE/GEDI/GEDI04_B_002").select("MU")
    gedi_l4b_valid_pixels = (
        gedi_l4b.reduceRegion(
            reducer=ee.Reducer.count(),
            geometry=aoi,
            scale=1000,
            maxPixels=1e9,
        )
        .get("MU")
        .getInfo()
    )
    print(f"GEDI L4B valid (non-masked) 1km AGBD pixels over KFP AOI: "
          f"{gedi_l4b_valid_pixels}")

    # --- Sentinel-2 surface reflectance (harmonized) ---
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate("2015-06-23", "2026-08-24")  # mission start to today
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
    )
    s2_count = s2.size().getInfo()
    print(f"\nSentinel-2 SR scenes over KFP AOI (<40% cloud, "
          f"2015-06-23 to present): {s2_count}")

    # Scene counts specifically within the second monitoring period,
    # which is the one we would use for an independent cross-check
    # against the 314,672-CER claim.
    s2_period2 = s2.filterDate(
        MONITORING_PERIOD_2_START.isoformat(),
        MONITORING_PERIOD_2_END.isoformat(),
    )
    s2_period2_count = s2_period2.size().getInfo()
    print(f"Sentinel-2 SR scenes within monitoring period 2 "
          f"({MONITORING_PERIOD_2_START} to {MONITORING_PERIOD_2_END}): "
          f"{s2_period2_count}")


# ---------------------------------------------------------------------------
# 3. NASA Earthdata check: OCO-2 / OCO-3 XCO2 granule availability
# ---------------------------------------------------------------------------
def run_earthaccess_checks():
    """Search NASA Earthdata (via earthaccess) for OCO-2 and OCO-3 Lite
    files whose spatial extent overlaps a regional box around KFP.

    Note: OCO-2/3 footprints are ~1.29 km x 2.25 km with sparse orbital
    tracks, so we search a generously padded regional box (roughly
    +/-1 degree, ~110 km) around KFP rather than the tight project
    boundary -- the toolkit's atmospheric stage was always scoped as a
    coarse regional cross-check, not project-scale detection, and this
    search radius reflects that honestly."""
    try:
        import earthaccess
    except ImportError:
        print("earthaccess is not installed locally. "
              "Run: pip install earthaccess")
        return

    try:
        earthaccess.login()
    except Exception as login_error:
        print(f"earthaccess.login() failed: {login_error}")
        print("Set up a free NASA Earthdata Login at "
              "https://urs.earthdata.nasa.gov/ and retry.")
        return

    regional_lon_min = KFP_LON_MIN_DEG - 1.0
    regional_lon_max = KFP_LON_MAX_DEG + 1.0
    regional_lat_min = KFP_LAT_MIN_DEG - 1.0
    regional_lat_max = KFP_LAT_MAX_DEG + 1.0

    for short_name, label in [
        ("OCO2_L2_Lite_FP", "OCO-2 Lite FP (v11)"),
        ("OCO3_L2_Lite_FP", "OCO-3 Lite FP (v11)"),
    ]:
        try:
            results = earthaccess.search_data(
                short_name=short_name,
                bounding_box=(
                    regional_lon_min,
                    regional_lat_min,
                    regional_lon_max,
                    regional_lat_max,
                ),
                temporal=(
                    MONITORING_PERIOD_1_START.isoformat(),
                    date.today().isoformat(),
                ),
                count=20,
            )
            print(f"\n{label}: {len(results)} granule(s) found "
                  f"(search capped at 20) over regional box "
                  f"[{regional_lon_min:.2f}, {regional_lat_min:.2f}, "
                  f"{regional_lon_max:.2f}, {regional_lat_max:.2f}]")
        except Exception as search_error:
            print(f"\n{label}: search failed -- {search_error}")


# ---------------------------------------------------------------------------
# 4. Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 70)
    print("KFP Stage 1 data access verification")
    print("=" * 70)

    check_gedi_latitude_coverage()

    print("\n" + "-" * 70)
    print("Google Earth Engine checks (GEDI L4A/L4B, Sentinel-2)")
    print("-" * 70)
    run_earth_engine_checks()

    print("\n" + "-" * 70)
    print("NASA Earthdata checks (OCO-2/OCO-3)")
    print("-" * 70)
    run_earthaccess_checks()

    print("\nDone. Paste this full output back so we can confirm data "
          "access before starting Stage 2 (biomass stage).")
    sys.exit(0)