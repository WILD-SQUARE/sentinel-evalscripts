"""
Water Quality Analysis — Trebujena, Cadiz (Guadalquivir)
========================================================
Standalone script that uses Google Earth Engine + Sentinel-2 L2A
to compute 3 water quality indicators for a 6-hectare territory:

  1. Water Turbidity  (NTU)   — Se2WaQ empirical model
  2. Water Thresholding       — MNDWI / NDWI / SWI multi-index
  3. Water Chlorophyll (NDCI) — Chl-a estimation

Usage:
  python water_quality_analysis.py

Requirements:
  pip install earthengine-api google-auth
  ee.Authenticate()  # first time only

Author: WildSquare
"""

import ee
import json
import sys
from datetime import datetime

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

# Google Earth Engine project
GEE_PROJECT = "n8n-workflow-482712"

# Territory: Trebujena, Cadiz — Guadalquivir marshlands (~6 ha)
# Approximate polygon covering wetland/river area near Trebujena
TERRITORY_NAME = "Guadalquivir - Trebujena, Cadiz"
TERRITORY_COORDS = [
    [-6.1780, 36.8680],
    [-6.1750, 36.8680],
    [-6.1750, 36.8650],
    [-6.1720, 36.8650],
    [-6.1720, 36.8620],
    [-6.1780, 36.8620],
    [-6.1780, 36.8680],
]

# Analysis period
DATE_START = "2024-06-01"
DATE_END   = "2025-03-01"

# Sentinel-2 config
CLOUD_MAX_PCT  = 15       # Max cloud cover percentage
SCALE          = 10       # Resolution in meters
MAX_PIXELS     = 1e7

# Water detection thresholds
MNDWI_THR = 0.1
NDWI_THR  = 0.2
SWI_THR   = 0.03

# Thumbnail size
THUMB_DIM = "512x512"

# Color palettes
TURBIDITY_PALETTE = ['4973F2', '82D35F', 'FEFD05', 'FD0004', '8E2026', 'D97CF5']
NDCI_PALETTE      = ['313695', '4575b4', 'e0f3f8', 'fee090', 'fdae61', 'f46d43', 'a50026']
WATER_PALETTE     = ['FFFFFF', '0033CC']
TRUECOLOR_BANDS   = ['B4', 'B3', 'B2']


# ─────────────────────────────────────────────────────────────────────────────
# INITIALIZE GEE
# ─────────────────────────────────────────────────────────────────────────────

def init_gee():
    """Initialize Google Earth Engine."""
    try:
        ee.Initialize(project=GEE_PROJECT)
        print(f"[OK] Google Earth Engine initialized (project: {GEE_PROJECT})")
    except Exception:
        print("[!] Authenticating with Google Earth Engine...")
        ee.Authenticate()
        ee.Initialize(project=GEE_PROJECT)
        print(f"[OK] Authenticated and initialized (project: {GEE_PROJECT})")


# ─────────────────────────────────────────────────────────────────────────────
# SPECTRAL INDICES
# ─────────────────────────────────────────────────────────────────────────────

def add_water_indices(image):
    """
    Add all water quality bands to a Sentinel-2 image:
      - NDWI:      (B03 - B08) / (B03 + B08)
      - MNDWI:     (B03 - B11) / (B03 + B11)
      - SWI:       (B05 - B11) / (B05 + B11)
      - NDCI:      (B05 - B04) / (B05 + B04)
      - TURBIDITY: 8.93 * (B03 / B01) - 6.39
      - CHL_A:     4.26 * (B03 / B01)^3.94
    """
    ndwi  = image.normalizedDifference(['B3', 'B8']).rename('NDWI')
    mndwi = image.normalizedDifference(['B3', 'B11']).rename('MNDWI')
    swi   = image.normalizedDifference(['B5', 'B11']).rename('SWI')
    ndci  = image.normalizedDifference(['B5', 'B4']).rename('NDCI')

    ratio_b3_b1 = image.select('B3').divide(image.select('B1'))
    turbidity = ratio_b3_b1.multiply(8.93).subtract(6.39).rename('TURBIDITY')

    chla = ratio_b3_b1.pow(3.94).multiply(4.26).rename('CHL_A')

    return image.addBands([ndwi, mndwi, swi, ndci, turbidity, chla])


def compute_water_mask(image):
    """
    Classify water pixels using multi-index thresholding:
      water = (MNDWI > 0.1) OR (NDWI > 0.2) OR (SWI > 0.03)
    Also uses SCL class 6 (water) as fallback.
    Returns a binary mask (1 = water, 0 = land).
    """
    mndwi_mask = image.select('MNDWI').gt(MNDWI_THR)
    ndwi_mask  = image.select('NDWI').gt(NDWI_THR)
    swi_mask   = image.select('SWI').gt(SWI_THR)

    # SCL water class
    scl_water = image.select('SCL').eq(6)

    # Combined: any index OR SCL
    water_mask = mndwi_mask.Or(ndwi_mask).Or(swi_mask).Or(scl_water).rename('WATER')
    return water_mask


def compute_cloud_mask(image):
    """
    Cloud mask from SCL: classes 8 (cloud medium) and 9 (cloud high).
    Returns binary mask (1 = cloud).
    """
    scl = image.select('SCL')
    return scl.eq(8).Or(scl.eq(9)).rename('CLOUD')


# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def analyze_territory(geom, date_start, date_end):
    """
    Run the full water quality analysis on a territory.
    Returns dict with statistics and thumbnail URLs.
    """
    print(f"\n{'='*60}")
    print(f"  TERRITORY: {TERRITORY_NAME}")
    print(f"  Period:    {date_start} -> {date_end}")
    print(f"{'='*60}\n")

    # ── Get best Sentinel-2 image ──
    collection = (
        ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(geom)
        .filterDate(date_start, date_end)
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_MAX_PCT))
        .sort('CLOUDY_PIXEL_PERCENTAGE')
    )

    count = collection.size().getInfo()
    print(f"[INFO] Found {count} Sentinel-2 images with <{CLOUD_MAX_PCT}% clouds")

    if count == 0:
        print("[ERROR] No images found. Try increasing CLOUD_MAX_PCT or date range.")
        return None

    # Use best (least cloudy) image
    image = collection.first()
    image_date = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd').getInfo()
    image_id = image.get('system:index').getInfo()
    cloud_pct = image.get('CLOUDY_PIXEL_PERCENTAGE').getInfo()
    print(f"[INFO] Selected image: {image_id}")
    print(f"[INFO] Date: {image_date} | Cloud cover: {cloud_pct:.1f}%")

    # ── Compute indices ──
    image_idx = add_water_indices(image)
    water_mask = compute_water_mask(image_idx)
    cloud_mask = compute_cloud_mask(image_idx)
    image_all = image_idx.addBands([water_mask, cloud_mask])

    # ── Statistics over the full territory ──
    print("\n[PROCESSING] Computing territory-wide statistics...")
    stats = image_all.select([
        'NDWI', 'MNDWI', 'SWI', 'NDCI', 'TURBIDITY', 'CHL_A', 'WATER', 'CLOUD'
    ]).reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=geom,
        scale=SCALE,
        maxPixels=MAX_PIXELS,
        bestEffort=True,
    ).getInfo()

    # ── Statistics over WATER pixels only ──
    print("[PROCESSING] Computing water-only statistics...")
    water_only = image_all.updateMask(water_mask)
    water_stats = water_only.select([
        'NDCI', 'TURBIDITY', 'CHL_A'
    ]).reduceRegion(
        reducer=ee.Reducer.mean().combine(
            ee.Reducer.minMax(), sharedInputs=True
        ),
        geometry=geom,
        scale=SCALE,
        maxPixels=MAX_PIXELS,
        bestEffort=True,
    ).getInfo()

    # ── Water area calculation ──
    print("[PROCESSING] Computing water surface area...")
    water_area_img = water_mask.multiply(ee.Image.pixelArea())
    water_area = water_area_img.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=geom,
        scale=SCALE,
        maxPixels=MAX_PIXELS,
        bestEffort=True,
    ).get('WATER').getInfo()

    total_area_img = ee.Image.pixelArea()
    total_area = total_area_img.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=geom,
        scale=SCALE,
        maxPixels=MAX_PIXELS,
        bestEffort=True,
    ).get('area').getInfo()

    water_area_ha = (water_area or 0) / 10000
    total_area_ha = (total_area or 0) / 10000
    water_pct = (water_area_ha / total_area_ha * 100) if total_area_ha > 0 else 0

    # ── Generate thumbnails ──
    print("[PROCESSING] Generating thumbnails...")
    geom_info = geom.getInfo()
    thumbnails = {}

    try:
        # True color
        thumbnails['true_color'] = image.select(TRUECOLOR_BANDS).getThumbUrl({
            'region': geom_info, 'dimensions': THUMB_DIM, 'format': 'png',
            'min': 0, 'max': 3000,
        })

        # Water mask
        thumbnails['water_mask'] = water_mask.getThumbUrl({
            'region': geom_info, 'dimensions': THUMB_DIM, 'format': 'png',
            'min': 0, 'max': 1, 'palette': WATER_PALETTE,
        })

        # Turbidity (water only)
        thumbnails['turbidity'] = water_only.select('TURBIDITY').getThumbUrl({
            'region': geom_info, 'dimensions': THUMB_DIM, 'format': 'png',
            'min': 0, 'max': 20, 'palette': TURBIDITY_PALETTE,
        })

        # NDCI (water only)
        thumbnails['ndci'] = water_only.select('NDCI').getThumbUrl({
            'region': geom_info, 'dimensions': THUMB_DIM, 'format': 'png',
            'min': -0.2, 'max': 0.4, 'palette': NDCI_PALETTE,
        })

        # Chl-a (water only)
        thumbnails['chlorophyll_a'] = water_only.select('CHL_A').getThumbUrl({
            'region': geom_info, 'dimensions': THUMB_DIM, 'format': 'png',
            'min': 0, 'max': 50, 'palette': TURBIDITY_PALETTE,
        })

        print("[OK] Thumbnails generated successfully")
    except Exception as e:
        print(f"[WARN] Thumbnail error: {e}")

    # ── Build results ──
    def safe_round(val, decimals=4):
        return round(val, decimals) if val is not None else None

    results = {
        "territory": TERRITORY_NAME,
        "analysis_date": datetime.now().isoformat(),
        "image": {
            "id": image_id,
            "date": image_date,
            "cloud_cover_pct": round(cloud_pct, 1),
            "satellite": "Sentinel-2 L2A (Harmonized)",
        },
        "area": {
            "total_ha": safe_round(total_area_ha, 2),
            "water_ha": safe_round(water_area_ha, 2),
            "water_pct": safe_round(water_pct, 1),
        },
        "indicators": {
            "water_thresholding": {
                "description": "Multi-index water body detection (MNDWI/NDWI/SWI)",
                "method": "pixel_water = (MNDWI > 0.1) OR (NDWI > 0.2) OR (SWI > 0.03)",
                "mean_MNDWI": safe_round(stats.get('MNDWI')),
                "mean_NDWI": safe_round(stats.get('NDWI')),
                "mean_SWI": safe_round(stats.get('SWI')),
                "water_fraction": safe_round(stats.get('WATER')),
                "cloud_fraction": safe_round(stats.get('CLOUD')),
                "thresholds": {
                    "MNDWI": MNDWI_THR,
                    "NDWI": NDWI_THR,
                    "SWI": SWI_THR,
                },
            },
            "water_turbidity": {
                "description": "Turbidity in NTU (Se2WaQ empirical model)",
                "formula": "Turb = 8.93 * (B03/B01) - 6.39",
                "unit": "NTU",
                "water_pixels_mean": safe_round(water_stats.get('TURBIDITY_mean')),
                "water_pixels_min": safe_round(water_stats.get('TURBIDITY_min')),
                "water_pixels_max": safe_round(water_stats.get('TURBIDITY_max')),
                "interpretation": classify_turbidity(water_stats.get('TURBIDITY_mean')),
            },
            "water_chlorophyll": {
                "description": "Chlorophyll-a via NDCI + empirical Chl_a model",
                "ndci": {
                    "formula": "NDCI = (B05 - B04) / (B05 + B04)",
                    "water_pixels_mean": safe_round(water_stats.get('NDCI_mean')),
                    "water_pixels_min": safe_round(water_stats.get('NDCI_min')),
                    "water_pixels_max": safe_round(water_stats.get('NDCI_max')),
                    "interpretation": classify_ndci(water_stats.get('NDCI_mean')),
                },
                "chlorophyll_a": {
                    "formula": "Chl_a = 4.26 * (B03/B01)^3.94",
                    "unit": "mg/m3",
                    "water_pixels_mean": safe_round(water_stats.get('CHL_A_mean')),
                    "water_pixels_min": safe_round(water_stats.get('CHL_A_min')),
                    "water_pixels_max": safe_round(water_stats.get('CHL_A_max')),
                    "trophic_state": classify_trophic(water_stats.get('CHL_A_mean')),
                },
            },
        },
        "thumbnails": thumbnails,
        "references": [
            "McFeeters (1996) - NDWI. Int. J. Remote Sensing, 17(7), 1425-1432",
            "Xu (2006) - MNDWI. Int. J. Remote Sensing, 27(14), 3025-3033",
            "Mishra & Mishra (2012) - NDCI. Remote Sens. Environ., 117, 394-406",
            "Potes et al. (2018) - Se2WaQ. Proc. IAHS, 380, 73-79",
            "Toming et al. (2016) - S2 water quality. Remote Sensing, 8(8), 640",
        ],
    }

    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLASSIFICATION HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def classify_turbidity(ntu):
    """Classify turbidity level from NTU value."""
    if ntu is None:
        return "No data"
    if ntu < 4:
        return "Agua limpia"
    elif ntu < 8:
        return "Turbidez baja"
    elif ntu < 12:
        return "Turbidez moderada"
    elif ntu < 16:
        return "Turbidez alta"
    elif ntu < 20:
        return "Turbidez muy alta"
    else:
        return "Turbidez extrema"


def classify_ndci(ndci):
    """Classify chlorophyll level from NDCI value."""
    if ndci is None:
        return "No data"
    if ndci < -0.2:
        return "Clorofila muy baja (agua clara)"
    elif ndci < 0.0:
        return "Clorofila baja"
    elif ndci < 0.2:
        return "Clorofila moderada"
    elif ndci < 0.4:
        return "Clorofila elevada"
    else:
        return "Clorofila alta / bloom algal"


def classify_trophic(chla):
    """Classify trophic state from Chl-a concentration (mg/m3)."""
    if chla is None:
        return "No data"
    if chla < 6:
        return "Oligotrofico"
    elif chla < 12:
        return "Mesotrofico"
    elif chla < 20:
        return "Meso-eutrofico"
    elif chla < 30:
        return "Eutrofico"
    elif chla < 50:
        return "Hipereutrofico"
    else:
        return "Bloom algal"


# ─────────────────────────────────────────────────────────────────────────────
# PRETTY PRINT
# ─────────────────────────────────────────────────────────────────────────────

def print_results(results):
    """Print results in a readable format to the console."""
    if not results:
        return

    ind = results['indicators']
    wt = ind['water_thresholding']
    tb = ind['water_turbidity']
    ch = ind['water_chlorophyll']

    print(f"\n{'='*60}")
    print(f"  RESULTS: {results['territory']}")
    print(f"  Image: {results['image']['date']} ({results['image']['cloud_cover_pct']}% clouds)")
    print(f"{'='*60}")

    print(f"\n  AREA")
    print(f"  {'Total:':<25} {results['area']['total_ha']} ha")
    print(f"  {'Water surface:':<25} {results['area']['water_ha']} ha ({results['area']['water_pct']}%)")

    print(f"\n  {'─'*56}")
    print(f"  1. WATER THRESHOLDING")
    print(f"  {'─'*56}")
    print(f"  {'Mean MNDWI:':<25} {wt['mean_MNDWI']}")
    print(f"  {'Mean NDWI:':<25} {wt['mean_NDWI']}")
    print(f"  {'Mean SWI:':<25} {wt['mean_SWI']}")
    print(f"  {'Water fraction:':<25} {wt['water_fraction']}")
    print(f"  {'Cloud fraction:':<25} {wt['cloud_fraction']}")

    print(f"\n  {'─'*56}")
    print(f"  2. WATER TURBIDITY")
    print(f"  {'─'*56}")
    print(f"  {'Mean (water):':<25} {tb['water_pixels_mean']} NTU")
    print(f"  {'Min:':<25} {tb['water_pixels_min']} NTU")
    print(f"  {'Max:':<25} {tb['water_pixels_max']} NTU")
    print(f"  {'Classification:':<25} {tb['interpretation']}")

    print(f"\n  {'─'*56}")
    print(f"  3. WATER CHLOROPHYLL")
    print(f"  {'─'*56}")
    print(f"  {'Mean NDCI (water):':<25} {ch['ndci']['water_pixels_mean']}")
    print(f"  {'NDCI range:':<25} [{ch['ndci']['water_pixels_min']}, {ch['ndci']['water_pixels_max']}]")
    print(f"  {'NDCI interpretation:':<25} {ch['ndci']['interpretation']}")
    print(f"  {'Mean Chl-a (water):':<25} {ch['chlorophyll_a']['water_pixels_mean']} mg/m3")
    print(f"  {'Chl-a range:':<25} [{ch['chlorophyll_a']['water_pixels_min']}, {ch['chlorophyll_a']['water_pixels_max']}] mg/m3")
    print(f"  {'Trophic state:':<25} {ch['chlorophyll_a']['trophic_state']}")

    print(f"\n  {'─'*56}")
    print(f"  THUMBNAILS")
    print(f"  {'─'*56}")
    for name, url in results.get('thumbnails', {}).items():
        print(f"  {name + ':':<25} {url[:80]}..." if url else f"  {name + ':':<25} N/A")

    print(f"\n{'='*60}\n")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    init_gee()

    # Build geometry
    geom = ee.Geometry.Polygon([TERRITORY_COORDS])

    # Run analysis
    results = analyze_territory(geom, DATE_START, DATE_END)

    if results:
        # Print to console
        print_results(results)

        # Save JSON
        output_path = "water_quality_results.json"
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        print(f"[OK] Results saved to {output_path}")
    else:
        print("[ERROR] Analysis failed. No results.")
        sys.exit(1)


if __name__ == "__main__":
    main()
