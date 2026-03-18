//VERSION=3
/*
 * WATER THRESHOLDING (Simple Water Bodies Mapping)
 * =================================================
 * Sentinel-2 evalscript for detecting and mapping water surfaces using
 * multiple spectral water indices with configurable thresholds.
 *
 * Indices used:
 *   - MNDWI (Modified NDWI):  (B03 - B11) / (B03 + B11)
 *   - NDWI  (Normalized):     (B03 - B08) / (B03 + B08)
 *   - SWI   (Sentinel Water): (B05 - B11) / (B05 + B11)  [S2 only]
 *
 * A pixel is classified as water if ANY of the indices exceeds its threshold.
 *
 * Based on SWBM by Mohor Gartner.
 * Supports: Sentinel-2 L2A, Sentinel-2 L1C, Landsat 8-9
 *
 * Author: WildSquare (adapted from SWBM)
 * License: CC BY-SA 4.0
 *
 * Usage:
 *   - Data source: Sentinel-2 L2A (recommended), L1C, or Landsat 8-9
 *   - Recommended: <15% cloud coverage
 *   - Adjust thresholds below to calibrate for your region
 */

// --- Configuration -----------------------------------------------------------

// Source satellite: "S2L2A" | "S2L1C" | "L8"
var source = "S2L2A";

// Index to export for histogram/statistics (0=MNDWI, 1=NDWI, 2=SWI)
var histogramIndex = 0;

// Water detection thresholds (calibrate per region)
var MNDWI_thr = 0.1;   // Modified NDWI threshold
var NDWI_thr = 0.2;    // NDWI threshold
var SWI_thr = 0.03;    // Sentinel Water Index threshold (S2 only)

// Water overlay color [R, G, B] (0-1 range)
var waterColor = [0.0, 0.2, 0.8]; // Deep blue

// --- Input bands based on source ---------------------------------------------

var inputList;
switch (source) {
  case "S2L2A":
    inputList = ["B02", "B03", "B04", "B05", "B08", "B11", "SCL", "dataMask"];
    break;
  case "S2L1C":
    inputList = ["B02", "B03", "B04", "B05", "B08", "B11", "dataMask"];
    break;
  case "L8":
    inputList = ["B02", "B03", "B04", "B05", "B06", "dataMask"];
    break;
}

// --- Setup -------------------------------------------------------------------

function setup() {
  return {
    input: inputList,
    output: [
      { id: "default", bands: 4 },
      { id: "index", bands: 1, sampleType: "FLOAT32" },
      { id: "eobrowserStats", bands: 2, sampleType: "FLOAT32" },
      { id: "dataMask", bands: 1 },
    ],
  };
}

// --- Helper functions --------------------------------------------------------

function index(a, b) {
  return (a - b) / (a + b);
}

/**
 * Water body detection for Sentinel-2.
 * Returns 1 if water, 0 if not water.
 */
function detectWaterS2(green, nir, swir, vre1) {
  try {
    var mndwi = index(green, swir);
    var ndwi = index(green, nir);
    var swi = index(vre1, swir);
    return mndwi > MNDWI_thr || ndwi > NDWI_thr || swi > SWI_thr ? 1 : 0;
  } catch (err) {
    return 0;
  }
}

/**
 * Water body detection for Landsat 8-9.
 * Returns 1 if water, 0 if not water.
 */
function detectWaterL8(green, nir, swir) {
  try {
    var mndwi = index(green, swir);
    var ndwi = index(green, nir);
    return mndwi > MNDWI_thr || ndwi > NDWI_thr ? 1 : 0;
  } catch (err) {
    return 0;
  }
}

/**
 * Cloud detection using SCL (Scene Classification Layer).
 */
function isCloud(scl) {
  return scl === 8 || scl === 9;
}

/**
 * Calculate the selected index value for statistics/histogram.
 */
function calcIndexValue(histIdx, green, nir, swir, vre1) {
  if (histIdx === 0) return index(green, swir);     // MNDWI
  if (histIdx === 1) return index(green, nir);       // NDWI
  if (histIdx === 2) return index(vre1, swir);       // SWI
  return NaN;
}

// --- Main evaluation ---------------------------------------------------------

function evaluatePixel(p) {
  var b = p.B02, g = p.B03, r = p.B04;
  var dataMask = p.dataMask;

  var nir, swir, vre1;
  vre1 = NaN;

  if (source === "S2L2A" || source === "S2L1C") {
    nir = p.B08;
    swir = p.B11;
    vre1 = p.B05;
  } else {
    // Landsat 8-9
    nir = p.B05;
    swir = p.B06;
  }

  // Water detection
  var w;
  if (source === "S2L2A" || source === "S2L1C") {
    w = detectWaterS2(g, nir, swir, vre1);
  } else {
    w = detectWaterL8(g, nir, swir);
  }

  // Index value for statistics
  var indexVal = dataMask === 1
    ? calcIndexValue(histogramIndex, g, nir, swir, vre1)
    : NaN;

  // Cloud check (S2L2A only)
  var cloud = false;
  if (source === "S2L2A") {
    cloud = isCloud(p.SCL);
  }

  // True color RGB
  var rgb = [r * 2.5, g * 2.5, b * 2.5, dataMask];

  // Output pixel
  var outPixel;
  if (cloud) {
    outPixel = rgb; // Show true color for cloudy pixels
  } else if (w === 1) {
    outPixel = [waterColor[0], waterColor[1], waterColor[2], dataMask];
  } else {
    outPixel = rgb;
  }

  return {
    default: outPixel,
    index: [indexVal],
    eobrowserStats: [w, cloud ? 1 : 0],
    dataMask: [dataMask],
  };
}
