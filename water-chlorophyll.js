//VERSION=3
/*
 * WATER CHLOROPHYLL (NDCI - Normalized Difference Chlorophyll Index)
 * ==================================================================
 * Sentinel-2 evalscript for estimating chlorophyll-a concentration in
 * water bodies using the NDCI spectral index.
 *
 * Formula:
 *   NDCI = (B05 - B04) / (B05 + B04)
 *
 * Where:
 *   B04 = Red band (665 nm)
 *   B05 = Red Edge 1 (705 nm)
 *
 * Interpretation:
 *   NDCI < -0.2  ->  Very low chlorophyll (deep blue, clear water)
 *   NDCI ~  0.0  ->  Low chlorophyll (light blue)
 *   NDCI ~  0.2  ->  Moderate chlorophyll (orange, productive water)
 *   NDCI >  0.4  ->  High chlorophyll / algal bloom (dark red)
 *
 * Additionally computes Chl_a empirical estimate (mg/m3) from Se2WaQ:
 *   Chl_a = 4.26 * (B03/B01)^3.94
 *
 * References:
 *   [1] S. Mishra, D.R. Mishra, "Normalized difference chlorophyll index:
 *       A novel model for remote estimation of chlorophyll-a concentration
 *       in turbid productive waters," Remote Sens. Environ., 2012.
 *   [2] M. Potes et al., "Use of Sentinel 2-MSI for water quality
 *       monitoring at Alqueva reservoir, Portugal," Proc. IAHS, 2018.
 *
 * Author: WildSquare (adapted from Sentinel Hub NDCI + Se2WaQ)
 * License: CC BY-SA 4.0
 *
 * Usage:
 *   - Data source: Sentinel-2 L2A
 *   - Best for: lakes, reservoirs, coastal waters
 *   - Recommended: <15% cloud coverage
 */

// --- Configuration -----------------------------------------------------------

// Display mode:
//   "ndci"  -> Color ramp based on NDCI index (-0.2 to 0.4)
//   "chla"  -> Color ramp based on Chl_a estimate (0 to 50 mg/m3)
var displayMode = "ndci";

// NDWI threshold to separate water from land
var NDWI_THRESHOLD = 0.0;

// NDCI color ramp: dark blue -> light blue -> orange -> dark red
var ndciRamp = [
  [-0.2, 0x313695], // Dark blue  - very low chlorophyll
  [-0.1, 0x4575b4], // Blue       - low chlorophyll
  [0.0, 0xe0f3f8],  // Light blue - minimal
  [0.1, 0xfee090],  // Light yellow
  [0.2, 0xfdae61],  // Orange     - moderate chlorophyll
  [0.3, 0xf46d43],  // Red-orange - elevated
  [0.4, 0xa50026],  // Dark red   - high chlorophyll / bloom
];

// Chl_a color scale (mg/m3)
var scaleChla = [0, 6, 12, 20, 30, 50];
var sc = 255;
var colorScaleChla = [
  [73 / sc, 111 / sc, 242 / sc],   // Blue   - oligotrophic
  [130 / sc, 211 / sc, 95 / sc],   // Green  - mesotrophic
  [254 / sc, 253 / sc, 5 / sc],    // Yellow - meso-eutrophic
  [253 / sc, 0 / sc, 4 / sc],      // Red    - eutrophic
  [142 / sc, 32 / sc, 38 / sc],    // Dark   - hyper-eutrophic
  [217 / sc, 124 / sc, 245 / sc],  // Purple - algal bloom
];

// --- Setup -------------------------------------------------------------------

function setup() {
  return {
    input: ["B01", "B02", "B03", "B04", "B05", "B08", "B11", "SCL", "dataMask"],
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

function isWater(scl) {
  return scl === 6;
}

function isCloud(scl) {
  return scl === 8 || scl === 9;
}

/**
 * Convert hex color to [r, g, b] array (0-1 range).
 */
function hexToRgb(hex) {
  return [
    ((hex >> 16) & 0xff) / 255,
    ((hex >> 8) & 0xff) / 255,
    (hex & 0xff) / 255,
  ];
}

/**
 * Interpolate color from a hex color ramp.
 * ramp: array of [value, hexColor] pairs, sorted by value.
 */
function colorFromRamp(value, ramp) {
  if (value <= ramp[0][0]) return hexToRgb(ramp[0][1]);
  if (value >= ramp[ramp.length - 1][0]) return hexToRgb(ramp[ramp.length - 1][1]);

  for (var i = 0; i < ramp.length - 1; i++) {
    if (value >= ramp[i][0] && value < ramp[i + 1][0]) {
      var t = (value - ramp[i][0]) / (ramp[i + 1][0] - ramp[i][0]);
      var c0 = hexToRgb(ramp[i][1]);
      var c1 = hexToRgb(ramp[i + 1][1]);
      return [
        c0[0] + t * (c1[0] - c0[0]),
        c0[1] + t * (c1[1] - c0[1]),
        c0[2] + t * (c1[2] - c0[2]),
      ];
    }
  }
  return hexToRgb(ramp[ramp.length - 1][1]);
}

// --- Main evaluation ---------------------------------------------------------

function evaluatePixel(p) {
  var dataMask = p.dataMask;

  // Water detection: NDWI + SCL
  var ndwi = index(p.B03, p.B08);
  var water = ndwi > NDWI_THRESHOLD || isWater(p.SCL);
  var cloud = isCloud(p.SCL);

  // NDCI: Normalized Difference Chlorophyll Index
  var ndci = index(p.B05, p.B04);

  // Chl_a empirical model (mg/m3) from Se2WaQ
  var chla = 4.26 * Math.pow(p.B03 / p.B01, 3.94);

  // Index value for statistics
  var indexVal = dataMask === 1 ? ndci : NaN;

  // True color fallback
  var trueColor = [p.B04 * 2.5, p.B03 * 2.5, p.B02 * 2.5, dataMask];

  // Output pixel
  var outPixel;
  if (cloud) {
    outPixel = [0.7, 0.7, 0.7, dataMask]; // Gray for clouds
  } else if (water) {
    var col;
    if (displayMode === "chla") {
      col = colorBlend(chla, scaleChla, colorScaleChla);
    } else {
      col = colorFromRamp(ndci, ndciRamp);
    }
    outPixel = [col[0], col[1], col[2], dataMask];
  } else {
    outPixel = trueColor;
  }

  return {
    default: outPixel,
    index: [indexVal],
    eobrowserStats: [water ? ndci : NaN, cloud ? 1 : 0],
    dataMask: [dataMask],
  };
}
