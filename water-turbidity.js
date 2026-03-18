//VERSION=3
/*
 * WATER TURBIDITY MAPPING
 * =======================
 * Sentinel-2 L2A evalscript for mapping water turbidity (NTU).
 *
 * Based on Se2WaQ (Sentinel-2 Water Quality) empirical model:
 *   Turb = 8.93 * (B03 / B01) - 6.39   [NTU]
 *
 * References:
 *   [1] M. Potes et al., "Use of Sentinel 2-MSI for water quality monitoring
 *       at Alqueva reservoir, Portugal," Proc. IAHS, vol. 380, pp. 73-79, 2018.
 *   [2] K. Toming et al., "First Experiences in Mapping Lake Water Quality
 *       Parameters with Sentinel-2 MSI Imagery," Remote Sens., vol. 8, 2016.
 *
 * Author: WildSquare (adapted from Se2WaQ by Nuno S. A. Pereira)
 * License: CC BY-NC-SA 4.0
 *
 * Usage in Sentinel Hub / EO Browser:
 *   - Data source: Sentinel-2 L2A
 *   - Recommended: low cloud coverage (<15%)
 *
 * Color scale (NTU):
 *   Blue   (0)  -> Clean water
 *   Green  (4)  -> Low turbidity
 *   Yellow (8)  -> Moderate turbidity
 *   Red    (12) -> High turbidity
 *   Dark   (16) -> Very high turbidity
 *   Purple (20) -> Extreme turbidity
 */

function setup() {
  return {
    input: ["B01", "B02", "B03", "B04", "B08", "B11", "SCL", "dataMask"],
    output: [
      { id: "default", bands: 4 },
      { id: "index", bands: 1, sampleType: "FLOAT32" },
      { id: "eobrowserStats", bands: 2, sampleType: "FLOAT32" },
      { id: "dataMask", bands: 1 },
    ],
  };
}

// --- Configuration -----------------------------------------------------------

// NDWI threshold to separate water from land
var NDWI_THRESHOLD = 0.0;

// Turbidity color scale (NTU values)
var scaleTurb = [0, 4, 8, 12, 16, 20];

// Corresponding colors: blue -> green -> yellow -> red -> dark red -> purple
var s = 255;
var colorScale = [
  [73 / s, 111 / s, 242 / s],   // Blue   - clean
  [130 / s, 211 / s, 95 / s],   // Green  - low
  [254 / s, 253 / s, 5 / s],    // Yellow - moderate
  [253 / s, 0 / s, 4 / s],      // Red    - high
  [142 / s, 32 / s, 38 / s],    // Dark   - very high
  [217 / s, 124 / s, 245 / s],  // Purple - extreme
];

// --- Helper functions --------------------------------------------------------

function isWater(scl) {
  // SCL class 6 = Water
  return scl === 6;
}

function isCloud(scl) {
  // SCL classes 8 = Cloud medium probability, 9 = Cloud high probability
  return scl === 8 || scl === 9;
}

function index(a, b) {
  return (a - b) / (a + b);
}

// --- Main evaluation ---------------------------------------------------------

function evaluatePixel(p) {
  var dataMask = p.dataMask;

  // NDWI for water detection
  var ndwi = index(p.B03, p.B08);

  // Turbidity empirical model (NTU)
  var turb = 8.93 * (p.B03 / p.B01) - 6.39;

  // True color fallback for non-water pixels
  var trueColor = [p.B04 * 2.5, p.B03 * 2.5, p.B02 * 2.5, dataMask];

  // Determine if pixel is water (using NDWI + SCL)
  var water = ndwi > NDWI_THRESHOLD || isWater(p.SCL);
  var cloud = isCloud(p.SCL);

  var indexVal = dataMask === 1 ? turb : NaN;

  // Output: show turbidity color for water, true color for land
  var outPixel;
  if (cloud) {
    outPixel = [0.7, 0.7, 0.7, dataMask]; // Gray for clouds
  } else if (water) {
    var col = colorBlend(turb, scaleTurb, colorScale);
    outPixel = [col[0], col[1], col[2], dataMask];
  } else {
    outPixel = trueColor;
  }

  return {
    default: outPixel,
    index: [indexVal],
    eobrowserStats: [water ? turb : NaN, cloud ? 1 : 0],
    dataMask: [dataMask],
  };
}
