# Polynomial Surface Metrology with Curvature Adaptation

A C++ toolkit for analysing 3D laser-scan point clouds (ASCII `.ply`) of curved parts. It fits a polynomial reference surface to the part, separates the nominal surface from features cut into it (grooves, holes), aligns the cloud to the part's axes and measures feature areas. A small Python generator creates synthetic test clouds with known geometry.

Visual Studio project folder: `PolinomilPlaneWithCurvatireDapation`.

---

## What it does

1. **Loads** an ASCII PLY cloud and groups points into scan lines, one per unique Y value, each sorted by X.
2. **Crops** the cloud in X to a region of interest.
3. **Fits a polynomial surface** `z = f(x, y)` of degree 2, 3 or 4 using least squares on normalised coordinates.
4. **Picks the polynomial degree automatically** by comparing the average Gaussian curvature of fits at different degrees.
5. **Classifies points** by distance to the fitted surface:
   - 🟢 **green**: on the nominal surface (within the tolerance band)
   - 🔴 **red**: off the surface (grooves, holes, defects)
6. **Aligns** the cloud to the part's X or Y axis (several strategies available).
7. **Measures** groove and hole areas by integrating widths across scan lines.
8. **Annotates** the cloud with text labels (for example, measured values) rendered as point rectangles.
9. **Saves** a coloured PLY on a background thread.

---

## Repository layout

| File | Purpose |
|---|---|
| `main.cpp` | Entry point: configuration, load, run pipeline, save |
| `PlyUtilities.h/.cpp` | PLY reading/writing, scan-line storage, X cropping, asynchronous save queue, output folder creation |
| `MyMetrologyEngine.h/.cpp` | Base engine: plane-based colouring, X-axis alignment, end-region colouring, groove and hole area measurement, text rendering, bottom-point removal |
| `HyperME.h/.cpp` | Derived engine: polynomial surface fitting, curvature-adaptive degree selection, flattening, several alignment methods, scan displacement correction |
| `shapeGenrator.py` | `ShapeGenrator`: builds a gridded surface from a function, adds grooves and noise |
| `plyGenrator.py` | `PlyGenartor`: extends the generator with PLY export and 3D plotting |
| `generate.py` (example script) | Builds a test surface with two grooves and writes `output.ply` |

---

## Requirements

**C++**
- Windows (folder creation uses PowerShell; paths use `\`)
- Visual Studio 2019/2022, x64, C++14 or later
- No external libraries

**Python** (for synthetic data)
- Python 3.8+
- `numpy`, `matplotlib`

```bash
pip install numpy matplotlib
```

---

## Quick start

### 1. Generate a test cloud

```python
from plyGenrator import PlyGenartor

shape = lambda x, y: 0.000005*x**4 - 0.002*x**2 + 0.00005*y**4 - 0.002*y**2
s = PlyGenartor(shape)            # default: x in [-10, 10], y in [-20, 20], 0.04 mm step
s.interpolateShape()
s.addGroove(center_x=-4.5, width=2.0, depth=2.0, wall_width=1.0, center_y=0, length=10.0)
s.addGroove(center_x= 4.5, width=2.0, depth=2.0, wall_width=1.0, center_y=0, length=10.0)
s.addNoise(std=0.04)              # std in mm (0.04 mm = 40 µm)
s.save_PLY(r'...\x64\Debug\ply\output.ply')
```

### 2. Configure and run the C++ pipeline

Edit the configuration block in `main.cpp`:

```cpp
string plyPath      = "./ply/output.ply";  // input cloud
string sampleName   = "test_sample";       // output subfolder name
string outputFolder = "./output";          // output root

float sectionY = 50.0f;   // Y range passed to the engine
float sectionX = 20.0f;   // X crop width (keeps |x| <= sectionX / 2)

PlyUtilities plyUtils(outputFolder, sampleName, 40);  // 40 -> 0.04 mm between scan lines
```

Build (x64) and run from the folder that contains `ply/`, for example `x64/Debug`.

### 3. Output

`./output/test_sample/colored_0.ply`, an ASCII PLY with `x y z red green blue` per vertex. Open it in CloudCompare or MeshLab to inspect the classification.

---

## Input format

ASCII PLY with `x y z` float vertices:

```
ply
format ascii 1.0
element vertex N
property float x
property float y
property float z
end_header
-10.0000 -20.0000 0.1234
...
```

Points with the same Y (within 0.0001) form one scan line. Units are assumed to be millimetres.

The third `PlyUtilities` constructor argument (`devisor`) sets the scan pitch: `scanPerMM = devisor / 1000`, so `40` means 0.04 mm between lines. This value is used as the Y step when integrating areas.

---

## Algorithm details

### Surface fitting (`HyperME::fitHyperSurface`)

Coordinates are centred and scaled to roughly [-1, 1], then a least-squares fit is solved with Gaussian elimination and partial pivoting over up to 15 basis terms:

| Degree | Terms |
|---|---|
| 4 | x⁴, y⁴, x³y, xy³, x²y² |
| 3 | x³, y³, x²y, xy² |
| 2 | x², y², xy |
| 1 | x, y |
| 0 | 1 |

Terms above the requested degree are masked out.

### Excluding features from the fit (`HyperME::getSubCloud`)

Points that lie deeper than a threshold below the current surface are dropped before refitting. The threshold varies linearly along Y between `zMin` and `zMax`, so grooves and holes don't pull the reference surface down.

### Curvature-adaptive degree selection

`colorByPlaneDistance()`:

1. Fits a degree-2 surface to all points.
2. For degree pairs (2, 3) and (3, 4), fits each degree and computes the mean |Gaussian curvature| over green points:

   K = (f_xx · f_yy − f_xy²) / (1 + f_x² + f_y²)²

3. Takes the ratio of the lower-degree curvature to the higher-degree curvature and compares it with α = 0.6:
   - ratio(2,3) ≤ 0.6 → degree 2
   - otherwise degree 3, upgraded to 4 if ratio(3,4) > 0.6
4. Runs two classification passes with the chosen degree:
   - pass 1: band ±0.2 mm
   - pass 2: band ±0.08 mm, and any point above the surface is also counted as green

`curvToBand()` (a log-linear map from curvature to band width, with a separate table for holes) is available for curvature-driven band selection.

### Alignment

| Method | Approach |
|---|---|
| `MyMetrologyEngine::alignXAxis` | Seed line from the lowest green points, refined 50× by regression on points below the line |
| `HyperME::alignXAxis` | Line through bottom edge, refined with shrinking residual windows (1.0 → 0.02 mm) |
| `HyperME::alignYAxis` | Regression on the rightmost green point per scan line (outer 15% of Y trimmed, top 25% by X kept) |
| `HyperME::alignByEigenVector` | Principal axis from the 2×2 covariance of green points |
| `HyperME::applyRotation` | Rotates by a given angle |

All methods rotate about Z, store the angle in `rotTeta` and set `dYScaling = cos(angle)` to correct the Y step in area integration. Points near the fitted reference line are marked pink (255, 105, 180).

### Area measurement

- **Grooves** (`getGrooveAreas`): finds groove centres on the scan line nearest the middle of the measurement zone, then walks scan lines through the zone and integrates groove width with the trapezoidal rule. Grooves that run off the edge of a scan line are discarded as incomplete.
- **Holes** (`getHoleAreas`): finds gaps between adjacent green points, locates the top and bottom edges of each gap, removes duplicates by centre and radius, and integrates width between those edges.

Both return `GrooveData { area, centerX, centerY }` and drop results below `minArea`.

### Other utilities

- `flattenPointCloud()`: subtracts the fitted surface from Z.
- `scanDisplacementCorrrection()` / `scanDisplacementCorrrectionHoles()`: estimate lateral scan drift from green-edge X positions.
- `colorEndRegions(blueRegion, keepRegion)`: colours red points outside the measurement zone blue.
- `writeTextToCloud(text, cx, cy)`: renders digits, `.`, `-` and spaces with a 5×7 bitmap font as a black-on-white point rectangle at the top Z of the cloud.
- `removeBottomPoints(percent)`: removes the lowest fraction of the Z range.
- `shiftInX(shift)`: translates the cloud in X.

---

## Colour legend

| Colour | RGB | Meaning |
|---|---|---|
| Green | (0, 255, 0) | Nominal surface |
| Red | (255, 0, 0) | Off-surface / feature |
| Blue | (0, 0, 255) | Red point outside the measurement zone |
| Pink | (255, 105, 180) | Alignment reference line |
| Black / white | | Text label |

---

## Known issues / TODO

- `PlyUtilities::SavePLY` is a stub: it sizes vectors with `0.4` and never decodes `g_Scans`, so it writes no points.
- `saveColoredPLY` and `SavePLY` join `path` and the filename without a separator (the save thread adds `\\`).
- `enqueueSave` pushes to `saveQueue` without locking `saveMutex`; that's fine in `main.cpp` because the thread is started right before, but it isn't thread-safe in general.
- `readFromPLY` groups scan lines with an O(N × lines) search; large clouds load slowly.
- `MyMetrologyEngine::computePlane()` is empty, so the base-class `colorByPlaneDistance` needs `planeN` / `planeD` set elsewhere.
- In `getGrooveAreas`, `noPointLeft` / `noPointRight` aren't reset between grooves, so one incomplete groove can mark later ones as corrupted.
- Raw pointer arrays (`coloredScan`, `lineCounts`) are never freed in the destructor, and `fitSurfceOfDegree` reallocates `coloredScan` on every call.
- The Python generator's comment says "3 micron noise" but `std=0.04` is 40 µm.
