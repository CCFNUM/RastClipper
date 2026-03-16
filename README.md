# RastClipper

A lightweight, header-minimal C++17 library for computing polygon overlap
fractions via scan-line rasterization and depth-buffer compositing.

Given a **subject** polygon and one or more **clip** polygons in 2-D,
RastClipper rasterizes them onto a shared bitmap and returns the fractional
area of the subject that each clip polygon covers, together with the
**pixel-weighted centroid** of each overlap region in the original input
coordinate system.  The approach is topology-agnostic: triangles,
quadrilaterals, pentagons, and arbitrary simple polygons can be freely
mixed on either side.

## Motivation

Exact geometric polygon clipping (e.g. Sutherland-Hodgman, Weiler-Atherton)
becomes fragile when many nearly-degenerate intersections arise, as is common
in non-conformal mesh interfaces for computational fluid dynamics.  Bitmap
rasterization sidesteps these robustness issues entirely: overlap is
determined by pixel counting, which is insensitive to vertex coincidences,
sliver polygons, and floating-point edge cases.

RastClipper mirrors the type names of the
[Clipper2](https://github.com/AngusJohnson/Clipper2) library (`PointD`,
`PathD`, `PathsD`), making it straightforward to swap between the two
backends via a compile-time flag.

## How it works

The rasterization pipeline has five stages:

```
 Subject polygon          Clip polygons
       |                       |
       v                       v
  [1] Viewport mapping — fit subject bbox into [-1,1] with 5% padding
       |                       |
       v                       v
  [2] Gross clip — trivial accept/reject against [-1,1] cube
       |                       |
       v                       v
  [3] Fine clip — Sutherland-Hodgman against [-1,1] (partial polygons)
       |                       |
       v                       v
  [4] Scan-line fill — edge-table rasterization with depth-buffer test
       |
       v
  [5] Pixel counting — 4-neighbor erosion, per-colour fractions & centroids
```

### 1. Viewport mapping

The subject polygon's bounding box defines a coordinate frame that maps
into the [-1, 1] normalised viewport.  A 5 % padding factor ensures that
the subject sits slightly inside the viewport boundary, so rounding never
pushes its vertices outside the rasterizable region.  Clip polygons are
mapped through the same transform and may extend beyond [-1, 1].

### 2 & 3. Polygon clipping

Before rasterization, each polygon is tested against the [-1, 1] cube:

- **Gross clip** classifies all vertices against the six bounding planes.
  If every vertex lies outside a single plane, the polygon is trivially
  rejected.  If all vertices are inside, it is trivially accepted.
- **Fine clip** applies the Sutherland-Hodgman algorithm for polygons that
  straddle the viewport boundary, producing a clipped polygon with new
  vertices on the boundary edges.

### 4. Scan-line fill

Normalised vertices are scaled to 1-based pixel coordinates on an
`N x N` bitmap (default `N = 128`).  The scan-line rasterizer builds
an **edge table** from the polygon's non-horizontal edges, recording
the x-intercept slope and depth gradient per edge.  Edges are shortened
at local minima/maxima to produce correct point-sampled fill (the
topmost and rightmost edges are omitted, avoiding double-counting at
shared vertices).

For each scan line, active-edge intersections are sorted by x and
processed in pairs.  Within each span, depth is linearly interpolated
and compared against the depth buffer:

| Mode   | Behaviour                              | Used for |
|--------|----------------------------------------|----------|
| ZFGEQU | Always write (unconditional)           | Subject  |
| ZFLESS | Write only if depth < existing depth   | Clips    |

The subject is drawn first at depth 1.0 (ZFGEQU), then each clip polygon
is drawn at depth 0.0 (ZFLESS).  A clip pixel is written only where the
subject has already been drawn (depth 1.0 > 0.0), naturally restricting
overlap to the subject's interior.

### 5. Pixel counting and centroid accumulation

Interior pixels are counted using **4-neighbor erosion**: a pixel is
counted only if all four of its direct neighbours (up, down, left, right)
are also non-background.  This strips one pixel of boundary rounding
artefact from every polygon edge, producing area fractions that converge
to the true geometric value as resolution increases.

Each clip polygon is assigned a unique colour.  The fraction for clip *i*
is the number of interior pixels coloured *i* divided by the total number
of interior pixels (subject + all clips).

In the same loop, the **pixel coordinates** (I, J) of each interior pixel
are accumulated per colour.  Dividing by the pixel count gives the average
pixel position, which is then converted back through the inverse viewport
mapping (pixel → normalised [-1, 1] → input coordinates) to produce the
overlap centroid in the caller's original 2-D coordinate system.

## API

```cpp
namespace RastClipper
{

struct PointD { double x, y; };
using  PathD  = std::vector<PointD>;
using  PathsD = std::vector<PathD>;

// Result of intersection fraction computation.
struct IntersectionResult
{
    std::vector<double> fractions; // area fractions, one per clip
    std::vector<PointD> centroids; // overlap centroids in input 2-D coords,
                                   // one per clip (zero if no overlap)
};

// Signed polygon area (shoelace formula).
double Area(const PathD& path);

// Compute overlap fractions and centroids of each clip polygon
// against the subject.
IntersectionResult IntersectionFractions(
    const PathD&  subject,
    const PathsD& clips,
    int resolution = 128);

}
```

### Parameters

| Parameter    | Default | Description                                       |
|--------------|---------|---------------------------------------------------|
| `subject`    | —       | The reference polygon (arbitrary vertex count).    |
| `clips`      | —       | One or more opposing polygons.                     |
| `resolution` | 128     | Bitmap side length. Higher values improve accuracy at the cost of O(N^2) memory and time. |

### Return value

An `IntersectionResult` containing:

- **`fractions`** — a `std::vector<double>` of size `clips.size()`.  Each
  entry is the fraction of the subject's rasterized interior area covered
  by the corresponding clip polygon.  Values are in [0, 1].  When clips do
  not overlap, they sum to at most 1.0.  When clips *do* overlap, earlier
  clips (lower index) take priority in the depth buffer; the later clip
  only gets pixels that the earlier one did not already claim.

- **`centroids`** — a `std::vector<PointD>` of size `clips.size()`.  Each
  entry is the pixel-weighted centroid of the overlap region between the
  subject and the corresponding clip polygon, expressed in the original
  input coordinate system.  If a clip has zero overlap (fraction = 0), its
  centroid is `{0, 0}`.

## Building

RastClipper is a single-source static library with no external
dependencies beyond the C++17 standard library.

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Running the tests

Tests are built automatically when RastClipper is the top-level project:

```bash
cd build
ctest          # or run the executable directly:
./test_rastclipper
```

### Integration as a submodule

When consumed as a subdirectory of another CMake project (e.g. via
`add_subdirectory`), the test target is not built.  Link against the
library with:

```cmake
add_subdirectory(external/RastClipper)
target_link_libraries(your_target PRIVATE RastClipper)
```

RastClipper also supports `find_package` after installation:

```cmake
find_package(RastClipper REQUIRED)
target_link_libraries(your_target PRIVATE RastClipper::RastClipper)
```

## Limitations

- **Simple polygons only.**  Self-intersecting polygons (e.g. a
  figure-eight) will produce incorrect area fractions.  Polygons with
  holes are not supported; use separate polygons instead.

- **Resolution-dependent accuracy.**  At the default resolution of 128,
  area fractions have an error of roughly O(1/N), i.e. ~1 %.  Increasing
  resolution to 256 or 512 improves accuracy at the cost of quadratic
  growth in memory and rasterization time.

- **Boundary erosion bias.**  The 4-neighbor erosion step discards one
  ring of pixels along every polygon edge.  For very small polygons
  (occupying only a few pixels), the erosion can discard a significant
  fraction of the polygon's area, biasing the result.  Increasing
  resolution mitigates this.

- **Clip priority on overlap.**  When two clip polygons overlap within
  the subject, only the one drawn first (lower index in the `clips`
  vector) claims the contested pixels.  The later clip's fraction is
  reduced accordingly.  This is by design (depth-buffer compositing),
  not a bug.

- **2-D only.**  All geometry is two-dimensional.  For 3-D mesh
  interfaces, the caller must project polygons into a local 2-D frame
  before invoking RastClipper.

## Licence

BSD-3-Clause.  See file headers for details.

Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and Arts.
