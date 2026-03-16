// File       : test_rastclipper.cpp
// Created    : Sun Mar 16 2026
// Author     : Mhamad Mahdi Alloush
// Description: Standalone tests for RastClipper library.
//              Returns 0 on success, non-zero on failure.
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "RastClipper/RastClipper.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

static int testCount = 0;
static int failCount = 0;

static void check(bool cond, const std::string& msg)
{
    ++testCount;
    if (!cond)
    {
        std::cerr << "  FAIL: " << msg << "\n";
        ++failCount;
    }
    else
    {
        std::cout << "  PASS: " << msg << "\n";
    }
}

static void checkClose(double actual,
                        double expected,
                        double tol,
                        const std::string& msg)
{
    ++testCount;
    const double err = std::abs(actual - expected);
    if (err > tol)
    {
        std::cerr << "  FAIL: " << msg << " (actual=" << actual
                  << ", expected=" << expected << ", err=" << err
                  << ", tol=" << tol << ")\n";
        ++failCount;
    }
    else
    {
        std::cout << "  PASS: " << msg << " (actual=" << actual
                  << ", expected=" << expected << ", err=" << err << ")\n";
    }
}

// -----------------------------------------------------------------------
// Test 1: Identical squares — full overlap
//
// Subject and clip are both the unit square [(0,0),(1,0),(1,1),(0,1)].
// The area fraction should be 1.0 (entire subject is covered).
// -----------------------------------------------------------------------
static void testIdenticalSquares()
{
    std::cout << "\n--- Test 1: Identical squares (full overlap) ---\n";

    using namespace RastClipper;

    PathD unitSquare = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};

    // Verify Area()
    const double area = Area(unitSquare);
    checkClose(area, 1.0, 1e-12, "Area of unit square = 1.0");

    // Full overlap: subject == clip
    PathsD clips = {unitSquare};
    const int NBIT = 256;
    auto result = IntersectionFractions(unitSquare, clips, NBIT);

    check(result.fractions.size() == 1, "Returned 1 fraction");
    check(result.centroids.size() == 1, "Returned 1 centroid");

    // Rasterization tolerance: the 4-neighbor erosion shaves ~1 pixel ring
    // from both subject and clip identically, so the ratio should still be
    // very close to 1.0.  At NBIT=256 the tolerance is ~2%.
    checkClose(result.fractions[0], 1.0, 0.02,
               "Fraction ~ 1.0 for identical squares");

    // Centroid of full overlap with unit square should be at (0.5, 0.5)
    checkClose(result.centroids[0].x, 0.5, 0.05,
               "Centroid x ~ 0.5 for identical squares");
    checkClose(result.centroids[0].y, 0.5, 0.05,
               "Centroid y ~ 0.5 for identical squares");
}

// -----------------------------------------------------------------------
// Test 2: Quarter overlap — partial intersection
//
// Subject = unit square [(0,0),(1,0),(1,1),(0,1)].
// Clip    = square shifted by (0.5, 0.5):
//           [(0.5,0.5),(1.5,0.5),(1.5,1.5),(0.5,1.5)].
// The overlap region is [(0.5,0.5),(1,0.5),(1,1),(0.5,1)] = 0.25 area.
// Subject area = 1.0, so expected fraction = 0.25.
// -----------------------------------------------------------------------
static void testQuarterOverlap()
{
    std::cout << "\n--- Test 2: Quarter overlap ---\n";

    using namespace RastClipper;

    PathD subject = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    PathD clip = {{0.5, 0.5}, {1.5, 0.5}, {1.5, 1.5}, {0.5, 1.5}};

    // Verify areas
    checkClose(Area(subject), 1.0, 1e-12, "Subject area = 1.0");
    checkClose(Area(clip), 1.0, 1e-12, "Clip area = 1.0");

    PathsD clips = {clip};
    const int NBIT = 256;
    auto result = IntersectionFractions(subject, clips, NBIT);

    check(result.fractions.size() == 1, "Returned 1 fraction");
    check(result.centroids.size() == 1, "Returned 1 centroid");

    // At NBIT=256, rasterization error is ~O(1/NBIT) ≈ 0.4%.
    // Use 5% tolerance to be safe against edge-erosion effects.
    checkClose(result.fractions[0], 0.25, 0.05,
               "Fraction ~ 0.25 for quarter overlap");

    // Overlap region is [(0.5,0.5),(1,0.5),(1,1),(0.5,1)].
    // Its centroid is (0.75, 0.75).
    checkClose(result.centroids[0].x, 0.75, 0.05,
               "Centroid x ~ 0.75 for quarter overlap");
    checkClose(result.centroids[0].y, 0.75, 0.05,
               "Centroid y ~ 0.75 for quarter overlap");
}

// -----------------------------------------------------------------------
// Test 3: Mixed topologies — triangle subject vs quad + pentagon + tri
//
// Subject = triangle [(0,0),(10,0),(5,10)],  area = 50.0
//   Left edge:  x = y/2       (from (0,0) to (5,10))
//   Right edge: x = 10 - y/2  (from (10,0) to (5,10))
//
// Three non-overlapping clips, all fully inside the subject:
//
// Clip 0 = quad [(1,1),(4,1),(4,3),(1,3)],  area = 6.0
//   At y=1: subject x in [0.5, 9.5]; quad x in [1, 4] — inside.
//   At y=3: subject x in [1.5, 8.5]; quad x in [1, 4] — inside.
//   Expected fraction = 6 / 50 = 0.12
//
// Clip 1 = pentagon [(5,1),(8,1),(8,3),(7,4),(5,3)],  area = 7.5
//   At y=1: [5,8] inside [0.5, 9.5].
//   At y=3: [5,8] inside [1.5, 8.5].
//   At y=4: x=7 inside [2, 8].
//   Non-overlapping with clip 0 (x ranges disjoint: [1,4] vs [5,8]).
//   Expected fraction = 7.5 / 50 = 0.15
//
// Clip 2 = triangle [(3,5),(5,5),(4,7)],  area = 2.0
//   At y=5: [3,5] inside [2.5, 7.5].
//   At y=7: x=4 inside [3.5, 6.5].
//   Non-overlapping with clips 0,1 (y ranges disjoint: [5,7] vs [1,4]).
//   Expected fraction = 2 / 50 = 0.04
//
// Total = 0.12 + 0.15 + 0.04 = 0.31
// -----------------------------------------------------------------------
static void testMixedTopologies()
{
    std::cout << "\n--- Test 3: Mixed topologies (tri vs quad+pent+tri) ---\n";

    using namespace RastClipper;

    PathD subject = {{0, 0}, {10, 0}, {5, 10}};
    checkClose(Area(subject), 50.0, 1e-12, "Triangle subject area = 50.0");

    PathD clipQuad = {{1, 1}, {4, 1}, {4, 3}, {1, 3}};
    PathD clipPent = {{5, 1}, {8, 1}, {8, 3}, {7, 4}, {5, 3}};
    PathD clipTri = {{3, 5}, {5, 5}, {4, 7}};

    checkClose(Area(clipQuad), 6.0, 1e-12, "Quad clip area = 6.0");
    checkClose(Area(clipPent), 7.5, 1e-12, "Pentagon clip area = 7.5");
    checkClose(Area(clipTri), 2.0, 1e-12, "Triangle clip area = 2.0");

    PathsD clips = {clipQuad, clipPent, clipTri};
    const int NBIT = 256;
    auto result = IntersectionFractions(subject, clips, NBIT);

    check(result.fractions.size() == 3, "Returned 3 fractions");
    check(result.centroids.size() == 3, "Returned 3 centroids");

    checkClose(result.fractions[0], 0.12, 0.02,
               "Quad fraction ~ 0.12");
    checkClose(result.fractions[1], 0.15, 0.02,
               "Pentagon fraction ~ 0.15");
    checkClose(result.fractions[2], 0.04, 0.02,
               "Small triangle fraction ~ 0.04");

    const double total =
        result.fractions[0] + result.fractions[1] + result.fractions[2];
    checkClose(total, 0.31, 0.03,
               "Sum of fractions ~ 0.31");

    // Clip 0 = quad [(1,1),(4,1),(4,3),(1,3)]: centroid at (2.5, 2.0)
    checkClose(result.centroids[0].x, 2.5, 0.3,
               "Quad centroid x ~ 2.5");
    checkClose(result.centroids[0].y, 2.0, 0.3,
               "Quad centroid y ~ 2.0");

    // Clip 1 = pentagon [(5,1),(8,1),(8,3),(7,4),(5,3)]: centroid at (6.6, 2.2)
    // (area-weighted centroid of the pentagon)
    checkClose(result.centroids[1].x, 6.6, 0.5,
               "Pentagon centroid x ~ 6.6");
    checkClose(result.centroids[1].y, 2.2, 0.5,
               "Pentagon centroid y ~ 2.2");

    // Clip 2 = triangle [(3,5),(5,5),(4,7)]: centroid at (4.0, 5.67)
    checkClose(result.centroids[2].x, 4.0, 0.3,
               "Triangle centroid x ~ 4.0");
    checkClose(result.centroids[2].y, 5.67, 0.3,
               "Triangle centroid y ~ 5.67");

    std::cout << "  INFO: fractions = [" << result.fractions[0] << ", "
              << result.fractions[1] << ", " << result.fractions[2] << "]\n";
    std::cout << "  INFO: centroids = [(" << result.centroids[0].x << ", "
              << result.centroids[0].y << "), (" << result.centroids[1].x
              << ", " << result.centroids[1].y << "), ("
              << result.centroids[2].x << ", " << result.centroids[2].y
              << ")]\n";
}

// -----------------------------------------------------------------------

int main()
{
    std::cout << "=== RastClipper Tests ===\n";

    testIdenticalSquares();
    testQuarterOverlap();
    testMixedTopologies();

    std::cout << "\n=== Summary: " << (testCount - failCount) << "/"
              << testCount << " passed ===\n";

    return failCount > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
