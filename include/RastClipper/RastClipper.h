// File       : RastClipper.h
// Created    : Sun Mar 16 2026
// Author     : Mhamad Mahdi Alloush
// Description: Rasterization-based polygon clipping library.
//              Mirrors Clipper2 type names (PointD, PathD, PathsD) so that
//              switching between Clipper2 and RastClipper is straightforward.
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef RAST_CLIPPER_H
#define RAST_CLIPPER_H

#include <vector>

namespace RastClipper
{

struct PointD
{
    double x, y;
};

using PathD = std::vector<PointD>;
using PathsD = std::vector<PathD>;

// Signed area of a polygon (shoelace formula).
// Positive for counter-clockwise winding.
double Area(const PathD& path);

// Given a subject polygon and M clip polygons (all 2-D), returns M area
// fractions where:
//   fraction[i] = rasterized_overlap(subject, clips[i])
//               / rasterized_interior_area(subject)
//
// The rasterization uses an NBIT x NBIT bitmap with depth-buffer compositing
// and 4-neighbor interior-pixel erosion.
//
// @param subject    The subject polygon (arbitrary vertex count).
// @param clips      Vector of clip polygons (arbitrary vertex count each).
// @param resolution Bitmap resolution (default 128).
// @return           Vector of area fractions, one per clip polygon.
std::vector<double> IntersectionFractions(const PathD& subject,
                                          const PathsD& clips,
                                          int resolution = 128);

} // namespace RastClipper

#endif // RAST_CLIPPER_H
