// File       : RastClipper.cpp
// Created    : Sun Mar 16 2026
// Author     : Mhamad Mahdi Alloush
// Description: Rasterization-based polygon clipping implementation.
//              Scan-line rasterization with depth-buffer compositing.
// Copyright (c) 2024 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "RastClipper/RastClipper.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace RastClipper
{

// -----------------------------------------------------------------------
// Shoelace area
// -----------------------------------------------------------------------

double Area(const PathD& path)
{
    const int n = static_cast<int>(path.size());
    if (n < 3)
        return 0.0;

    double area = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const int j = (i + 1) % n;
        area += path[i].x * path[j].y;
        area -= path[j].x * path[i].y;
    }
    return 0.5 * area;
}

namespace
{

// Depth test modes
constexpr int ZFLESS = 2;
constexpr int ZFGEQU = 7;

// Internal vertex: normalised x, y in [-1,1] plus depth
struct Vertex3
{
    double x, y, depth;
};

// Pixel-space vertex (1-based) with original depth
struct PixelVertex
{
    int px, py;
    double depth;
};

// Edge table entry for scan-line fill
struct EdgeEntry
{
    int yMin, yMax;
    double x, dx, d, dd;
};

// -----------------------------------------------------------------------
// Gross clipping: trivial accept/reject against [-1,1]^3
// Returns: 0 = reject, numPts = accept, 1 = partial
// Trivial accept/reject against the [-1,1]^3 viewport
// -----------------------------------------------------------------------

int grossClip(const Vertex3* pts, int numPts)
{
    int lc = 0, rc = 0, bc = 0, tc = 0, nc = 0, fc = 0;

    for (int i = 0; i < numPts; ++i)
    {
        if (pts[i].x < -1.0)
            ++lc;
        else if (pts[i].x > 1.0)
            ++rc;

        if (pts[i].y < -1.0)
            ++bc;
        else if (pts[i].y > 1.0)
            ++tc;

        if (pts[i].depth < -1.0)
            ++nc;
        else if (pts[i].depth > 1.0)
            ++fc;
    }

    if (lc == numPts || rc == numPts || bc == numPts || tc == numPts ||
        nc == numPts || fc == numPts)
        return 0;

    if (lc + rc + bc + tc + nc + fc == 0)
        return numPts;

    return 1;
}

// -----------------------------------------------------------------------
// Fine clipping: Sutherland-Hodgman against [-1,1]^3
// Clips a closed polygon to the [-1,1]^3 viewport
// -----------------------------------------------------------------------

std::vector<Vertex3> fineClip(const Vertex3* pts, int numPts)
{
    constexpr double EPS = 1.0001;
    constexpr double NEPS = -1.0001;

    std::vector<Vertex3> out(pts, pts + numPts);

    for (int dim = 0; dim < 3; ++dim)
    {
        std::vector<Vertex3> local;
        const int n = static_cast<int>(out.size());
        if (n == 0)
            break;

        auto getC = [dim](const Vertex3& v) -> double
        { return dim == 0 ? v.x : (dim == 1 ? v.y : v.depth); };

        for (int j = 0; j < n; ++j)
        {
            const int k = (j + 1) % n;
            const double cj = getC(out[j]);
            const double ck = getC(out[k]);

            if (cj < NEPS)
            {
                if (ck < NEPS)
                {
                    continue;
                }
                else
                {
                    double delta = (-1.0 - ck) / (cj - ck);
                    Vertex3 is;
                    is.x = out[k].x + (out[j].x - out[k].x) * delta;
                    is.y = out[k].y + (out[j].y - out[k].y) * delta;
                    is.depth =
                        out[k].depth + (out[j].depth - out[k].depth) * delta;
                    if (dim == 0)
                        is.x = -1.0;
                    else if (dim == 1)
                        is.y = -1.0;
                    else
                        is.depth = -1.0;
                    local.push_back(is);
                }
            }
            else if (cj > EPS)
            {
                if (ck > EPS)
                {
                    continue;
                }
                else
                {
                    double delta = (1.0 - ck) / (cj - ck);
                    Vertex3 is;
                    is.x = out[k].x + (out[j].x - out[k].x) * delta;
                    is.y = out[k].y + (out[j].y - out[k].y) * delta;
                    is.depth =
                        out[k].depth + (out[j].depth - out[k].depth) * delta;
                    if (dim == 0)
                        is.x = 1.0;
                    else if (dim == 1)
                        is.y = 1.0;
                    else
                        is.depth = 1.0;
                    local.push_back(is);
                }
            }
            else
            {
                local.push_back(out[j]);
            }

            if (ck < NEPS)
            {
                double delta = (-1.0 - cj) / (ck - cj);
                Vertex3 is;
                is.x = out[j].x + (out[k].x - out[j].x) * delta;
                is.y = out[j].y + (out[k].y - out[j].y) * delta;
                is.depth = out[j].depth + (out[k].depth - out[j].depth) * delta;
                if (dim == 0)
                    is.x = -1.0;
                else if (dim == 1)
                    is.y = -1.0;
                else
                    is.depth = -1.0;
                local.push_back(is);
            }
            else if (ck > EPS)
            {
                double delta = (1.0 - cj) / (ck - cj);
                Vertex3 is;
                is.x = out[j].x + (out[k].x - out[j].x) * delta;
                is.y = out[j].y + (out[k].y - out[j].y) * delta;
                is.depth = out[j].depth + (out[k].depth - out[j].depth) * delta;
                if (dim == 0)
                    is.x = 1.0;
                else if (dim == 1)
                    is.y = 1.0;
                else
                    is.depth = 1.0;
                local.push_back(is);
            }
        }

        out = std::move(local);
    }

    return out;
}

// -----------------------------------------------------------------------
// Scan line: fill pixel spans with depth test
// All pixel coordinates are 1-based.
// Bitmap indexed as (x-1) + (y-1)*IRES.
// -----------------------------------------------------------------------

void scanLine(std::vector<int>& bitcol,
              std::vector<float>& bitdep,
              int IRES,
              const int* ints,
              const double* intd,
              int numint,
              int xy,
              int zmode,
              int colour)
{
    for (int pos = 0; pos < numint - 1; pos += 2)
    {
        int s1 = ints[pos];
        int s2 = ints[pos + 1];
        double d1 = intd[pos];
        double d2 = intd[pos + 1];

        double ddep;
        if (s1 == s2)
        {
            ddep = 0.0;
            s2 = s2 + 1;
        }
        else
        {
            ddep = (d2 - d1) / static_cast<double>(s2 - s1);
            if (s2 == IRES)
                s2 = s2 + 1;
        }

        if (zmode == ZFLESS)
        {
            for (int i = s1; i < s2; ++i)
            {
                const double depth = d1 + static_cast<double>(i - s1) * ddep;
                const int idx = (i - 1) + (xy - 1) * IRES;
                if (depth < static_cast<double>(bitdep[idx]))
                {
                    bitcol[idx] = colour;
                    bitdep[idx] = static_cast<float>(depth);
                }
            }
        }
        else
        {
            for (int i = s1; i < s2; ++i)
            {
                const double depth = d1 + static_cast<double>(i - s1) * ddep;
                const int idx = (i - 1) + (xy - 1) * IRES;
                bitcol[idx] = colour;
                bitdep[idx] = static_cast<float>(depth);
            }
        }
    }
}

// -----------------------------------------------------------------------
// Delete duplicate consecutive pixel-space vertices
// -----------------------------------------------------------------------

void delDuplVx(std::vector<PixelVertex>& verts, int& numPt)
{
    if (numPt <= 1)
        return;

    int write = 0;
    for (int i = 0; i < numPt; ++i)
    {
        const int next = (i + 1) % numPt;
        if (verts[i].px != verts[next].px || verts[i].py != verts[next].py)
        {
            verts[write] = verts[i];
            ++write;
        }
    }
    numPt = write;
}

// -----------------------------------------------------------------------
// Scan polygon: edge table construction + scan-line fill
// Rasterizes a polygon via scan-line fill
// -----------------------------------------------------------------------

void scanPolygon(std::vector<int>& bitcol,
                 std::vector<float>& bitdep,
                 int IRES,
                 const Vertex3* verts,
                 int numPt,
                 int zmode,
                 int colour)
{
    if (numPt < 3)
        return;

    // Scale from normalised [-1,1] to 1-based pixel coords, clamped to [1,IRES]
    std::vector<PixelVertex> pv(numPt);
    int xMin = IRES + 1, xMax = 0;
    int yMin = IRES + 1, yMax = 0;
    double dxMin = 0.0, dxMax = 0.0;

    for (int i = 0; i < numPt; ++i)
    {
        int px = static_cast<int>(
            std::lround((verts[i].x + 1.0) * (IRES - 1) * 0.5 + 1.0));
        int py = static_cast<int>(
            std::lround((verts[i].y + 1.0) * (IRES - 1) * 0.5 + 1.0));
        px = std::max(1, std::min(IRES, px));
        py = std::max(1, std::min(IRES, py));

        pv[i] = {px, py, verts[i].depth};

        if (px < xMin)
        {
            xMin = px;
            dxMin = verts[i].depth;
        }
        if (px > xMax)
        {
            xMax = px;
            dxMax = verts[i].depth;
        }
        if (py < yMin)
            yMin = py;
        if (py > yMax)
            yMax = py;
    }

    // Delete duplicate vertices
    delDuplVx(pv, numPt);

    // Degenerate: single scan line or too few vertices
    if (numPt < 3 || yMin == yMax)
    {
        int ints[2] = {xMin, xMax};
        double intd[2] = {dxMin, dxMax};
        scanLine(bitcol, bitdep, IRES, ints, intd, 2, yMin, zmode, colour);
        return;
    }

    // Build edge table, ignoring horizontal edges
    std::vector<EdgeEntry> edgeTable;
    edgeTable.reserve(numPt);
    bool shLast = false;

    for (int i = 0; i < numPt; ++i)
    {
        const int j = (i + 1) % numPt;
        const int k = (i - 1 + numPt) % numPt;

        const int yi = pv[i].py;
        const int yj = pv[j].py;
        const int yk = pv[k].py;

        if (yi < yj)
        {
            EdgeEntry e;
            e.yMax = yj;
            e.dx = static_cast<double>(pv[j].px - pv[i].px) /
                   static_cast<double>(yj - yi);
            e.dd = (pv[j].depth - pv[i].depth) / static_cast<double>(yj - yi);

            if (yi <= yk)
            {
                e.yMin = yi;
                e.x = static_cast<double>(pv[i].px);
                e.d = pv[i].depth;
            }
            else
            {
                e.yMin = yi + 1;
                e.x = static_cast<double>(pv[i].px) + e.dx;
                e.d = pv[i].depth + e.dd;
            }
            edgeTable.push_back(e);
        }
        else if (yi > yj)
        {
            EdgeEntry e;
            e.yMin = yj;
            e.x = static_cast<double>(pv[j].px);
            e.dx = static_cast<double>(pv[i].px - pv[j].px) /
                   static_cast<double>(yi - yj);
            e.d = pv[j].depth;
            e.dd = (pv[i].depth - pv[j].depth) / static_cast<double>(yi - yj);

            if (yi > yk || yi == yMax)
            {
                e.yMax = yi;
            }
            else
            {
                e.yMax = yi - 1;
            }
            edgeTable.push_back(e);
        }
        else
        {
            // Horizontal edge
            if (yi > yk && yi < yMax)
            {
                if (!edgeTable.empty())
                {
                    edgeTable.back().yMax = yi - 1;
                }
                else
                {
                    shLast = true;
                }
            }
        }
    }

    if (shLast && !edgeTable.empty())
    {
        edgeTable.back().yMax = pv[0].py - 1;
    }

    // Sort edge table by yMin
    std::sort(edgeTable.begin(),
              edgeTable.end(),
              [](const EdgeEntry& a, const EdgeEntry& b)
    { return a.yMin < b.yMin; });

    const int numET = static_cast<int>(edgeTable.size());
    if (numET == 0)
        return;

    // Walk scan lines, find intersections, fill spans
    std::vector<int> ints(numET + 2);
    std::vector<double> intd(numET + 2);

    for (int scanY = yMin; scanY <= yMax; ++scanY)
    {
        int numint = 0;
        for (int ej = 0; ej < numET; ++ej)
        {
            if (scanY < edgeTable[ej].yMin)
                break;
            if (scanY > edgeTable[ej].yMax)
                continue;

            ints[numint] = static_cast<int>(
                std::lround(edgeTable[ej].x +
                            static_cast<double>(scanY - edgeTable[ej].yMin) *
                                edgeTable[ej].dx));
            intd[numint] = edgeTable[ej].d +
                           static_cast<double>(scanY - edgeTable[ej].yMin) *
                               edgeTable[ej].dd;
            ++numint;
        }

        // Sort intersection points by x (insertion sort)
        for (int a = 1; a < numint; ++a)
        {
            int tmpI = ints[a];
            double tmpD = intd[a];
            int b = a - 1;
            while (b >= 0 && ints[b] > tmpI)
            {
                ints[b + 1] = ints[b];
                intd[b + 1] = intd[b];
                --b;
            }
            ints[b + 1] = tmpI;
            intd[b + 1] = tmpD;
        }

        scanLine(bitcol,
                 bitdep,
                 IRES,
                 ints.data(),
                 intd.data(),
                 numint,
                 scanY,
                 zmode,
                 colour);
    }
}

// -----------------------------------------------------------------------
// Draw polygon: gross clip -> fine clip -> scan fill
// Orchestrates clipping and scan-line fill
// -----------------------------------------------------------------------

void drawPolygon(std::vector<int>& bitcol,
                 std::vector<float>& bitdep,
                 int IRES,
                 const Vertex3* points,
                 int numPts,
                 int zmode,
                 int colour)
{
    if (numPts < 3)
        return;

    const int gc = grossClip(points, numPts);
    if (gc == 0)
        return;

    if (gc < numPts)
    {
        auto clipped = fineClip(points, numPts);
        const int cn = static_cast<int>(clipped.size());
        if (cn > 2)
        {
            scanPolygon(
                bitcol, bitdep, IRES, clipped.data(), cn, zmode, colour);
        }
    }
    else
    {
        scanPolygon(bitcol, bitdep, IRES, points, numPts, zmode, colour);
    }
}

} // anonymous namespace

// -----------------------------------------------------------------------
// IntersectionFractions
// -----------------------------------------------------------------------

IntersectionResult
IntersectionFractions(const PathD& subject, const PathsD& clips, int resolution)
{
    const int IRES = resolution;
    const int IRES2 = IRES * IRES;
    const int M = static_cast<int>(clips.size());

    IntersectionResult result;
    result.fractions.resize(M, 0.0);
    result.centroids.resize(M, {0.0, 0.0});

    if (subject.size() < 3 || M == 0)
        return result;

    // Subject bounding box
    double uMin = subject[0].x, uMax = subject[0].x;
    double vMin = subject[0].y, vMax = subject[0].y;
    for (const auto& p : subject)
    {
        uMin = std::min(uMin, p.x);
        uMax = std::max(uMax, p.x);
        vMin = std::min(vMin, p.y);
        vMax = std::max(vMax, p.y);
    }

    // Viewport mapping with 5% padding
    const double uCenter = 0.5 * (uMin + uMax);
    const double vCenter = 0.5 * (vMin + vMax);
    const double uHalfRange = 1.05 * 0.5 * (uMax - uMin);
    const double vHalfRange = 1.05 * 0.5 * (vMax - vMin);

    if (uHalfRange < 1e-30 || vHalfRange < 1e-30)
        return result;

    auto toNorm = [&](double u, double v, double depth) -> Vertex3
    { return {(u - uCenter) / uHalfRange, (v - vCenter) / vHalfRange, depth}; };

    // Allocate bitmaps: (x-1) + (y-1)*IRES, 1-based coords
    // Colour: 0 = background, 1 = subject, i+2 = clip[i]
    std::vector<int> bitcol(IRES2, 0);
    std::vector<float> bitdep(IRES2, -1.0f);

    // Draw subject at depth=1.0, unconditional (ZFGEQU)
    {
        const int n = static_cast<int>(subject.size());
        std::vector<Vertex3> sv(n);
        for (int i = 0; i < n; ++i)
            sv[i] = toNorm(subject[i].x, subject[i].y, 1.0);

        drawPolygon(bitcol, bitdep, IRES, sv.data(), n, ZFGEQU, 1);
    }

    // Draw each clip at depth=0.0, ZFLESS
    for (int ci = 0; ci < M; ++ci)
    {
        const int n = static_cast<int>(clips[ci].size());
        if (n < 3)
            continue;

        std::vector<Vertex3> cv(n);
        for (int i = 0; i < n; ++i)
            cv[i] = toNorm(clips[ci][i].x, clips[ci][i].y, 0.0);

        drawPolygon(bitcol, bitdep, IRES, cv.data(), n, ZFLESS, ci + 2);
    }

    // Count interior pixels with 4-neighbor erosion and accumulate
    // pixel-weighted centroids per clip colour.
    // 1-based loop: J=2..IRES-1 (y), I=2..IRES-1 (x)
    int icount = 0;
    std::vector<int> colourCounts(M, 0);
    std::vector<double> centroidSumI(M, 0.0);
    std::vector<double> centroidSumJ(M, 0.0);

    for (int j = 2; j <= IRES - 1; ++j)
    {
        for (int i = 2; i <= IRES - 1; ++i)
        {
            const int idx = (i - 1) + (j - 1) * IRES;
            if (bitcol[idx] == 0)
                continue;

            if (bitcol[(i - 2) + (j - 1) * IRES] == 0)
                continue;
            if (bitcol[i + (j - 1) * IRES] == 0)
                continue;
            if (bitcol[(i - 1) + (j - 2) * IRES] == 0)
                continue;
            if (bitcol[(i - 1) + j * IRES] == 0)
                continue;

            ++icount;
            const int col = bitcol[idx];
            if (col >= 2 && col - 2 < M)
            {
                const int ci = col - 2;
                ++colourCounts[ci];
                centroidSumI[ci] += static_cast<double>(i);
                centroidSumJ[ci] += static_cast<double>(j);
            }
        }
    }

    if (icount > 0)
    {
        const double inv = 1.0 / static_cast<double>(icount);
        for (int i = 0; i < M; ++i)
        {
            result.fractions[i] = static_cast<double>(colourCounts[i]) * inv;

            if (colourCounts[i] > 0)
            {
                // Average pixel position (1-based)
                const double avgI =
                    centroidSumI[i] / static_cast<double>(colourCounts[i]);
                const double avgJ =
                    centroidSumJ[i] / static_cast<double>(colourCounts[i]);

                // Pixel (1-based) -> normalised [-1,1]
                const double normU =
                    2.0 * (avgI - 1.0) / static_cast<double>(IRES - 1) - 1.0;
                const double normV =
                    2.0 * (avgJ - 1.0) / static_cast<double>(IRES - 1) - 1.0;

                // Normalised [-1,1] -> input coordinates
                result.centroids[i].x = uCenter + normU * uHalfRange;
                result.centroids[i].y = vCenter + normV * vHalfRange;
            }
        }
    }

    return result;
}

} // namespace RastClipper
