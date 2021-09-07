//
// Copyright (c) 2020-2021 Huang Qinjin (huangqinjin@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)
//
#include "ICP2D.hpp"
#include <cassert>

using namespace ICP2D;


Sim2D ICP2D::solve(const PointSet& src, const PointSet& dst, const WeightVector& w)
{
    const std::size_t n = src.size();
    assert(n == src.size());
    assert(n == dst.size());
    assert(n == w.size());

    double ws = w.sum(); // sum of weights
    Point cs = Point::Zero(); // centroid of src
    Point cd = Point::Zero(); // centroid of dst
    for (std::size_t i = 0; i < n; ++i)
    {
        cs += w[i] * src[i];
        cd += w[i] * dst[i];
    }
    cs /= ws;
    cd /= ws;

    double a = 0;
    Eigen::Matrix2d W = Eigen::Matrix2d::Zero();
    for (std::size_t i = 0; i < n; ++i)
    {
        W += w[i] * (src[i] - cs) * (dst[i] - cd).transpose();
        a += w[i] * (src[i] - cs).squaredNorm();
    }

    Point t(W(0, 0) + W(1, 1), W(0, 1) - W(1, 0));

    Sim2D T;
    T.s = t.norm() / a;
    T.r = std::atan2(t.y(), t.x());
    T.x = -cd.x();
    T.y = -cd.y();

    t = T.transform() * cs;
    T.x = -t.x();
    T.y = -t.y();
    return T;
}

double ICP2D::error(const Sim2D& T, const PointSet& src, const PointSet& dst, const WeightVector& w)
{
    const std::size_t n = src.size();
    assert(n == src.size());
    assert(n == dst.size());
    assert(n == w.size());

    double e = 0;
    Transform transform = T.transform();

    for (std::size_t i = 0; i < n; ++i)
    {
        e += 0.5 * w[i] * (transform * src[i] - dst[i]).squaredNorm();
    }

    return e;
}
