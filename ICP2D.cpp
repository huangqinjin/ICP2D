//
// Copyright (c) 2020-2021 Huang Qinjin (huangqinjin@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)
//
#include "ICP2D.hpp"
#include <cassert>
#include <new>
#include <random>
#include <numeric>
#include <algorithm>

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

std::size_t ICP2D::Sampler::population(std::size_t num, size_t max) noexcept
{
    ++max;
    if (num > max) return 0;
    if (num * 2 > max) num = max - num;
    if (num == 0) return 1;

    std::size_t result = max;
    for (std::size_t i = 2; i <= num; ++i)
    {
        result *= (max - i + 1);
        result /= i;
    }
    return result;
}

Sampler* ICP2D::Sampler::random(std::size_t num, std::size_t max) noexcept
{
    struct impl : Sampler
    {
        const std::size_t num;
        const std::size_t max;
        std::mt19937_64 eng;

        impl(std::size_t num, std::size_t max) noexcept
            : num(num), max(max), eng(0) {}
        std::size_t size() const noexcept override { return num; }
        std::size_t maximum() const noexcept override { return max; }

        void sample(std::size_t* output) noexcept override
        {
            // C++17 std::sample
            for (std::size_t i = 0, n = num, m = max + 1; n != 0; ++i)
            {
                if (std::uniform_int_distribution<std::size_t>(0, --m)(eng) < n)
                {
                    *output++ = i;
                    --n;
                }
            }
        }
    };

    return new (std::nothrow) impl(num, max);
}

Sampler* ICP2D::Sampler::ordered(std::size_t num, std::size_t max) noexcept
{
    struct impl : Sampler
    {
        const std::size_t num;
        const std::size_t max;

        impl(std::size_t num, std::size_t max) noexcept
             : num(num), max(max) { init(); }
        std::size_t size() const noexcept override { return num; }
        std::size_t maximum() const noexcept override { return max; }
        void init() noexcept { std::iota(data(), data() + num, std::size_t(0)); }
        std::size_t* data() noexcept { return reinterpret_cast<std::size_t*>(this + 1); }

        void sample(std::size_t* output) noexcept override
        {
            std::copy_n(data(), num, output);
            for (std::size_t i = num; i != 0;)
            {
                --i;
                if (data()[i] < max - (num - i - 1))
                {
                    std::iota(data() + i, data() + num, data()[i] + 1);
                    return;
                }
            }
            init();
        }
    };

    void* p = ::operator new(sizeof(impl) + sizeof(std::size_t) * num, std::nothrow);
    return p ? new (p) impl(num, max) : nullptr;
}

static void header(FILE* out, const BoundingBox& view, double size)
{
    Eigen::Vector2i
            origin = (view.min)().array().floor().cast<int>(),
            fov = (view.diagonal() + Point::Ones()).array().ceil().cast<int>();

    size = size / fov.maxCoeff();
    long w = std::lround(size * fov.x());
    long h = std::lround(size * fov.y());

    std::fprintf(out, R"SVG(<svg width="%+012ld" height="%+012ld" viewBox="%+012d %+012d %+012d %+012d" xmlns="http://www.w3.org/2000/svg">)SVG",
                 w, h, origin.x(), origin.y(), fov.x(), fov.y());
    std::fputc('\n', out);
    std::fprintf(out, R"SVG(<rect x="%+012d" y="%+012d" width="%+012d" height="%+012d"/>)SVG",
                 origin.x(), origin.y(), fov.x(), fov.y());
    std::fputc('\n', out);
}

SVG::SVG() noexcept : out(nullptr) {}
SVG::~SVG() { if (out) std::fclose(out); }

void SVG::open(const char* file, const Point& scale) noexcept
{
    T = Sim2D::Identity();
    this->scale = scale;
    view = BoundingBox(Point::Zero(), Point::Ones());
    if (out) std::fclose(out);
    out = std::fopen(file, "w");
    header(out, view, 1);
    view.setEmpty();

    std::fprintf(out, R"SVG(<g transform="scale(%f %f)" text-anchor="middle" dominant-baseline="central">)SVG",
                 scale.x(), scale.y());
    std::fputc('\n', out);
}

void SVG::close() noexcept
{
    std::fputs(R"SVG(</g>)SVG", out);
    std::fputc('\n', out);
    std::fputs(R"SVG(</svg>)SVG", out);
    std::fputc('\n', out);
    std::rewind(out);

    view.transform(
        Eigen::Translation2d(scale.cwiseProduct(view.center())) *
        Eigen::Scaling(scale * 1.2) *
        Eigen::Translation2d(-view.center()));

    header(out, view, 1000);
    std::fclose(out);
    out = nullptr;
}

void SVG::push(const Sim2D& transform) noexcept
{
    T = transform;
    std::fprintf(out, R"SVG(<g transform="translate(%f %f) rotate(%f) scale(%f)">)SVG",
                 T.x, T.y, T.r / 3.141592654 * 180.0, T.s);
    std::fputc('\n', out);
}

void SVG::pop() noexcept
{
    std::fputs(R"SVG(</g>)SVG", out);
    std::fputc('\n', out);
    T = Sim2D::Identity();
}

void SVG::draw(const PointSet& points, double radius, const char* color) noexcept
{
    Transform transform = T.transform();

    for (const Point& p : points)
    {
        std::fprintf(out, R"SVG(<circle cx="%f" cy="%f" r="%f" fill-opacity="0" stroke-width="%f" stroke="%s"/>)SVG",
                     p.x(), p.y(), radius, radius * 0.1, color);
        std::fputc('\n', out);
        view.extend(transform * p);
    }
}
