//
// Copyright (c) 2020-2021 Huang Qinjin (huangqinjin@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)
//
#ifndef ICP2D_HPP
#define ICP2D_HPP

#if defined(_MSC_VER)
#  define ICP2D_EXPORT __declspec(dllexport)
#  define ICP2D_IMPORT __declspec(dllimport)
#else
#  define ICP2D_EXPORT __attribute__((visibility("default")))
#  define ICP2D_IMPORT __attribute__((visibility("default")))
#endif

#if defined(ICP2D_EXPORTS)
#  define ICP2D_API ICP2D_EXPORT
#elif ICP2D_SHARED_LIBRARY
#  define ICP2D_API ICP2D_IMPORT
#else
#  define ICP2D_API
#endif


#include <cstdio>
#include <vector>
#include <ostream>
#include <Eigen/Geometry>

namespace ICP2D
{
    using Scaling = Eigen::UniformScaling<double>;
    using Rotation = Eigen::Rotation2Dd;
    using Translation = Eigen::Translation2d;
    using Transform = Eigen::Affine2d;

    using Point = Eigen::Vector2d;
    using PointSet = std::vector<Point, Eigen::aligned_allocator<Point>>;
    using WeightVector = Eigen::VectorXd;
    using BoundingBox = Eigen::AlignedBox2d;

    struct Sim2D
    {
        double s;
        double r;
        double x;
        double y;

        Scaling scaling() const noexcept
        {
            return Scaling(s);
        }

        Rotation rotation() const noexcept
        {
            return Rotation(r);
        }

        Translation translation() const noexcept
        {
            return Translation(x, y);
        }

        Transform transform() const noexcept
        {
            return translation() * rotation() * scaling();
        }

        static Sim2D Identity()
        {
            Sim2D T;
            T.s = 1;
            T.r = 0;
            T.x = 0;
            T.y = 0;
            return T;
        }

        Sim2D operator*(const Sim2D& other) const noexcept
        {
            Eigen::Vector2d t = transform() * other.translation().vector();
            Sim2D T;
            T.s = s * other.s;
            T.r = r + other.r;
            T.x = t.x();
            T.y = t.y();
            return T;
        }
    };

    template<class CharT, class Traits>
    std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const Sim2D& T)
    {
        return os
            << '{'
            << 's' << ':' << T.s << ',' << ' '
            << 'r' << ':' << T.r << ',' << ' '
            << 'x' << ':' << T.x << ',' << ' '
            << 'y' << ':' << T.y
            << '}'
            ;
    }

    ICP2D_API Sim2D solve(const PointSet& src, const PointSet& dst, const WeightVector& w);
    ICP2D_API double error(const Sim2D& T, const PointSet& src, const PointSet& dst, const WeightVector& w);

    class ICP2D_API SVG
    {
        FILE* out;
        Sim2D T;
        Point scale;
        BoundingBox view;

    public:
        SVG() noexcept;
        ~SVG();
        void open(const char* file, const Point& scale = Point::Ones()) noexcept;
        void close() noexcept;
        void push(const Sim2D& transform) noexcept;
        void pop() noexcept;
        void draw(const PointSet& points, double radius, const char* color) noexcept;
    };
}


#endif //ICP2D_HPP
