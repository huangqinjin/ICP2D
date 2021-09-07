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


namespace ICP2D
{

}


#endif //ICP2D_HPP
