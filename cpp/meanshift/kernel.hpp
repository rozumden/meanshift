// Copyright (c) 2017 James Pritts, Denys Rozumnyi
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in 
// all copies or substantial portions of the Software.
//
// The Software is provided "as is", without warranty of any kind.
#ifndef __MEANSHIFT_KERNEL_HPP__
#define __MEANSHIFT_KERNEL_HPP__ 

#include "meanshift/manifold/chart.hpp"
#include "meanshift/profile.hpp"

template <typename Chart,typename ...T>
struct Kernel;

template <typename T>
struct kernel_traits;

template <typename _Chart,typename ...T>
struct kernel_traits<Kernel<_Chart,T...> >
{
  typedef _Chart Chart;
  typedef typename chart_traits<_Chart>::EmbeddedVector Vector;
  typedef typename chart_traits<_Chart>::Scalar Scalar;
};

template<typename Chart,typename Shape>
struct Kernel<Chart,Profile<Shape> >
{
  static std::tuple<typename chart_traits<Chart>::Scalar,
		    typename chart_traits<Chart>::TangentVector>
  K(const typename kernel_traits<Kernel>::Vector& x, 
    const typename kernel_traits<Kernel>::Vector& x_i,
    const typename kernel_traits<Kernel>::Scalar h) 
  {
    const typename chart_traits<Chart>::Scalar d = 
      chart_traits<Chart>::TangentVector::RowsAtCompileTime;

    typename chart_traits<Chart>::TangentVector u = Chart::logm(x,x_i);
    typename chart_traits<Chart>::Scalar dist_sq = u.adjoint()*u; 
    typename chart_traits<Chart>::Scalar w_i = 
      Profile<Shape>::k(dist_sq/(h*h))/pow(h,d);

    return std::make_tuple(w_i,w_i*u);
  }

  static std::tuple<typename chart_traits<Chart>::Scalar,
		    typename chart_traits<Chart>::TangentVector>
  G(const typename kernel_traits<Kernel>::Vector& x, 
    const typename kernel_traits<Kernel>::Vector& x_i,
    const typename kernel_traits<Kernel>::Scalar h)
  {
    const typename chart_traits<Chart>::Scalar d = 
      chart_traits<Chart>::TangentVector::RowsAtCompileTime;

    typename chart_traits<Chart>::TangentVector u = Chart::logm(x,x_i);
    typename chart_traits<Chart>::Scalar dist_sq = u.adjoint()*u;
    typename chart_traits<Chart>::Scalar w_i = Profile<Shape>::g(dist_sq/(h*h))/pow(h,d+2);

    return std::make_tuple(w_i,w_i*u);
  }
};


#endif
