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
#ifndef __MEANSHIFT_VECTOR_SPACE_HPP__
#define __MEANSHIFT_VECTOR_SPACE_HPP__

template<typename _Scalar, size_t N>
struct Rn;

template<typename _Scalar, size_t _N>
struct chart_traits<Rn<_Scalar, _N> >
{
  typedef Eigen::Matrix<_Scalar,_N,1> TangentVector;
  typedef Eigen::Matrix<_Scalar,_N,1> EmbeddedVector;
  typedef _Scalar Scalar;
};

template<typename _Scalar, size_t N> 
struct Rn
{
  static typename chart_traits<Rn>::EmbeddedVector 
  expm(const typename chart_traits<Rn>::EmbeddedVector& x, 
       const typename chart_traits<Rn>::TangentVector& u)
  {
    return x+u;
  }

  static typename chart_traits<Rn>::TangentVector 
  logm(const typename chart_traits<Rn>::EmbeddedVector& u, 
       const typename chart_traits<Rn>::EmbeddedVector& v) 
  {
    return v-u;
  }
};

template <typename _Scalar>
using R1 = Rn<_Scalar,1>;

template <typename _Scalar>
using R2 = Rn<_Scalar,2>;

template <typename _Scalar>
using R3 = Rn<_Scalar,3>;

template <typename _Scalar>
using R4 = Rn<_Scalar,4>;

template <typename _Scalar>
using R5 = Rn<_Scalar,5>;

template <typename _Scalar>
using R6 = Rn<_Scalar,6>;

#endif
