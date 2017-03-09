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
#ifndef __MEXPLUSPLUS_EIGEN_VECTOR_HPP__
#define __MEXPLUSPLUS_EIGEN_VECTOR_HPP__

template <typename Scalar,int M>
using VectorN = Eigen::Matrix<Scalar,M,1>; 

template <typename Scalar,int M>
mxArray*
eigen_to_matlab(const std::vector<VectorN<Scalar,M> >& u)
{
  size_t n = u.size();
  if (n > 0) {
    size_t m = u[0].rows();
    mexplus::MxArray numeric(mexplus::MxArray::Numeric<Scalar>(m,n));
    size_t j = 0;
    for (const auto& u_i : u) {
      for (size_t i=0; i<m; i++)  // rows
	numeric.set(i,j,u_i.derived().coeff(i,0));
      j++;
    }
    return numeric.release();
  }
}

template <typename Scalar,int M>
std::vector<VectorN<Scalar,M> >
matlab_to_eigen(const mxArray* array)
{
  mexplus::MxArray numeric(array);
  size_t m = numeric.rows();
  size_t n = numeric.cols();
    
  std::vector<VectorN<Scalar,M> > u;
  u.resize(n);
  size_t j = 0;
  for (auto& u_i : u) {
    for (size_t i = 0; i < m; i++) 
      u_i[i] = numeric.at<Scalar>(i,j);    
    j++;
  }

  return u;
}

namespace mexplus  {
  template <>
  mxArray* 
  MxArray::from(const std::vector<Eigen::Vector2f>& u)
  {
    return eigen_to_matlab(u);
  }

  template <>
  mxArray* 
  MxArray::from(const std::vector<Eigen::Vector2d>& u)
  {
    return eigen_to_matlab(u);
  }

  template <>
  mxArray* 
  MxArray::from(const std::vector<Eigen::Vector3f>& u)
  {
    return eigen_to_matlab(u);
  }

  template <>
  mxArray* 
  MxArray::from(const std::vector<Eigen::Vector3d>& u)
  {
    return eigen_to_matlab(u);
  }

  template <>
  mxArray* 
  MxArray::from(const std::vector<Eigen::Vector4f>& u)
  {
    return eigen_to_matlab(u);
  }

  template <>
  mxArray* 
  MxArray::from(const std::vector<Eigen::Vector4d>& u)
  {
    return eigen_to_matlab(u);
  }

  template <>
  mxArray* 
  MxArray::from(const std::vector<VectorN<double,6> >& u)
  {
    return eigen_to_matlab(u);
  }

  template <>
  mxArray* 
  MxArray::from(const std::vector<VectorN<float,6> >& u)
  {
    return eigen_to_matlab(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Eigen::Vector2f>* u) 
  {
    *u = matlab_to_eigen<typename Eigen::Vector2f::Scalar,2>(array);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Eigen::Vector2d>* u) 
  {
    *u = matlab_to_eigen<typename Eigen::Vector2d::Scalar,2>(array);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Eigen::Vector3f>* u) 
  {
    *u = matlab_to_eigen<typename Eigen::Vector3f::Scalar,3>(array);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Eigen::Vector3d>* u) 
  {
    *u = matlab_to_eigen<typename Eigen::Vector3d::Scalar,3>(array);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Eigen::Vector4f>* u) 
  {
    *u = matlab_to_eigen<typename Eigen::Vector3f::Scalar,4>(array);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Eigen::Vector4d>* u) 
  {
    *u = matlab_to_eigen<typename Eigen::Vector3d::Scalar,4>(array);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<VectorN<float,6> >* u) 
  {
    *u = matlab_to_eigen<float,6>(array);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<VectorN<double,6> >* u) 
  {
    *u = matlab_to_eigen<double,6>(array); 
  }

}

#endif
