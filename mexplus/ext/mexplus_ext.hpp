#ifndef MEANSHIFT_MEXPLUS_EXT_HPP_
#define MEANSHIFT_MEXPLUS_EXT_HPP_

#include <mexplus.h>
#include <eigen3/Eigen/Sparse>

#include "sparse_matlab_eigen.hpp"

namespace mexplus {
  // Define template specializations for doubles.
  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<double, Eigen::ColMajor>& u) 
  {
    return eigen2matlab<double, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<double, Eigen::ColMajor>* u) 
  {
    matlab2eigen<double, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<double, Eigen::RowMajor>& u) 
  {
    return eigen2matlab<double,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<double, Eigen::RowMajor>* u) 
  {
    matlab2eigen<double,Eigen::RowMajor>(array,u);
  }

  // Define template specializations for ints.
  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<int, Eigen::ColMajor>& u) 
  {
    return eigen2matlab<int, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<int, Eigen::ColMajor>* u) 
  {
    matlab2eigen<int, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<int, Eigen::RowMajor>& u) 
  {
    return eigen2matlab<int,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<int, Eigen::RowMajor>* u) 
  {
    matlab2eigen<int,Eigen::RowMajor>(array,u);
  }

  // Define template specializations for long long.
  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<long long, Eigen::ColMajor>& u) 
  {
    return eigen2matlab<long long, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<long long, Eigen::ColMajor>* u) 
  {
    matlab2eigen<long long, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<long long, Eigen::RowMajor>& u) 
  {
    return eigen2matlab<long long,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<long long, Eigen::RowMajor>* u) 
  {
    matlab2eigen<long long,Eigen::RowMajor>(array,u);
  }

}

#endif
