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
#ifndef MEANSHIFT_SPARSE_MATLAB_EIGEN_HPP_
#define MEANSHIFT_SPARSE_MATLAB_EIGEN_HPP_

#include <eigen3/Eigen/Sparse>

template <typename Scalar, int Majority>
void 
matlab_to_eigen(const mxArray* array, Eigen::SparseMatrix<Scalar,Majority>* u) 
{
  mwSize m = mxGetM(array);
  mwSize n = mxGetN(array);
  mwSize nnz = mxGetNzmax(array);
  u->resize(m,n);
  
  mwIndex *ir,*jc;
  double *pr;
  pr = mxGetPr(array);
  ir = mxGetIr(array);
  jc = mxGetJc(array);
  
  // Filling a sparse matrix
  typedef Eigen::Triplet<Scalar> T;
  std::vector<T> tripletList;
  tripletList.reserve(nnz);
  int j = 0;
  for(int k=0; k<nnz; k++)
  {
      while (jc[j+1] <= k) {
        j++;
        if(j > n+1) {
          j = -1;
          break;
        }
      }
      if (j == -1) { break; }
      tripletList.push_back(T(ir[k],j,pr[k]));
  }
  u->setFromTriplets(tripletList.begin(), tripletList.end());
}

template <typename Scalar, int Majority>
mxArray*  
eigen_to_matlab(const Eigen::SparseMatrix<Scalar,Majority>& u) 
{
  mwSize m,n,nzmax;
  m = u.rows();
  n = u.cols();
  nzmax = u.nonZeros();
  mxArray *array = mxCreateSparse(m,n,nzmax,mxREAL);

  Eigen::SparseMatrix<Scalar,Eigen::ColMajor> ucol(u);

  mwIndex *ir,*jc;
  double *pr;
  pr = mxGetPr(array);
  ir = mxGetIr(array);
  jc = mxGetJc(array);
  int nel = 0; 
  for (int k=0; k<ucol.outerSize(); ++k) {
    jc[k] = nel;
    for (typename Eigen::SparseMatrix<Scalar,Eigen::ColMajor>::InnerIterator it(ucol,k); it; ++it) {
      pr[nel] = it.value();
      ir[nel] = it.row();
      nel++;
    }
  }
  jc[n] = nel;

  return array;
}

namespace mexplus  {

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<double, Eigen::ColMajor>& u) 
  {
    return ::eigen_to_matlab<double, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<double, Eigen::ColMajor>* u) 
  {
    matlab_to_eigen<double, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<double, Eigen::RowMajor>& u) 
  {
    return eigen_to_matlab<double,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<double, Eigen::RowMajor>* u) 
  {
    matlab_to_eigen<double,Eigen::RowMajor>(array,u);
  }

  // Define template specializations for ints.
  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<int, Eigen::ColMajor>& u) 
  {
    return eigen_to_matlab<int, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<int, Eigen::ColMajor>* u) 
  {
    matlab_to_eigen<int, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<int, Eigen::RowMajor>& u) 
  {
    return eigen_to_matlab<int,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<int, Eigen::RowMajor>* u) 
  {
    matlab_to_eigen<int,Eigen::RowMajor>(array,u);
  }

  // Define template specializations for long long.
  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<long long, Eigen::ColMajor>& u) 
  {
    return eigen_to_matlab<long long, Eigen::ColMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<long long, Eigen::ColMajor>* u) 
  {
    matlab_to_eigen<long long, Eigen::ColMajor>(array,u);
  }

  template <>
  mxArray* 
  MxArray::from(const Eigen::SparseMatrix<long long, Eigen::RowMajor>& u) 
  {
    return eigen_to_matlab<long long,Eigen::RowMajor>(u);
  }

  template <>
  void 
  MxArray::to(const mxArray* array, Eigen::SparseMatrix<long long, Eigen::RowMajor>* u) 
  {
    matlab_to_eigen<long long,Eigen::RowMajor>(array,u);
  }

} // namespace mexplus


#endif
