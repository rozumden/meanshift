#ifndef MEANSHIFT_SPARSE_MATLAB_EIGEN_HPP_
#define MEANSHIFT_SPARSE_MATLAB_EIGEN_HPP_

#include <mexplus.h>
#include <eigen3/Eigen/Sparse>

template <typename Scalar, int Majority>
void 
matlab2eigen(const mxArray* array, Eigen::SparseMatrix<Scalar,Majority>* u) 
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
eigen2matlab(const Eigen::SparseMatrix<Scalar,Majority>& u) 
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

#endif