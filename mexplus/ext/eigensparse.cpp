#include <stdio.h>
#include <mex.h>
#include <array>
#include <algorithm>
#include <math.h>
#include "matrix.h"

#include <mexplus.h>

#include "mexplus_ext.hpp"

#include <eigen3/Eigen/Sparse>

using namespace mexplus;

MEX_DEFINE(new) (int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) 
{
    InputArguments input(nrhs, prhs, 3);
    OutputArguments output(nlhs, plhs, 1);

    Eigen::SparseMatrix<double, Eigen::RowMajor> sm1, sm2, sm3;

    mexplus::MxArray::to<Eigen::SparseMatrix<double,Eigen::RowMajor>>(prhs[0],&sm1);
    mexplus::MxArray::to<Eigen::SparseMatrix<double,Eigen::RowMajor>>(prhs[1],&sm2);
    mexplus::MxArray::to<Eigen::SparseMatrix<double,Eigen::RowMajor>>(prhs[2],&sm3);

    auto tempop1 = sm1 + sm2;
    auto tempop2 = tempop1 - sm3;

    mxArray* array = mexplus::MxArray::from(tempop2.eval());
    output.set(0, array);
}


MEX_DISPATCH
