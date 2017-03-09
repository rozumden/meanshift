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
/** MEX dispatch helper library.
 *
 * Copyright 2013 Kota Yamaguchi.
 */

#include "mexplus/dispatch.h"

#define EXPECT(condition) if (!(condition)) \
    mexErrMsgTxt(#condition " not true.")
#define PRINTF(...) mexPrintf(__VA_ARGS__); \
    mexCallMATLAB(0, NULL, 0, NULL, "drawnow")

namespace {

MEX_DEFINE(foo) (int nlhs,
                 mxArray* plhs[],
                 int nrhs,
                 const mxArray* prhs[]) {
}

MEX_DEFINE(bar) (int nlhs,
                 mxArray* plhs[],
                 int nrhs,
                 const mxArray* prhs[]) {
}

}  // namespace

MEX_DISPATCH
