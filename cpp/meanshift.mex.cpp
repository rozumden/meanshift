#include <stdio.h>
#include <mex.h>
#include <array>
#include <algorithm>
#include <iostream>
#include <math.h>

#include "mexplus/arguments.h"
#include "mexplus/dispatch.h"
#include "mexplusplus/eigen/eigen.hpp"

#include "meanshift/common/std_ext.hpp"

#include "meanshift/meanshift.hpp"
#include "meanshift/mex_factory.hpp"

using namespace mexplus;

typedef mexplus::Session<MeanshiftMexIface> MexSession;

MEX_DEFINE(new)(int nlhs, mxArray* plhs[],
		int nrhs, const mxArray* prhs[]) 
{
  try {
    InputArguments input(nrhs, prhs, 0, 9, 
			 "scalar","manifold","profile","min_support", 
			 "bandwidth_type","bandwidth","max_bandwidth",
			 "pilot_knn");
    OutputArguments output(nlhs, plhs, 1);
  
    std::string profile = input.get<std::string>("profile","epanechnikov");
    std::transform(begin(profile),end(profile),
		   begin(profile),::tolower);

    std::string manifold = input.get<std::string>("manifold","s2");
    std::transform(begin(manifold),end(manifold),
		   begin(manifold),::tolower);

    std::string bandwidth_type = input.get<std::string>("bandwidth_type","fixed");
    std::transform(begin(bandwidth_type),end(bandwidth_type),
		   begin(bandwidth_type),::tolower);

    std::string scalar = input.get<std::string>("scalar","double");
    std::transform(begin(scalar),end(scalar),
		   begin(scalar),::tolower);

    output.set(0,
	       MexSession::create(factory[std::make_tuple(scalar,
							  manifold,
							  profile,
							  bandwidth_type)](input)));
  }
  catch(const std::exception& t) {
    std::cout << "Unable to parse parameters" << std::endl;
  }
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[],
		    int nrhs, const mxArray* prhs[]) 
{
  InputArguments input(nrhs, prhs, 1);
  OutputArguments output(nlhs, plhs, 0);
  MexSession::destroy(input.get(0));
}

MEX_DEFINE(fit_and_predict)(int nlhs, mxArray* plhs[],
			    int nrhs, const mxArray* prhs[])
{
  InputArguments input(nrhs-1, prhs, 1);
  OutputArguments output(nlhs, plhs, 6);
  MeanshiftMexIface *meanshift = MexSession::get(input.get(0));
  meanshift->fit_and_predict(plhs,prhs+1);
}

//MEX_DEFINE(est_density)(int nlhs, mxArray* plhs[],
//      int nrhs, const mxArray* prhs[])
//{
//  InputArguments input(nrhs-2, prhs, 1);
//  OutputArguments output(nlhs, plhs, 1);
//  MeanshiftMexIface *meanshift = MexSession::get(input.get(0));
//  meanshift->est_density(plhs,prhs+1);
//}
//
//MEX_DEFINE(est_dist)(int nlhs, mxArray* plhs[],
//      int nrhs, const mxArray* prhs[])
//{
//  InputArguments input(nrhs-3, prhs, 1);
//  OutputArguments output(nlhs, plhs, 1);
//  MeanshiftMexIface *meanshift = MexSession::get(input.get(0));
//  meanshift->est_dist(plhs,prhs+1);
//}

MEX_DISPATCH

