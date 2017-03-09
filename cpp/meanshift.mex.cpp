#include <stdio.h>
#include <mex.h>
#include <array>
#include <algorithm>
#include <math.h>

#include <mexplus.h>

#include "meanshift/mexplus_ext.hpp"
#include "meanshift/std_ext.hpp"

#include "meanshift/meanshift.hpp"
#include "meanshift/mex_factory.hpp"

using namespace mexplus;

typedef mexplus::Session<MeanshiftMexIface> MexSession;

MEX_DEFINE(new) (int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) 
{
  InputArguments input(nrhs, prhs, 0, 8, 
    "profile","manifold","bandwidth",
    "bandwidth_parameter","scalar",
    "min_support","min_ratio","max_bw");
  OutputArguments output(nlhs, plhs, 1);
  
  std::string profile = input.get<std::string>("profile","epanechnikov");
  std::transform(profile.begin(),profile.end(),profile.begin(),::tolower);

  std::string manifold = input.get<std::string>("manifold","s2");
  std::transform(manifold.begin(),manifold.end(),manifold.begin(),::tolower);

  std::string bandwidth = input.get<std::string>("bandwidth","fixed");
  std::transform(bandwidth.begin(),bandwidth.end(),bandwidth.begin(),::tolower);

  std::string scalar = input.get<std::string>("scalar","double");
  std::transform(scalar.begin(),scalar.end(),scalar.begin(),::tolower);

  output.set(0,
	     MexSession::create(factory[std::make_tuple(scalar,
        profile,manifold,bandwidth)](input)));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[],
		    int nrhs, const mxArray* prhs[]) 
{
  InputArguments input(nrhs, prhs, 1);
  OutputArguments output(nlhs, plhs, 0);
  MexSession::destroy(input.get(0));
}

MEX_DEFINE(train)(int nlhs, mxArray* plhs[],
          int nrhs, const mxArray* prhs[])
{
  InputArguments input(nrhs-1, prhs, 1);
  OutputArguments output(nlhs, plhs, 4);
  MeanshiftMexIface *meanshift = MexSession::get(input.get(0));
  meanshift->train(plhs,prhs+1);
}

MEX_DEFINE(est_density)(int nlhs, mxArray* plhs[],
      int nrhs, const mxArray* prhs[])
{
  InputArguments input(nrhs-2, prhs, 1);
  OutputArguments output(nlhs, plhs, 1);
  MeanshiftMexIface *meanshift = MexSession::get(input.get(0));
  meanshift->est_density(plhs,prhs+1);
}

MEX_DEFINE(est_dist)(int nlhs, mxArray* plhs[],
      int nrhs, const mxArray* prhs[])
{
  InputArguments input(nrhs-3, prhs, 1);
  OutputArguments output(nlhs, plhs, 1);
  MeanshiftMexIface *meanshift = MexSession::get(input.get(0));
  meanshift->est_dist(plhs,prhs+1);
}


// MEX_DEFINE(est_modes)(int nlhs, mxArray* plhs[],
// 		      int nrhs, const mxArray* prhs[])
// {
//   InputArguments input(nrhs-1, prhs, 1);
//   OutputArguments output(nlhs, plhs, 2);
//   MeanshiftMexIface *meanshift = MexSession::get(input.get(0));
//   meanshift->est_modes(plhs,prhs+1);
// }

// MEX_DEFINE(clust_modes)(int nlhs, mxArray* plhs[],
//           int nrhs, const mxArray* prhs[])
// {
//   InputArguments input(nrhs-2, prhs, 1);
//   OutputArguments output(nlhs, plhs, 2);
//   MeanshiftMexIface *meanshift = MexSession::get(input.get(0));
//   meanshift->clust_modes(plhs,prhs+1);
// }


MEX_DISPATCH
