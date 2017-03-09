#ifndef MEANSHIFT_MEANSHIFT_HPP_
#define MEANSHIFT_MEANSHIFT_HPP_

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

#include "meanshift/profile.hpp"
#include "meanshift/chart.hpp"
#include "meanshift/kernel.hpp"

template <typename T, typename Compare>
std::vector<int> sort_permutation(
    std::vector<T> const& vec,
    Compare compare)
{
    std::vector<int> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](int i, int j){ return compare(vec[i], vec[j]); });
    return p;
}


template <typename T>
std::vector<T> apply_permutation(
    std::vector<T> const& vec,
    std::vector<int> const& p)
{
    std::vector<T> sorted_vec(p.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](int i){ return vec[i]; });
    return sorted_vec;
}



template<typename Kernel>
typename Kernel::Scalar
calc_likelihood(const typename Kernel::EmbeddedVector& u_0, 
		const std::vector<typename Kernel::EmbeddedVector>& v,
		Kernel& kernel)
{  
  typename Kernel::Scalar likelihood = 0.0;
  kernel.init_current(u_0);

  for (auto v_i : v)
    likelihood += kernel(u_0,v_i);

  likelihood /= v.size();

  return likelihood;
}

template<typename Kernel>
std::tuple<typename Kernel::EmbeddedVector,typename Kernel::Scalar>
calc_meanshift(const typename Kernel::EmbeddedVector& u0, 
	       const std::vector<typename Kernel::EmbeddedVector>& v, 
	       Kernel& kernel) 
{
  typedef typename Kernel::EmbeddedVector EmbeddedVector;
  typedef typename Kernel::TangentVector TangentVector;
  typedef typename Kernel::Scalar Scalar;
    
  Scalar renorm = 0;
  TangentVector cum_delta(0.0,0.0);
  
  Shadow<Kernel> shadow(kernel);
  kernel.init_current(u0);
  
  for (auto v_i : v) {
    TangentVector delta;
    Scalar w;
    std::tie(delta,w) = shadow.calc_weighted_tangent(u0,v_i);
    if (w > 0.0) {
      cum_delta += w*delta;
      renorm += w;
    }
  }

  EmbeddedVector u1 = u0;
  
  if (renorm > 0) {
    cum_delta /= renorm;			
    u1 = Kernel::Metric::expm(u0,cum_delta);
  }
  
  return std::make_tuple(u1,Kernel::Metric::sq_norm(cum_delta));
}


template<typename Kernel> 
std::tuple<std::vector<typename Kernel::EmbeddedVector> ,std::vector<int>>
clust_modes(std::vector<typename Kernel::EmbeddedVector>& modes, 
    std::vector<typename Kernel::Scalar>& likelihoods, 
    Kernel& kernel, int min_support, double min_ratio)
{
  typedef typename Kernel::Scalar Scalar;
  typedef typename Kernel::EmbeddedVector EmbeddedVector;
  typedef typename Kernel::TangentVector TangentVector;

  auto p = sort_permutation(likelihoods,
  [](Scalar const& a, Scalar const& b){ return a > b; });
  modes = apply_permutation(modes, p);

  std::vector<EmbeddedVector> clust;
  std::vector<int> lbl(modes.size());
  std::fill(lbl.begin(), lbl.end(), 0);

  std::vector<int> clust_sizes;
  if (modes.size() > 0) {
    clust.push_back(modes[0]);
    clust_sizes.push_back(1);
    lbl[0] = 1;
  }

  kernel.init_data(modes);
  for (int i = 1; i < modes.size(); i++)
  {
    kernel.init_current(modes[i]);
    for (int j = 0; j < clust.size(); j++)
    {
      // Scalar spacing = Kernel::Metric::sq_dist(clust[j],modes[i])/kernel.sq_bw(modes[o]);
      Scalar distance = Kernel::Metric::sq_dist2(clust[j],modes[i]);
      Scalar spacing2 = distance/kernel.h2_;
      if (spacing2 < min_ratio) 
      {
        lbl[i] = j + 1;
        clust_sizes[j]++;
        break;
      }
    }

    if (lbl[i] == 0) 
    {
      clust.push_back(modes[i]);
      lbl[i] = clust.size();
      clust_sizes.push_back(1);
    }
  }

  std::vector<EmbeddedVector> supported_clust;
  supported_clust.push_back(clust[0]);
  for (int i = 1; i < clust.size(); ++i)
  {
    if(clust_sizes[i] > min_support) 
    {
      supported_clust.push_back(clust[i]);
    } else {
      std::replace(lbl.begin(), lbl.end(), i+1, 0);
    }
  }
  auto rollback_permutation = sort_permutation(p,
  [](int const& a, int const& b){ return a < b; });
  lbl = apply_permutation(lbl, rollback_permutation);

  return std::make_tuple(supported_clust,lbl);
}


template<typename Kernel> 
std::tuple<std::vector<typename Kernel::EmbeddedVector> ,std::vector<typename Kernel::Scalar>>
est_modes(std::vector<typename Kernel::EmbeddedVector>& u, 
	  Kernel& kernel, 
	  const double rel_err = 1.0e-7, int max_iter = 200) 
{
  typedef typename Kernel::EmbeddedVector EmbeddedVector;
  typedef typename Kernel::TangentVector TangentVector;
  typedef typename Kernel::Scalar Scalar;

  std::vector<Scalar> likelihoods;
  std::vector<EmbeddedVector> modes;
  
  likelihoods.reserve(u.size());
  modes.reserve(u.size());			

  kernel.init_data(u);
  for (auto u_i : u) {
    int iter = 0;
    Scalar sq_dist = std::numeric_limits<Scalar>::max();
    do {
      std::tie(u_i,sq_dist) = calc_meanshift(u_i,u,kernel);
    }
    while (sq_dist > rel_err && iter++ < max_iter);
    likelihoods.push_back(calc_likelihood(u_i,u,kernel));
    modes.push_back(u_i);
  }
  return std::make_tuple(modes,likelihoods);
}

template<typename Kernel>
std::vector<typename Kernel::Scalar> 
est_density(const std::vector<typename Kernel::EmbeddedVector>& modes,
	    const std::vector<typename Kernel::EmbeddedVector>& u,
	    Kernel& kernel)
{
  std::vector<typename Kernel::Scalar> density;

  density.reserve(modes.size());

  kernel.init_data(u);
  for (auto modes_i : modes) 
    density.push_back(calc_likelihood(modes_i,u,kernel));
  
  return density;
}

template<typename Kernel>
std::vector<typename Kernel::Scalar> 
est_dist(const std::vector<typename Kernel::EmbeddedVector>& u,
      const std::vector<typename Kernel::EmbeddedVector>& clust,
      std::vector<int> lbl, Kernel& kernel)
{
  std::vector<typename Kernel::Scalar> dist;

  dist.reserve(u.size());
  typename Kernel::Scalar d;

  for (int i = 0; i < u.size(); ++i)
  {
    if (lbl[i] > 0) {
      d = Kernel::Metric::sq_dist2(clust[lbl[i]-1],u[i]);
    } else {
      d = 0;
    }
    dist.push_back(d);
  }
  return dist;
}

#endif
