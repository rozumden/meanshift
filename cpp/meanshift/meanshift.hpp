#ifndef __MEANSHIFT_MEANSHIFT_HPP__
#define __MEANSHIFT_MEANSHIFT_HPP__

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <tuple>

#include <boost/range/iterator_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range/combine.hpp>
#include <boost/range.hpp>

#include "meanshift/profile.hpp"
#include "meanshift/manifold/chart.hpp"
#include "meanshift/kernel.hpp"
#include "meanshift/metric.hpp"
#include "meanshift/bandwidth.hpp"

template<typename Kernel,typename Bandwidth>
typename kernel_traits<Kernel>::Scalar
est_kernel_density(const typename kernel_traits<Kernel>::Vector& x, 
		   const std::vector<typename kernel_traits<Kernel>::Vector>& X,
		   const Kernel& kernel,
		   Bandwidth& bandwidth)
{  
  typename kernel_traits<Kernel>::Scalar f_h = 0.0;
  typename kernel_traits<Kernel>::Scalar f_h_i = 0.0; 
    
  unsigned int i = 0;
    
  for (const auto& x_i : X) {
    const typename kernel_traits<Kernel>::Scalar h = bandwidth.predict(x_i,i++);    
    std::tie(f_h_i,std::ignore) = kernel.K(x,x_i,h);
    f_h += f_h_i;
  }
    
  f_h /= X.size();
    
  return f_h;
}

template<typename Kernel,typename Bandwidth>
std::tuple<typename kernel_traits<Kernel>::Vector,
	   typename kernel_traits<Kernel>::Scalar>
shift_point(const typename kernel_traits<Kernel>::Vector& x, 
	    const std::vector<typename kernel_traits<Kernel>::Vector>& X, 
	    Kernel& kernel,
	    Bandwidth& bandwidth) 
{
  typedef typename kernel_traits<Kernel>::Chart Chart;
  
  typedef typename chart_traits<Chart>::TangentVector TangentVector;
  typedef typename chart_traits<Chart>::Scalar Scalar;
    
  Scalar w = 0.0;
  Scalar w_i = 0.0;  
  
  TangentVector mu = TangentVector::Zero();
  TangentVector mu_i = TangentVector::Zero();
  
  unsigned int i = 0;
  for (const typename kernel_traits<Kernel>::Vector& x_i : X) {
    const typename chart_traits<Chart>::Scalar h = bandwidth.predict(x_i,i++);
    std::tie(w_i,mu_i) = kernel.G(x,x_i,h);
    w += w_i;
    mu += mu_i;
  }
  
  typename kernel_traits<Kernel>::Vector y(x);
  
  if (w > 0.0) {
    mu /= w;
    y = Chart::expm(x,mu);
  } 

  Scalar norm_mu = Metric<Chart>::norm(mu);  
  return std::make_tuple(y,norm_mu);
} 

template<typename Kernel, typename Bandwidth> 
std::tuple<std::vector<typename kernel_traits<Kernel>::Vector>, 
	   std::vector<typename kernel_traits<Kernel>::Scalar>,
	   std::vector<std::vector<typename kernel_traits<Kernel>::Vector> > >
est_modes(const std::vector<typename kernel_traits<Kernel>::Vector>& X_fit,
	  const std::vector<typename kernel_traits<Kernel>::Vector>& X_predict, 
	  Kernel& kernel,Bandwidth& bandwidth, 
	  const double rel_err = 1.0e-6,
	  int max_iter = 300) 
{
  typedef typename kernel_traits<Kernel>::Scalar Scalar;
  typedef typename kernel_traits<Kernel>::Vector Vector;
  
  std::vector<Scalar> likelihoods;
  std::vector<Vector> modes;
  std::vector<std::vector<Vector> > y_kt;
  
  y_kt.reserve(X_predict.size());
  likelihoods.reserve(X_predict.size());
  modes.reserve(X_predict.size());			

  for (const auto& x_i : X_predict) {
    y_kt.emplace_back();
    std::vector<Vector>& y_t = y_kt.back();
    y_t.reserve(max_iter);
    Scalar dist = std::numeric_limits<Scalar>::max();
    Vector y_i(x_i);
    unsigned int iter = 0;
    //    while (!peturbed) {
    do {
      std::tie(y_i,dist) = shift_point(y_i,X_fit,kernel,bandwidth);
      y_t.push_back(y_i);
    }
    while (dist > rel_err && iter++ < max_iter);
    //      y_i = peturb(y_i);
    //}
    Scalar likelihood = est_kernel_density(y_i,X_fit,kernel,bandwidth);
    modes.push_back(y_i);
    likelihoods.push_back(likelihood);
  }
  return std::make_tuple(modes,likelihoods,y_kt);
}

template<typename Kernel, typename Bandwidth> 
std::tuple<std::vector<typename kernel_traits<Kernel>::Vector>, 
	   std::vector<unsigned int> >
clust_modes(const std::vector<typename kernel_traits<Kernel>::Vector>& X1,
	    const std::vector<typename kernel_traits<Kernel>::Vector>& modes0, 
	    std::vector<typename kernel_traits<Kernel>::Scalar>& likelihoods, 
	    Kernel& kernel,Bandwidth& bandwidth,const unsigned int min_support) 
{
  using namespace std::placeholders;

  typedef typename kernel_traits<Kernel>::Scalar Scalar;
  typedef typename kernel_traits<Kernel>::Vector Vector;
  typedef typename kernel_traits<Kernel>::Chart Chart;

  const unsigned int N = likelihoods.size();

  std::vector<unsigned int> idx(N);
  std::iota(begin(idx), end(idx), static_cast<unsigned int>(0));
  std::sort(begin(idx), end(idx),
	    [&](const unsigned int i, const unsigned int j) 
	    { 
	      return likelihoods[i] > likelihoods[j]; 
	    });

  std::vector<typename kernel_traits<Kernel>::Vector> modes(N);
  std::transform(begin(idx), end(idx), begin(modes),
		 [&](const unsigned int i) 
		 { 
		   return modes0[i]; 
		 });
  
  std::vector<unsigned int> lbl;
  std::vector<Vector> clust;
  std::vector<unsigned int> clust_sizes;
  std::vector<Scalar> dd;
  
  lbl.reserve(N);
  clust.reserve(N);
  clust_sizes.reserve(N);
  dd.reserve(N);

  if (modes.size() > 0) {
    clust.push_back(modes[0]);
    clust_sizes.push_back(1);
    lbl.push_back(1);

    for (int i = 1; i < modes.size(); i++) {
      std::transform(begin(clust),end(clust),std::back_inserter(dd),
		     std::bind(&Metric<Chart>::dist,_1,modes[i]));
      typename std::vector<Scalar>::const_iterator iter = 
	std::min_element(begin(dd),end(dd));
      const Scalar h = bandwidth.predict(modes[i],idx[i]);
      const Scalar d = *iter;
      if (d <= h) {
	unsigned int j = 
	  std::distance(typename std::vector<Scalar>::const_iterator(begin(dd)),iter);
	lbl.push_back(j+1);
	clust_sizes[j]++;
      } else {
	clust.push_back(modes[i]);
	lbl.push_back(clust.size());
	clust_sizes.push_back(1);
      }
    
      dd.clear();
    }
  }
  
  std::vector<Vector> supported_clust;
  supported_clust.reserve(N);

  for (int i = 0; i < clust.size(); ++i) {
    if (clust_sizes[i] >= min_support) 
      supported_clust.push_back(clust[i]);
    else 
      std::replace(lbl.begin(), lbl.end(), i+1, 0);
  }
  
  std::vector<unsigned int> idx2(N);
  std::iota(begin(idx2), end(idx2), static_cast<unsigned int>(0));
  std::sort(begin(idx2), end(idx2),
	    [&](const unsigned int i, const unsigned int j) 
	    { 
	      return idx[i] < idx[j]; 
	    });
  
  std::vector<unsigned int> labels(N,0);
  std::transform(begin(idx2), end(idx2), begin(labels),
		 [&](int i)
		 { return lbl[i]; });

  return std::make_tuple(supported_clust,labels);
}

template<typename EmbeddedVector,typename Kernel,typename Bandwidth>
std::tuple<std::vector<EmbeddedVector>,
	   std::vector<unsigned int>,
	   std::vector<EmbeddedVector>, 
	   std::vector<typename EmbeddedVector::Scalar>,
	   std::vector<std::vector<EmbeddedVector> > >
predict(const std::vector<EmbeddedVector>& X_fit,
	const std::vector<EmbeddedVector>& X_predict,
	Kernel kernel,Bandwidth& bandwidth,
	const unsigned int min_support) 
{
  std::vector<EmbeddedVector> modes,clust;
  std::vector<typename EmbeddedVector::Scalar> likelihoods;
  std::vector<unsigned int> labels;
  std::vector<std::vector<EmbeddedVector> > tracks;
  
  std::tie(modes,likelihoods,tracks) = est_modes(X_fit,X_predict,kernel,bandwidth);
  std::tie(clust,labels) = 
    clust_modes(X_fit,modes,likelihoods,kernel,bandwidth,min_support); 
  
  return std::make_tuple(clust,labels,modes,likelihoods,tracks); 
}


#endif
