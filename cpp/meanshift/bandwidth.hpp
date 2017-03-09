#ifndef __MEANSHIFT_BANDWIDTH_HPP__
#define __MEANSHIFT_BANDWIDTH_HPP__

#include "common/std_ext.hpp"

template<typename Scalar>
struct FixedBandwidth
{
  FixedBandwidth(const Scalar h) : h_(h)
  { }
  
  template<typename Vector>
  Scalar predict(const Vector& x_i,const size_t i)
  { 
    return h_;
  }
  
  Scalar h_;
};

template<typename Scalar>
struct SampleVariableBandwidth
{
  SampleVariableBandwidth() 
  { }

  SampleVariableBandwidth(const std::vector<Scalar>& h)
  {
    std::copy(begin(h),end(h),std::back_inserter(h_));
  }

  void set_bandwidth(const std::vector<Scalar>& h)
  {
    h_.clear();
    std::copy(begin(h),end(h),std::back_inserter(h_));
  }

  template<typename Vector>
  Scalar predict(const Vector& x_i,const size_t i)
  { 
    return h_[i];
  }
  
  std::vector<Scalar> h_;
};


template <typename EmbeddedVector, 
	  typename Metric>
std::vector<typename EmbeddedVector::Scalar>
est_pilot_bandwidths(const std::vector<EmbeddedVector>& X,
		     const unsigned int kth_nn,
		     Metric metric) 
{
  using namespace std::placeholders;
    
  typedef typename EmbeddedVector::Scalar Scalar;

  std::vector<Scalar> d2;
  std::vector<Scalar> h0,h1;

  d2.reserve(X.size());
  h0.reserve(X.size());

  for (const auto& x_i : X) {
    std::transform(begin(X),end(X),std::back_inserter(d2),
		   std::bind(&Metric::dist_sq,x_i,_1));
    std::nth_element(begin(d2),begin(d2)+kth_nn,end(d2),std::less<Scalar>());
    const Scalar h_knn = d2[kth_nn]+10*std::numeric_limits<Scalar>::epsilon();
    h0.push_back(sqrt(h_knn));
    d2.clear();
  } 
  
  std::copy(begin(h0),end(h0),back_inserter(h1));
  unsigned int nth = static_cast<unsigned int>(ceil(0.3*X.size()))-1;
  std::nth_element(begin(h1),begin(h1)+nth,end(h1));
  std::replace_if(begin(h0),end(h0),
		  std::bind(std::greater<Scalar>(),_1,h1[nth]),h1[nth]);

  return h0;
}

template <typename EmbeddedVector, 
	  typename Kernel,
	  typename Bandwidth>
std::vector<typename EmbeddedVector::Scalar>
est_sample_bandwidths(const std::vector<EmbeddedVector>& X,
		      Bandwidth& pilot_bandwidth,
		      Kernel kernel) 
{
  typename EmbeddedVector::Scalar lglambda = 
    static_cast<typename EmbeddedVector::Scalar>(0.0);
  std::vector<typename EmbeddedVector::Scalar> f;
  f.reserve(X.size());
  for (const auto& x_i : X) {
    f.push_back(est_kernel_density(x_i,X,kernel,pilot_bandwidth));
    lglambda += log(f.back());
  }
  lglambda /= X.size();
  
  typename EmbeddedVector::Scalar lambda = exp(lglambda);
  std::vector<typename EmbeddedVector::Scalar> hh;
  hh.reserve(X.size());
  unsigned int i = 0;
  for (const auto& x_i : X) {
    typename EmbeddedVector::Scalar h0 = pilot_bandwidth.predict(x_i,i);
    typename EmbeddedVector::Scalar h = h0*sqrt(lambda/f[i]);
    hh.push_back(h);
    i++;
  }

  return hh;
}

#endif
