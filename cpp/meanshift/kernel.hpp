#ifndef MEANSHIFT_KERNEL_HPP_
#define MEANSHIFT_KERNEL_HPP_

#include <meanshift/metric.hpp>
#include <unordered_map>
#include <functional>
#include <math.h>

template <typename Kernel>
class Shadow;

template<typename Profile, typename Metric, typename Bandwidth>
class Kernel;


template<typename Scalar>
class FixedBandwidth
{  
public:	
  FixedBandwidth(const Scalar& h) : h2_(h*h) { } 
public:
  Scalar h2_;
};

template<typename Scalar>
class BalloonVariableBandwidth
{  
public:  
  BalloonVariableBandwidth(const int knn) : k(knn), h2_(0) { }
public:
  Scalar h2_;
  int k;
};

template<typename Scalar>
class SampleVariableBandwidth
{  
public:  
  SampleVariableBandwidth(const Scalar h0, const Scalar m) : 
   h0_(h0), h2_(h0*h0), max_bw(m) { }
public:
  Scalar h0_;
  Scalar h2_;
  Scalar max_bw;
};

template<typename _Profile, typename _Metric>
class Kernel<_Profile,_Metric,SampleVariableBandwidth<typename _Metric::Scalar> > : public SampleVariableBandwidth<typename _Metric::Scalar>
{
  template <typename T>
  friend class Shadow;

public: 
  typedef _Profile Profile;
  typedef _Metric Metric;
  typedef typename Metric::EmbeddedVector EmbeddedVector;
  typedef typename Metric::TangentVector TangentVector;
  typedef typename Metric::Scalar Scalar;


public:
  Kernel(Scalar h0, const Scalar max_bw) : 
    lambda(1),SampleVariableBandwidth<typename Metric::Scalar>(h0,max_bw) 
  { }
    
  Scalar operator()(const EmbeddedVector& x, const EmbeddedVector& y)
  {
    Scalar sq_bw = this->sq_bw(y);
    return Profile::k(Metric::sq_dist(x,y)/sq_bw)/sq_bw; 
  }  

  Scalar calc_likelihood(const EmbeddedVector& u_0, 
    const std::vector<EmbeddedVector>& v)
  {  
    Scalar likelihood = 0.0;

    for (auto v_i : v)
      likelihood += (*this)(u_0,v_i);

    likelihood /= v.size();

    return likelihood;
  }

  Scalar sq_bw(const EmbeddedVector& data_point)
  {
    return this->map[data_point];
  }

  void init_current(const EmbeddedVector& current) {}

  void init_data(const std::vector<EmbeddedVector>& data)
  {
    for (auto d : data)
    {
      map[d] = this->h0_*this->h0_;
    }
    std::vector<Scalar> l;
    this->lambda = 1;
    for (auto d : data) 
    {
      this->lambda *= pow(this->calc_likelihood(d,data),1.0/data.size());
      l.push_back(this->calc_likelihood(d,data));
    }
    for (int i = 0; i < data.size(); ++i)
    {
      Scalar bw = (this->h0_)*sqrt(this->lambda /l[i]);
      bw = std::min(bw,this->max_bw);
      map[data[i]] = bw*bw;
    }
  }
private:
  Scalar lambda;
  std::unordered_map<EmbeddedVector, Scalar> map;
};

template<typename _Profile, typename _Metric>
class Kernel<_Profile,_Metric,BalloonVariableBandwidth<typename _Metric::Scalar> > : public BalloonVariableBandwidth<typename _Metric::Scalar>
{
  template <typename T>
  friend class Shadow;

public: 
  typedef _Profile Profile;
  typedef _Metric Metric;
  typedef typename Metric::EmbeddedVector EmbeddedVector;
  typedef typename Metric::TangentVector TangentVector;
  typedef typename Metric::Scalar Scalar;

public:
  Kernel(int k) : 
    BalloonVariableBandwidth<typename Metric::Scalar>(k) 
  { }
    
  Scalar operator()(const EmbeddedVector& x, const EmbeddedVector& y)
  {
    Scalar sq_bw = this->sq_bw(y);
    return Profile::k(Metric::sq_dist(x,y)/sq_bw)/sq_bw; 
  }  

  Scalar sq_bw(const EmbeddedVector& data_point)
  {
    return this->h2_;
  }

  void init_current(const EmbeddedVector& current)
  {
    Scalar min_bw = 0.1, max_bw = 0.5;
    std::vector<Scalar> dist;
    for (int i = 0; i < data_->size(); ++i)
      dist.push_back(Metric::sq_dist(current,data_->at(i)));
    int n = std::min((size_t)(this->k - 1), data_->size()-1);
    std::nth_element(dist.begin(), dist.begin()+n, dist.end());
    this->h2_ = std::min(std::max(dist[n],min_bw*min_bw),max_bw*max_bw);
  }

  void init_data(const std::vector<EmbeddedVector>& data) {
    data_ = &data;
  }
private:
  const std::vector<EmbeddedVector>* data_;
};

template<typename _Profile, typename _Metric>
class Kernel<_Profile,_Metric,FixedBandwidth<typename _Metric::Scalar> > : public FixedBandwidth<typename _Metric::Scalar>
{
  template <typename T>
  friend class Shadow;

public: 
  typedef _Profile Profile;
  typedef _Metric Metric;
  typedef typename Metric::EmbeddedVector EmbeddedVector;
  typedef typename Metric::TangentVector TangentVector;
  typedef typename Metric::Scalar Scalar;

public:
  Kernel(const Scalar& h) : 
    FixedBandwidth<typename Metric::Scalar>(h)
  {  }
    
  Scalar operator()(const EmbeddedVector& x, const EmbeddedVector& y)
  {
    Scalar sq_bw = this->sq_bw(y);
    return Profile::k(Metric::sq_dist(x,y)/sq_bw)/sq_bw; 
  }  

  Scalar sq_bw(const EmbeddedVector& data_point)
  {
    return this->h2_;
  }

  void init_current(const EmbeddedVector& current) {}
  void init_data(const std::vector<EmbeddedVector>& data) {}
};

template<typename Kernel>
class Shadow
{
public:
  typedef typename Kernel::Profile Profile;
  typedef typename Kernel::EmbeddedVector EmbeddedVector;
  typedef typename Kernel::TangentVector TangentVector;
  typedef typename Kernel::Scalar Scalar;
  typedef typename Kernel::Metric Metric;

  Shadow(const Kernel& kernel) : kernel_(kernel)
  {  }
    
  Scalar operator()(const EmbeddedVector& x, const EmbeddedVector& y) 
  {
    Scalar sq_bw = kernel_.sq_bw(y);
    return Profile::g(Metric::sq_dist(x,y)/sq_bw)/sq_bw/sq_bw; 
  }
    
    
  std::tuple<TangentVector,Scalar> calc_weighted_tangent(const EmbeddedVector& x, const EmbeddedVector& y) 
  {
    TangentVector u = Metric::Chart::logm(x,y);
    Scalar sq_bw = kernel_.sq_bw(y);
    Scalar weight = Profile::g(Metric::sq_norm(u)/sq_bw)/sq_bw/sq_bw;
    return std::make_tuple(u,weight);
  }
    
private:
  Kernel kernel_;
};



#endif
