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
template<typename Scalar,   
   typename Profile,
   template<typename S> class Chart>
class KernelAdapter<Scalar,Profile,Chart,BalloonVariableBandwidth> 
{
public:
  KernelAdapter(InputArguments input) : 
    kernel_(get_k(input)),
    min_support(get_min_support(input)),
    min_ratio(get_min_ratio(input))
  { }
private:
  int get_k(InputArguments input) 
  { return input.get<int>("bandwidth_parameter",20); }
  Scalar get_min_support(InputArguments input) 
  { return input.get<int>("min_support",10); }
  Scalar get_min_ratio(InputArguments input) 
  { return input.get<Scalar>("min_ratio",1); }
public:
  Kernel<Profile,Metric<Chart<Scalar> >,BalloonVariableBandwidth<Scalar> > kernel_;
  const int min_support;
  const double min_ratio;
};


template<typename Scalar,   
   typename Profile,
   template<typename S> class Chart>
class KernelAdapter<Scalar,Profile,Chart,SampleVariableBandwidth> 
{
public:
  KernelAdapter(InputArguments input) : 
    kernel_(get_h0(input),get_max_bw(input)),
    min_support(get_min_support(input)),
    min_ratio(get_min_ratio(input))
  { }
private:
  Scalar get_h0(InputArguments input) 
  { return input.get<Scalar>("bandwidth_parameter",0.2); }
  Scalar get_max_bw(InputArguments input) 
  { return input.get<Scalar>("max_bw",0.4); }
  Scalar get_min_support(InputArguments input) 
  { return input.get<int>("min_support",10); }
  Scalar get_min_ratio(InputArguments input) 
  { return input.get<Scalar>("min_ratio",1); }
public:
  Kernel<Profile,Metric<Chart<Scalar> >,SampleVariableBandwidth<Scalar> > kernel_;
  const int min_support;
  const double min_ratio;
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
    Scalar bw_sq = this->bw_sq(y);
    return Profile::k(Metric::sq_dist(x,y)/bw_sq)/bw_sq; 
  }  

  Scalar bw_sq(const EmbeddedVector& data_point)
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


template<typename Chart, typename Profile>
std::vector<typename Chart::Scalar> 
est_density(const std::vector<typename Chart::EmbeddedVector>& modes,
	    const std::vector<typename Chart::EmbeddedVector>& u,
	    const Kernel<Chart,Profile>& kernel)
{
  typedef typename Chart::Scalar Scalar;
  std::vector<Scalar> density;

  density.reserve(modes.size());

  kernel.init_data(u);
  for (auto modes_i : modes) 
    density.push_back(est_kernel_density(modes_i,u,kernel));
  
  return density;
}

template<typename Chart,typename Profile>
std::vector<typename Chart::Scalar> 
est_dist(const std::vector<typename Chart::EmbeddedVector>& u,
	 const std::vector<typename Chart::EmbeddedVector>& clust,
	 std::vector<unsigned int> lbl, Kernel<Chart,Profile>& kernel)
{
  typedef typename Chart::Scalar Scalar;
  typedef typename Chart::EmbeddedVector EmbeddedVector;

  std::vector<Scalar> dist;
  dist.reserve(u.size());
  for (const auto & zipped : boost::combine(u,lbl)) {
    int label;
    EmbeddedVector u_k;  
    boost::tie(u_k,label) = zipped;
    Scalar d = 0;
    if (label > 0)
      d = Metric<Chart>::dist_sq(clust[label-1],u_k);
    dist.push_back(d);
  }
  
  return dist;
}
