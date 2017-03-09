#ifndef MEANSHIFT_MEX_FACTORY_HPP_
#define MEANSHIFT_MEX_FACTORY_HPP_

#include <string>
#include <unordered_map>

#include <boost/function.hpp>
#include <boost/functional/factory.hpp>
#include <boost/container/flat_map.hpp>

#include "meanshift/meanshift.hpp"

#include <mexplus.h>

using namespace mexplus;

class MeanshiftMexIface
{
 public:
  virtual void est_density(mxArray *plhs[], const mxArray *prhs[]) = 0;
  // virtual void est_modes(mxArray *plhs[], const mxArray *prhs[]) = 0;
  // virtual void clust_modes(mxArray *plhs[], const mxArray *prhs[]) = 0;
  virtual void train(mxArray *plhs[], const mxArray *prhs[]) = 0;
  virtual void est_dist(mxArray *plhs[], const mxArray *prhs[]) = 0;

  virtual ~MeanshiftMexIface() {};
};

template <typename Scalar, 
	  typename Profile,
	  template<typename S> class Chart,
	  template<typename S> class Bandwidth>
class KernelAdapter;


template<typename Scalar,	  
	 typename Profile,
	 template<typename S> class Chart>
class KernelAdapter<Scalar,Profile,Chart,FixedBandwidth> 
{
public:
  KernelAdapter(InputArguments input) : 
    kernel_(get_bandwidth(input)),
    min_support(get_min_support(input)),
    min_ratio(get_min_ratio(input))
  { }
private:
  Scalar get_bandwidth(InputArguments input) 
  { return input.get<Scalar>("bandwidth_parameter",0.3); }
  Scalar get_min_support(InputArguments input) 
  { return input.get<int>("min_support",10); }
  Scalar get_min_ratio(InputArguments input) 
  { return input.get<Scalar>("min_ratio",1); }
public:
  Kernel<Profile,Metric<Chart<Scalar> >,FixedBandwidth<Scalar> > kernel_;
  const int min_support;
  const double min_ratio;
};


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
 

template <typename Scalar, typename Profile, 
	  template<typename S> class Chart,
	  template<typename S> class Bandwidth>
class MeanshiftMex : public MeanshiftMexIface,
		     private KernelAdapter<Scalar,Profile,Chart,Bandwidth>		     
{
public:
  typedef typename Chart<Scalar>::EmbeddedVector EmbeddedVector;
  
  MeanshiftMex(InputArguments input) : 
    KernelAdapter<Scalar,Profile,Chart,Bandwidth>(input) { }

  virtual void train(mxArray *plhs[], const mxArray *prhs[]) override
  { 
    std::vector<EmbeddedVector> u;
    mexplus::MxArray::to<std::vector<EmbeddedVector> >(prhs[0],&u);
    int min_support = 10;
    double min_ratio = 1;

    std::vector<EmbeddedVector> modes;
    std::vector<Scalar> likelihoods;
    std::tie(modes,likelihoods) = ::est_modes(u,this->kernel_);

    std::vector<EmbeddedVector> clust;
    std::vector<int> lbl;
    std::tie(clust,lbl) = ::clust_modes(modes,likelihoods,
       this->kernel_,this->min_support,this->min_ratio); 
    if (clust.size() > 0) {
      plhs[0] = mexplus::MxArray::from(clust);
    } else {
      std::vector<int> dummy;
      plhs[0] = mexplus::MxArray::from(dummy);
    }
    plhs[1] = mexplus::MxArray::from(lbl);
    plhs[2] = mexplus::MxArray::from(modes);
    plhs[3] = mexplus::MxArray::from(likelihoods);
  }

  virtual void est_dist(mxArray *plhs[], const mxArray *prhs[]) override
  { 
    std::vector<EmbeddedVector> u;
    mexplus::MxArray::to<std::vector<EmbeddedVector> >(prhs[0],&u);

    std::vector<EmbeddedVector> clust;
    mexplus::MxArray::to<std::vector<EmbeddedVector> >(prhs[1],&clust);

    std::vector<int> lbl;
    mexplus::MxArray::to<std::vector<int> >(prhs[2],&lbl);

    std::vector<Scalar> dist = ::est_dist(u,clust,lbl,this->kernel_);

    plhs[0] = mexplus::MxArray::from(dist);
  }

  // virtual void est_modes(mxArray *plhs[], const mxArray *prhs[]) override
  // { 
  //   std::vector<EmbeddedVector> u;
  //   mexplus::MxArray::to<std::vector<EmbeddedVector> >(prhs[0],&u);
  //   std::vector<EmbeddedVector> modes;
  //   std::vector<Scalar> likelihoods;
  //   std::tie(modes,likelihoods) = ::est_modes(u,this->kernel_); 
  //   plhs[0] = mexplus::MxArray::from(modes);
  //   plhs[1] = mexplus::MxArray::from(likelihoods);
  // }

  // virtual void clust_modes(mxArray *plhs[], const mxArray *prhs[]) override
  // { 
  //   std::vector<EmbeddedVector> u;
  //   std::vector<EmbeddedVector> modes;
  //   std::vector<Scalar> likelihoods;
  //   mexplus::MxArray::to<std::vector<EmbeddedVector> >(prhs[0],&modes);
  //   mexplus::MxArray::to<std::vector<Scalar> >(prhs[1],&likelihoods);
  //   std::vector<EmbeddedVector> clust;
  //   std::vector<int> lbl;
  //   std::tie(clust,lbl) = ::clust_modes(modes,likelihoods,this->kernel_); 
  //   if (clust.size() > 0) {
  //     plhs[0] = mexplus::MxArray::from(clust);
  //   } else {
  //     std::vector<int> dummy;
  //     plhs[0] = mexplus::MxArray::from(dummy);
  //   }
  //   plhs[1] = mexplus::MxArray::from(lbl);
  // }

  virtual void est_density(mxArray *plhs[], const mxArray *prhs[]) override
  {
    typedef Kernel<Profile,Metric<Chart<Scalar> >,Bandwidth<Scalar> > Kernel_;

    std::vector<typename Kernel_::EmbeddedVector> modes,u;

    mexplus::MxArray::to<std::vector<typename Kernel_::EmbeddedVector> >(prhs[0],&modes);
    mexplus::MxArray::to<std::vector<typename Kernel_::EmbeddedVector> >(prhs[1],&u);
    std::vector<Scalar> density = ::est_density(modes,u,this->kernel_);
    plhs[0] = mexplus::MxArray::from(density);
  }

};

typedef MeanshiftMex<double,Gaussian,S2,FixedBandwidth> DoubleGaussianS2Fixed;
typedef MeanshiftMex<double,Epanechnikov,S2,FixedBandwidth> DoubleEpanechnikovS2Fixed;
typedef MeanshiftMex<float,Gaussian,S2,FixedBandwidth> SingleGaussianS2Fixed;
typedef MeanshiftMex<float,Epanechnikov,S2,FixedBandwidth> SingleEpanechnikovS2Fixed;

typedef MeanshiftMex<double,Gaussian,S2,BalloonVariableBandwidth> DoubleGaussianS2Balloon;
typedef MeanshiftMex<double,Epanechnikov,S2,BalloonVariableBandwidth> DoubleEpanechnikovS2Balloon;
typedef MeanshiftMex<float,Gaussian,S2,BalloonVariableBandwidth> SingleGaussianS2Balloon;
typedef MeanshiftMex<float,Epanechnikov,S2,BalloonVariableBandwidth> SingleEpanechnikovS2Balloon;

typedef MeanshiftMex<double,Gaussian,S2,SampleVariableBandwidth> DoubleGaussianS2Sample;
typedef MeanshiftMex<double,Epanechnikov,S2,SampleVariableBandwidth> DoubleEpanechnikovS2Sample;
typedef MeanshiftMex<float,Gaussian,S2,SampleVariableBandwidth> SingleGaussianS2Sample;
typedef MeanshiftMex<float,Epanechnikov,S2,SampleVariableBandwidth> SingleEpanechnikovS2Sample;

typedef std::tuple<std::string, std::string, std::string, std::string> FactoryKey;

struct FactoryHash : public std::unary_function<FactoryKey,std::size_t>
{
 std::size_t operator()(const FactoryKey& key) const
 {
return std::hash<std::string>()(std::get<0>(key)) ^  
       std::hash<std::string>()(std::get<1>(key)) ^
       std::hash<std::string>()(std::get<2>(key)) ^
       std::hash<std::string>()(std::get<3>(key));
}
};

typedef boost::function<MeanshiftMexIface*(InputArguments input)> FactoryProduct;

std::unordered_map<FactoryKey,FactoryProduct,FactoryHash> factory = 
  { { std::make_tuple("double","epanechnikov","s2","fixed"), std::bind(boost::factory<DoubleEpanechnikovS2Fixed*>(),std::placeholders::_1) },
    { std::make_tuple("double","gaussian","s2","fixed"), std::bind(boost::factory<DoubleGaussianS2Fixed*>(),std::placeholders::_1) },
    { std::make_tuple("single","epanechnikov","s2","fixed"), std::bind(boost::factory<SingleEpanechnikovS2Fixed*>(),std::placeholders::_1) },
    { std::make_tuple("single","gaussian","s2","fixed"), std::bind(boost::factory<SingleGaussianS2Fixed*>(),std::placeholders::_1) } ,

    { std::make_tuple("double","epanechnikov","s2","balloon"), std::bind(boost::factory<DoubleEpanechnikovS2Balloon*>(),std::placeholders::_1) },
    { std::make_tuple("double","gaussian","s2","balloon"), std::bind(boost::factory<DoubleGaussianS2Balloon*>(),std::placeholders::_1) },
    { std::make_tuple("single","epanechnikov","s2","balloon"), std::bind(boost::factory<SingleEpanechnikovS2Balloon*>(),std::placeholders::_1) },
    { std::make_tuple("single","gaussian","s2","balloon"), std::bind(boost::factory<SingleGaussianS2Balloon*>(),std::placeholders::_1) },

    { std::make_tuple("double","epanechnikov","s2","sample"), std::bind(boost::factory<DoubleEpanechnikovS2Sample*>(),std::placeholders::_1) },
    { std::make_tuple("double","gaussian","s2","sample"), std::bind(boost::factory<DoubleGaussianS2Sample*>(),std::placeholders::_1) },
    { std::make_tuple("single","epanechnikov","s2","sample"), std::bind(boost::factory<SingleEpanechnikovS2Sample*>(),std::placeholders::_1) },
    { std::make_tuple("single","gaussian","s2","sample"), std::bind(boost::factory<SingleGaussianS2Sample*>(),std::placeholders::_1) }  }; 


#endif
