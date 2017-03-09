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
#ifndef _MEANSHIFT_MEX_FACTORY_HPP_
#define _MEANSHIFT_MEX_FACTORY_HPP_

#include <string>
#include <unordered_map>

#include <boost/function.hpp>

#include <boost/container/flat_map.hpp>
#include <boost/functional/factory.hpp>
#include <boost/functional/value_factory.hpp>

#include "mexplus/arguments.h"
#include "mexplus/dispatch.h"

#include "mexplusplus/eigen/vector.hpp"
#include "meanshift/manifold/vector_space.hpp"

#include "meanshift/meanshift.hpp"
#include "meanshift/bandwidth.hpp"
#include "meanshift/metric.hpp"

using namespace mexplus;

struct MeanshiftMexIface
{
public:
  virtual void fit(mxArray *plhs[], const mxArray *prhs[]) = 0;
  virtual void predict(mxArray *plhs[], const mxArray *prhs[]) = 0;
  virtual void fit_and_predict(mxArray *plhs[], const mxArray *prhs[]) = 0;
  
  virtual ~MeanshiftMexIface() {};
};

template <typename T>
struct BandwidthMex;
  
template <typename Scalar>
struct BandwidthMex<FixedBandwidth<Scalar> >
{
  BandwidthMex(const InputArguments& input) : 
    bandwidth_(input.get<Scalar>("bandwidth",1.0))
  { }

  template <typename EmbeddedVector, typename Metric, typename Kernel>
  void init(const std::vector<EmbeddedVector>& X, Metric metric, Kernel kernel)
  { 
    N_ = X.size();
  }
 
  std::vector<Scalar> get_bandwidths()
  {
    std::vector<Scalar> bandwidths(N_);
    std::fill(begin(bandwidths),end(bandwidths),bandwidth_.h_);
    return bandwidths;
  }
 
  FixedBandwidth<Scalar> bandwidth_;
  Scalar h_;
  unsigned int N_;
};

template <typename Scalar>
struct BandwidthMex<SampleVariableBandwidth<Scalar> >
{
  BandwidthMex(const InputArguments& input) : 
    pilot_knn_(input.get<unsigned int>("pilot_knn",5))
  { }

  template <typename EmbeddedVector, typename Metric, typename Kernel>
  void init(const std::vector<EmbeddedVector>& X, 
	    Metric metric, 
	    Kernel kernel)
  {
    std::vector<Scalar> h0 = est_pilot_bandwidths(X,pilot_knn_,metric);
    SampleVariableBandwidth<Scalar> pilot_bandwidth(h0);
    //    std::vector<Scalar> h = est_sample_bandwidths(X,pilot_bandwidth,kernel);
    bandwidth_.set_bandwidth(h0);
    //bandwidth_.set_bandwidth(h);
  }

  std::vector<Scalar> get_bandwidths()
  {
    return bandwidth_.h_;
  }

  SampleVariableBandwidth<Scalar> bandwidth_;
  unsigned int pilot_knn_;
};

template<typename Scalar,	  
	 template<typename T> class Chart,
	 typename Shape,
	 template <typename T> class Bandwidth>
class MeanshiftMex : public MeanshiftMexIface,
		     public BandwidthMex<Bandwidth<Scalar> >
{
private:
  unsigned int get_min_support(const InputArguments& input) 
  { 
    return input.get<unsigned int>("min_support",10); 
  }
  
  std::string get_bandwidth_type(const InputArguments& input) 
  { 
    return input.get<std::string>("bandwidth_type","fixed"); 
  }
  
  unsigned int get_pilot_knn(const InputArguments& input) 
  { 
    return input.get<unsigned int>("pilot_knn",5); 
  }
  
public:
  typedef typename chart_traits<Chart<Scalar> >::EmbeddedVector EmbeddedVector;
  
  MeanshiftMex(const InputArguments& input) :   
    BandwidthMex<Bandwidth<Scalar> >(input),
    min_support_(this->get_min_support(input)),
    bandwidth_type_(this->get_bandwidth_type(input))
  {  
    
  }

  virtual void fit(mxArray *plhs[], const mxArray *prhs[]) override 
  {
    X_fit_ = 
      std::shared_ptr<std::vector<EmbeddedVector> >(new std::vector<EmbeddedVector>());     mexplus::MxArray::to<std::vector<EmbeddedVector> >(prhs[0],X_fit_.get());
    BandwidthMex<Bandwidth<Scalar> >::init(*X_fit_,metric_,kernel_);
  }

  virtual void predict(mxArray *plhs[], const mxArray *prhs[]) override
  { 
    if (X_fit_ == nullptr) 
      MEXPLUS_ERROR("Call fit before predict!");
    
    std::vector<EmbeddedVector> X_predict;
    mexplus::MxArray::to<std::vector<EmbeddedVector> >(prhs[0],&X_predict);
    
    std::vector<EmbeddedVector> clust,modes;
    std::vector<std::vector<EmbeddedVector> > y_k;
    std::vector<unsigned int> labels;
    std::vector<Scalar> likelihoods;
    
    std::tie(clust,labels,modes,likelihoods,y_k) = 
      ::predict(*X_fit_,X_predict, kernel_,this->bandwidth_,min_support_);
    
    MxArray cell_array(MxArray::Cell(1, X_predict.size()));
    unsigned int i = 0;
    for (const std::vector<EmbeddedVector>& y_i : y_k)
      cell_array.set(i++,y_i);
    
    plhs[0] = mexplus::MxArray::from(clust);
    plhs[1] = mexplus::MxArray::from(labels);
    plhs[2] = mexplus::MxArray::from(modes);
    plhs[3] = mexplus::MxArray::from(likelihoods);
    plhs[4] = mexplus::MxArray::from(this->get_bandwidths());
    plhs[5] = cell_array.release();
  }
  
  virtual void fit_and_predict(mxArray *plhs[], const mxArray *prhs[]) override
  { 
    this->fit(plhs,prhs);
    this->predict(plhs,prhs);
  }

  const unsigned int min_support_;
  std::string bandwidth_type_;
  
  Kernel<Chart<Scalar>, Profile<Shape> > kernel_;
  Metric<Chart<Scalar> > metric_;

  std::shared_ptr<std::vector<EmbeddedVector> > X_fit_;
};

typedef MeanshiftMex<float,R2,Gaussian,FixedBandwidth> SingleR2GaussianFixed;
typedef MeanshiftMex<float,R2,Epanechnikov,FixedBandwidth> SingleR2EpanechnikovFixed;
typedef MeanshiftMex<double,R2,Gaussian,FixedBandwidth> DoubleR2GaussianFixed;
typedef MeanshiftMex<double,R2,Epanechnikov,FixedBandwidth> DoubleR2EpanechnikovFixed;

typedef MeanshiftMex<float,R2,Gaussian,SampleVariableBandwidth> SingleR2GaussianSample;
typedef MeanshiftMex<float,R2,Epanechnikov,SampleVariableBandwidth> SingleR2EpanechnikovSample;
typedef MeanshiftMex<double,R2,Gaussian,SampleVariableBandwidth> DoubleR2GaussianSample;
typedef MeanshiftMex<double,R2,Epanechnikov,SampleVariableBandwidth> DoubleR2EpanechnikovSample;

typedef MeanshiftMex<float,R3,Gaussian,FixedBandwidth> SingleR3GaussianFixed;
typedef MeanshiftMex<float,R3,Epanechnikov,FixedBandwidth> SingleR3EpanechnikovFixed;
typedef MeanshiftMex<double,R3,Gaussian,FixedBandwidth> DoubleR3GaussianFixed;
typedef MeanshiftMex<double,R3,Epanechnikov,FixedBandwidth> DoubleR3EpanechnikovFixed;

typedef MeanshiftMex<float,R6,Gaussian,FixedBandwidth> SingleR6GaussianFixed;
typedef MeanshiftMex<float,R6,Epanechnikov,FixedBandwidth> SingleR6EpanechnikovFixed;
typedef MeanshiftMex<double,R6,Gaussian,FixedBandwidth> DoubleR6GaussianFixed;
typedef MeanshiftMex<double,R6,Epanechnikov,FixedBandwidth> DoubleR6EpanechnikovFixed;

//typedef MeanshiftMex<float,Gaussian,S2,FixedBandwidth> SingleGaussianS2Fixed;
//typedef MeanshiftMex<float,Epanechnikov,S2,FixedBandwidth> SingleEpanechnikovS2Fixed;
//typedef MeanshiftMex<double,Gaussian,S2,FixedBandwidth> DoubleGaussianS2Fixed;
//typedef MeanshiftMex<double,Epanechnikov,S2,FixedBandwidth> DoubleEpanechnikovS2Fixed;
//
//typedef MeanshiftMex<double,Gaussian,S2,BalloonVariableBandwidth> DoubleGaussianS2Balloon;
//typedef MeanshiftMex<double,Epanechnikov,S2,BalloonVariableBandwidth> DoubleEpanechnikovS2Balloon;
//typedef MeanshiftMex<float,Gaussian,S2,BalloonVariableBandwidth> SingleGaussianS2Balloon;
//typedef MeanshiftMex<float,Epanechnikov,S2,BalloonVariableBandwidth> SingleEpanechnikovS2Balloon;
//
//typedef MeanshiftMex<double,Gaussian,S2,SampleVariableBandwidth> DoubleGaussianS2Sample;
//typedef MeanshiftMex<double,Epanechnikov,S2,SampleVariableBandwidth> DoubleEpanechnikovS2Sample;
//typedef MeanshiftMex<float,Gaussian,S2,SampleVariableBandwidth> SingleGaussianS2Sample;
//typedef MeanshiftMex<float,Epanechnikov,S2,SampleVariableBandwidth> SingleEpanechnikovS2Sample;

typedef std::tuple<std::string, 
		   std::string, 
		   std::string, 
		   std::string> FactoryKey;

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

typedef boost::function<MeanshiftMexIface*(const InputArguments& input)> FactoryProduct;

std::unordered_map<FactoryKey,FactoryProduct,FactoryHash> factory = {
    { std::make_tuple("single","r2","epanechnikov","fixed"), 
      std::bind(boost::factory<SingleR2EpanechnikovFixed*>(),std::placeholders::_1) },
    { std::make_tuple("single","r2","gaussian","fixed"), 
      std::bind(boost::factory<SingleR2GaussianFixed*>(),std::placeholders::_1) },
    { std::make_tuple("double","r2","epanechnikov","fixed"), 
      std::bind(boost::factory<DoubleR2EpanechnikovFixed*>(),std::placeholders::_1) },
    { std::make_tuple("double","r2","gaussian","fixed"), 
      std::bind(boost::factory<DoubleR2GaussianFixed*>(),std::placeholders::_1) },

    { std::make_tuple("single","r2","epanechnikov","sample"), 
      std::bind(boost::factory<SingleR2EpanechnikovSample*>(),std::placeholders::_1) },
    { std::make_tuple("single","r2","gaussian","sample"), 
      std::bind(boost::factory<SingleR2GaussianSample*>(),std::placeholders::_1) },
    { std::make_tuple("double","r2","epanechnikov","sample"), 
      std::bind(boost::factory<DoubleR2EpanechnikovSample*>(),std::placeholders::_1) },
    { std::make_tuple("double","r2","gaussian","sample"), 
      std::bind(boost::factory<DoubleR2GaussianSample*>(),std::placeholders::_1) },

    { std::make_tuple("single","r3","epanechnikov","fixed"), 
      std::bind(boost::factory<SingleR3EpanechnikovFixed*>(),std::placeholders::_1) },
    { std::make_tuple("single","r3","gaussian","fixed"), 
      std::bind(boost::factory<SingleR3GaussianFixed*>(),std::placeholders::_1) } ,
    { std::make_tuple("double","r3","epanechnikov","fixed"), 
      std::bind(boost::factory<DoubleR3EpanechnikovFixed*>(),std::placeholders::_1) },
    { std::make_tuple("double","r3","gaussian","fixed"), 
      std::bind(boost::factory<DoubleR3GaussianFixed*>(),std::placeholders::_1) }, 

    { std::make_tuple("single","r6","epanechnikov","fixed"), 
      std::bind(boost::factory<SingleR6EpanechnikovFixed*>(),std::placeholders::_1) },
    { std::make_tuple("single","r6","gaussian","fixed"), 
      std::bind(boost::factory<SingleR6GaussianFixed*>(),std::placeholders::_1) } ,
    { std::make_tuple("double","r6","epanechnikov","fixed"), 
      std::bind(boost::factory<DoubleR6EpanechnikovFixed*>(),std::placeholders::_1) },
    { std::make_tuple("double","r6","gaussian","fixed"), 
      std::bind(boost::factory<DoubleR6GaussianFixed*>(),std::placeholders::_1) } 

    //  { std::make_tuple("double","epanechnikov","s2","fixed"), 
    //    std::bind(boost::factory<DoubleEpanechnikovS2Fixed*>(),std::placeholders::_1) },
    //  { std::make_tuple("double","gaussian","s2","fixed"), 
    //    std::bind(boost::factory<DoubleGaussianS2Fixed*>(),std::placeholders::_1) },
    //  { std::make_tuple("single","epanechnikov","s2","fixed"), 
    //    std::bind(boost::factory<SingleEpanechnikovS2Fixed*>(),std::placeholders::_1) },
    //  { std::make_tuple("single","gaussian","s2","fixed"), 
    //    std::bind(boost::factory<SingleGaussianS2Fixed*>(),std::placeholders::_1) } ,
    //
    //  { std::make_tuple("double","epanechnikov","s2","balloon"), 
    //    std::bind(boost::factory<DoubleEpanechnikovS2Balloon*>(),std::placeholders::_1) },
    //  { std::make_tuple("double","gaussian","s2","balloon"), 
    //    std::bind(boost::factory<DoubleGaussianS2Balloon*>(),std::placeholders::_1) },
    //  { std::make_tuple("single","epanechnikov","s2","balloon"), 
    //    std::bind(boost::factory<SingleEpanechnikovS2Balloon*>(),std::placeholders::_1) },
    //  { std::make_tuple("single","gaussian","s2","balloon"), 
    //    std::bind(boost::factory<SingleGaussianS2Balloon*>(),std::placeholders::_1) },
    //
    //  { std::make_tuple("double","epanechnikov","s2","sample"), 
    //    std::bind(boost::factory<DoubleEpanechnikovS2Sample*>(),std::placeholders::_1) },
    //  { std::make_tuple("double","gaussian","s2","sample"), 
    //    std::bind(boost::factory<DoubleGaussianS2Sample*>(),std::placeholders::_1) },
    //  { std::make_tuple("single","epanechnikov","s2","sample"), 
    //    std::bind(boost::factory<SingleEpanechnikovS2Sample*>(),std::placeholders::_1) },
    //  { std::make_tuple("single","gaussian","s2","sample"), 
    //    std::bind(boost::factory<SingleGaussianS2Sample*>(),std::placeholders::_1) }  
}; 

#endif
