#ifndef MEANSHIFT_RIEMANNIAN_METRIC_HPP_
#define MEANSHIFT_RIEMANNIAN_METRIC_HPP_

template<typename _Chart> 
class Metric : public _Chart 
{
public:
  typedef typename _Chart::TangentVector TangentVector;
  typedef typename _Chart::EmbeddedVector EmbeddedVector;
  typedef typename _Chart::Scalar Scalar;
  typedef _Chart Chart;

public:
  static Scalar norm(const TangentVector& delta) 
  { 
    return sqrt(delta.adjoint()*delta);
  }

  static Scalar sq_norm(const TangentVector& delta) 
  { 
    return delta.adjoint()*delta;
  }

//  template <>
//  Scalar dist<Chart>(const EmbeddedVector& u, const EmbeddedVector& v) const
//  {
//    return Chart::dist(u,v);
//  }

  // static Scalar dist(const EmbeddedVector& u, const EmbeddedVector& v)
  // {
  //   TangentVector x = Chart::logm(u,v);
  //   return norm(x);
  // }

 // template <>
 // Scalar sq_dist<Chart>(const EmbeddedVector& u, const EmbeddedVector& v) const
 // { 
 //   return T::dist(u,v);
 // }

  static Scalar sq_dist2(const EmbeddedVector& u, const EmbeddedVector& v) 
  { 
    TangentVector x = Chart::logm(u,v);
    return sq_norm(x);
  }  
};

#endif
