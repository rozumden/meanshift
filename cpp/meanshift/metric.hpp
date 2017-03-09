#ifndef __MEANSHIFT_METRIC_HPP__
#define __MEANSHIFT_METRIC_HPP__

template<typename Chart> 
struct Metric
{
  static 
  typename chart_traits<Chart>::Scalar 
  norm(const typename chart_traits<Chart>::TangentVector& u) 
  { 
    return sqrt(norm_sq(u));
  }

  static 
  typename chart_traits<Chart>::Scalar 
  norm_sq(const typename chart_traits<Chart>::TangentVector& u) 
  { 
    return u.adjoint()*u;
  }

  static 
  typename chart_traits<Chart>::Scalar 
  dist(const typename chart_traits<Chart>::EmbeddedVector& x, 
       const typename chart_traits<Chart>::EmbeddedVector& y) 
  { 
    return sqrt(dist_sq(x,y));
  }  

  static 
  typename chart_traits<Chart>::Scalar 
  dist_sq(const typename chart_traits<Chart>::EmbeddedVector& x, 
	  const typename chart_traits<Chart>::EmbeddedVector& y) 
  { 
    typename chart_traits<Chart>::TangentVector u = Chart::logm(x,y);
    return norm_sq(u);
  }  
};

#endif
