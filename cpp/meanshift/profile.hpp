#ifndef MEANSHIFT_PROFILE_HPP_
#define MEANSHIFT_PROFILE_HPP_

struct Epanechnikov
{
  template <typename T> 
  static T k(const T& x)
  {
    return fabs(x) > 1.0 ? 0.0 : 1.0-x;
  }

  template <typename T> 
  static T g(const T& x) 
  {
    return fabs(x) > 1.0 ? 0.0 : 1.0;
  }

};

struct Gaussian
{
  template <typename T> 
  static T k(const T& x)
  {
    return fabs(x) > 10.0 ? 0.0 : exp(-0.5*x);
  }

  template <typename T> 
  static T g(const T& x) 
  {
    return fabs(x) > 3.0 ? 0.0 : exp(-0.5*x);
  }
};

#endif
