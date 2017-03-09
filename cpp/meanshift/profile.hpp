#ifndef __MEANSHIFT_PROFILE_HPP__
#define __MEANSHIFT_PROFILE_HPP__

template <typename Derived>
struct Profile 
{
  typedef double T;
  //  template <typename T> 
  static T k(const T& x) 
  {
    return Derived::k(x);
  }

  //  template <typename T> 
  static T g(const T& x) 
  {
    return Derived::g(x);
  }
};

struct Epanechnikov : public Profile<Epanechnikov>
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

struct Gaussian : public Profile<Gaussian>
{
  template <typename T> 
  static T k(const T& x)
  {
    return fabs(x) > 3.0 ? 0.0 : exp(-0.5*x);
  }

  template <typename T> 
  static T g(const T& x) 
  {
    return Gaussian::k(x);
  }
};

#endif
