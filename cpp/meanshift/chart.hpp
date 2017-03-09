#ifndef MEANSHIFT_CHART_HPP
#define MEANSHIFT_CHART_HPP

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include <math.h>

template <typename T>
T sacos(const T& rho)
{
  T theta = acos(rho);
  if (rho < 0)
    return theta-M_PI;

  return theta;
}

typedef enum matrix_action {identity,inverse} MatrixAction;

template<typename EmbeddedVector>
Eigen::Matrix<typename EmbeddedVector::Scalar,3,3> Ru(const EmbeddedVector& u, const MatrixAction action) 
{
  typedef typename EmbeddedVector::Scalar Scalar;
  Scalar flip = 1.0;
  if (action == MatrixAction::inverse)
    flip = -1.0;

  if ((u-EmbeddedVector(-1,0,0)).norm() < std::numeric_limits<Scalar>::epsilon()) 
    return -1.0*Eigen::Matrix<Scalar,3,3>::Identity();  
  
  if ((u-EmbeddedVector(1,0,0)).norm() < std::numeric_limits<Scalar>::epsilon()) 
    return Eigen::Matrix<Scalar,3,3>::Identity();  

  EmbeddedVector ru(0,-u[2],u[1]);
  ru /= ru.norm();
  Eigen::AngleAxis<Scalar> R(flip*acos(u[0]),ru);
  
  return R.matrix();
}

template<typename _Scalar> class S2 
{
public:
  typedef _Scalar Scalar;
  typedef Eigen::Matrix<Scalar,2,1> TangentVector;
  typedef Eigen::Matrix<Scalar,3,1> EmbeddedVector;

  static const Scalar toleps;

public:
  static Scalar dist(EmbeddedVector u, EmbeddedVector v) 
  {
    u /= u.norm();
    v /= v.norm();
    Scalar s = u.adjoint()*v;
    Scalar val;
    if (fabs(s) < 1)  {
      val = sacos<Scalar>(s);
    } else {
      val = 0;
    }
    return val;
  }

  static Scalar sq_dist(const EmbeddedVector& u, const EmbeddedVector& v)
  { 
    Scalar val = dist(u,v);
    
    return val*val;
  };

  static EmbeddedVector expm(const EmbeddedVector& u, const TangentVector& x) 
  {
    Scalar nx = x.norm();
    EmbeddedVector v0;
    if (nx > std::numeric_limits<Scalar>::epsilon()) {
      v0[0] = cos(nx);
      v0[1] = (sin(nx)/nx)*x[0];
      v0[2] = (sin(nx)/nx)*x[1];
    } 
    else {
      v0 = EmbeddedVector(1.0,0.0,0.0);
    }

    Eigen::Matrix<Scalar,3,3> R = Ru(u,MatrixAction::identity);
    return R*v0;
  }

  static TangentVector logm(const EmbeddedVector& u, const EmbeddedVector& v) 
  {
    Eigen::Matrix<Scalar,3,3> R = Ru(u,MatrixAction::inverse);
    EmbeddedVector ve = R*v;

    Scalar norm = (ve.segment(1,2)).norm();
    if (norm > std::numeric_limits<Scalar>::epsilon())
      return sacos(ve[0])*ve.segment(1,2)/norm;
    else
      return TangentVector(0,0);
  }
};

#endif
