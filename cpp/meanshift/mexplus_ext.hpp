#ifndef MEANSHIFT_MEXPLUS_EXT_HPP_
#define MEANSHIFT_MEXPLUS_EXT_HPP_

#include <mexplus.h>
#include <eigen3/Eigen/Dense>

template <typename Scalar>
using Vector3 = Eigen::Matrix<Scalar,3,1>; 

namespace mexplus {
  // Define two template specializations.
  template <>
  mxArray* 
  MxArray::from(const std::vector<Vector3<double> >& u) 
  {
    //    typedef typename Derived::Scalar Scalar;
    typedef double Scalar;
    size_t n = u.size();
    if (n > 0) {
      size_t m = u[0].rows();
      MxArray numeric(MxArray::Numeric<Scalar>(m,n));
      size_t j = 0;
      for (auto u_i : u) {
      	for  (size_t i=0;i<m;i++)  // rows
      	  numeric.set(i,j,u_i.derived().coeff(i,0));
      	j++;
      }
      return numeric.release();
    }
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Vector3<double> >* u) 
  {
    //typedef typename Derived::Scalar Scalar;
    // Write your conversion code. For example,
    typedef double Scalar;

    MxArray numeric(array);
    size_t m = numeric.rows();
    size_t n = numeric.cols();
    u->reserve(m);
    for (size_t j = 0; j < n; j++) {
      Vector3<double> u_i;
      for (size_t i = 0; i < m; i++) 
	      u_i[i] = numeric.at<double>(i,j);    
      u->push_back(u_i);
    }
  }

  template <>
  mxArray* 
  MxArray::from(const std::vector<Vector3<float> >& u) 
  {
    //    typedef typename Derived::Scalar Scalar;
    typedef float Scalar;
    size_t n = u.size();
    if (n > 0) {
      size_t m = u[0].rows();
      MxArray numeric(MxArray::Numeric<Scalar>(m,n));
      size_t j = 0;
      for (auto u_i : u) {
        for  (size_t i=0;i<m;i++)  // rows
          numeric.set(i,j,u_i.derived().coeff(i,0));
        j++;
      }
      return numeric.release();
    }
  }

  template <>
  void 
  MxArray::to(const mxArray* array, std::vector<Vector3<float> >* u) 
  {
    //typedef typename Derived::Scalar Scalar;
    // Write your conversion code. For example,
    typedef float Scalar;

    MxArray numeric(array);
    size_t m = numeric.rows();
    size_t n = numeric.cols();
    u->reserve(m);
    for (size_t j = 0; j < n; j++) {
      Vector3<float> u_i;
      for (size_t i = 0; i < m; i++) 
        u_i[i] = numeric.at<float>(i,j);    
      u->push_back(u_i);
    }
  }
} // namespace mexplus

#endif
