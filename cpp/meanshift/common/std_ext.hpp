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
#ifndef MEANSHIFT_STD_EXT_HPP_
#define MEANSHIFT_STD_EXT_HPP_

#include <eigen3/Eigen/Dense>

namespace std
{
  template <> struct hash<Eigen::Matrix<double,3,1>>
  {
    typedef size_t result_type;
    typedef Eigen::Matrix<double,3,1> argument_type;

    size_t operator()(Eigen::Matrix<double,3,1> const & uri) const noexcept
    {
        return std::hash<double>()(uri[0]) ^
               std::hash<double>()(uri[1]) ^
               std::hash<double>()(uri[2]);
    }
  };

  template <> struct hash<Eigen::Matrix<float,3,1>>
  {
    typedef size_t result_type;
    typedef Eigen::Matrix<float,3,1> argument_type;

    size_t operator()(Eigen::Matrix<float,3,1> const & uri) const noexcept
    {
        return std::hash<float>()(uri[0]) ^
               std::hash<float>()(uri[1]) ^
               std::hash<float>()(uri[2]);
    }
  };
} // namespace std

#endif
