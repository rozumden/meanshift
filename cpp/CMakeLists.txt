cmake_minimum_required(VERSION 3.6)
project(mexplus)
find_package(Matlab REQUIRED)
matlab_add_mex(NAME meanshift SRC meanshift.mex.cpp OUTPUT_NAME meanshift)
target_compile_options(meanshift PRIVATE -std=c++11 -fno-omit-frame-pointer)
target_include_directories(meanshift PRIVATE . ~/src/meanshift/cpp ~/src/mexplus/include ~/src/mexplusplus/cpp)
install(TARGETS meanshift LIBRARY DESTINATION mex)