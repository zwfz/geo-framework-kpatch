cmake_minimum_required(VERSION 3.1.0)

project(geom)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|(Apple)?Clang")
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

add_library(geom STATIC vector.cc matrix3x3.cc bspline.cc trimesh.cc)
