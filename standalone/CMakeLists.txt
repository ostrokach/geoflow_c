
cmake_minimum_required (VERSION 2.6)

set( CMAKE_VERBOSE_MAKEFILE on )

project (geoflow_standalone)

#FIND_PACKAGE(Eigen)

include_directories( /usr/local/Cellar/eigen/3.2.1/include/eigen3/ )

add_executable(osm main.cpp
     cpbconcz2.cpp
     solvationsub.cpp
     modules.cpp
     subfunction.cpp
     pbsolvercz.cpp
     surfconcz.cpp
)

