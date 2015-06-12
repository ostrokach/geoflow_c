
cmake_minimum_required (VERSION 2.6)

set( CMAKE_VERBOSE_MAKEFILE on )

project (geoflow_standalone)

SET(geoflow_standalone_CMAKE_DIR "${geoflow_standalone_SOURCE_DIR}/cmake")

include_directories( /usr/local/Cellar/eigen/3.2.4/include/eigen3/ 
                     /usr/local/Cellar/boost/1.58.0/include )
#include_directories( /Users/d3k084/PNNL_rep/apbs_finish-geoflow/apbs-pdb2pqr/apbs/include/ )

#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${geoflow_standalone_CMAKE_DIR}")
#FIND_PACKAGE(Eigen3 REQUIRED)
#include_directories( ${EIGEN3_INCLUDE_DIR} )

add_executable(osm main.cpp
     cpbconcz2.cpp
     solvationsub.cpp
     modules.cpp
     subfunction.cpp
     pbsolvercz.cpp
     surfconcz.cpp
)

add_executable(GeometricFlow 
     runGeometricFlow.C
     GeometricFlow.C
)
