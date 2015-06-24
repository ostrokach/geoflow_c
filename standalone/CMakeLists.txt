
cmake_minimum_required (VERSION 2.6)

set( CMAKE_VERBOSE_MAKEFILE on )

project (geoflow_standalone)

SET(geoflow_standalone_CMAKE_DIR "${geoflow_standalone_SOURCE_DIR}/cmake")

set(GOOGLE_LIBS -L$(GMOCK_HOME)/osx -L$(GMOCK_HOME)/osx/gtest 
           -lgtest -lgmock)

set(GOOGLE_INCLUDE -I$(GMOCK_HOME)/include 
	   -I$(GMOCK_HOME)/gtest/include -std=c++11)

include_directories( /usr/local/Cellar/eigen/3.2.4/include/eigen3/ 
                     /usr/local/Cellar/boost/1.58.0/include 
#		     ${GOOGLE_INCLUDE}
                    )

#target_link_libraries(osm
#                      ${GOOGLE_LIBS}
#			)

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
     runGeometricFlow.cpp
     GeometricFlow.cpp
     Atom.cpp
)
