cmake_minimum_required(VERSION 3.2)


##
## build with APBS 
##
#add_items(
#    SOURCES
#    cpbconcz2.cpp
#    modules.cpp
#    pbsolvercz.cpp
#    solvationsub.cpp
#    subfunction.cpp
#    surfconcz.cpp
#)
# 
#add_items(
#    EXTERNAL_HEADERS
#    cpbconcz2.h
#)
#
#add_sublibrary(geoflow) 
#
## Enable C++11 features
#if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"
#       OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#	add_definitions(-std=c++0x)
#endif()


#
#  build standalone
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET(geoflow_standalone_CMAKE_DIR "${geoflow_standalone_SOURCE_DIR}/cmake")

set(GOOGLE_LIBS -L$(GMOCK_HOME)/osx -L$(GMOCK_HOME)/osx/gtest 
           -lgtest -lgmock)

set(GOOGLE_INCLUDE -I$(GMOCK_HOME)/include 
	   -I$(GMOCK_HOME)/gtest/include -std=c++11)

include_directories( /usr/local/Cellar/eigen/3.2.4/include/eigen3/ 
                     /usr/local/Cellar/boost/1.58.0/include 
#		     ${GOOGLE_INCLUDE}
                    )

add_library(GeometricFlowLib 
     GeometricFlowWrap.cpp
     GeometricFlow.cpp
     Atom.cpp
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

add_executable(GeometricFlowWrap
     runGeometricFlowWrap.c
)

target_link_libraries(GeometricFlowWrap GeometricFlowLib)
