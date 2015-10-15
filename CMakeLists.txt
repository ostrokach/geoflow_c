cmake_minimum_required(VERSION 3.1.0)

enable_testing()

set(CMAKE_VERBOSE_MAKEFILE on)

include(ExternalProject)

project(geoflow_c CXX)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

# TODO: Go download boost if it's not found.
find_package(Boost)
include_directories(${Boost_INCLUDE_DIR})

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/src/cmake)
FIND_PACKAGE(Eigen3)

if (NOT EIGEN3_INCLUDE_DIR)
	# I can't seem to get this to work properly.
	# CMake complains about the add_dependencies line.
	#  -k*
	ExternalProject_add( eigen3
		URL https://bitbucket.org/eigen/eigen/get/3.2.5.tar.gz
		TIMEOUT 120
		LOG_DOWNLOAD ON)
	ExternalProject_Get_Property(eigen3 source_dir)
	include_directories(${source_dir})
	add_dependencies(GeometricFlow, eigen3)
else()
	include_directories(${EIGEN3_INCLUDE_DIR})
endif()

add_subdirectory(src)
add_subdirectory(test)
