cmake_minimum_required(VERSION 3.1.0)

project (geoflow_c)

add_executable(GeometricFlow 
     runGeometricFlow.cpp
     GeometricFlow.cpp
     Atom.cpp
)

target_link_libraries(GeometricFlow)
