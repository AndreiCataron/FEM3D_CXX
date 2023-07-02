cmake_minimum_required(VERSION 3.25)
project(FEM)

set(CMAKE_CXX_STANDARD 17)

set(EXPRTK_DIR /Users/andrei/CLionProjects/FEM/exprtk)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

add_executable(FEM main.cpp ${EXPRTK_DIR}/exprtk.hpp FEM3D.hpp FEM3D.cpp)

target_link_libraries(FEM PRIVATE "${FOO_LIB}")

target_include_directories(FEM PRIVATE "${FOO_PATH}")

target_include_directories(FEM PRIVATE ${EXPRTK_DIR})

include_directories(${EXPRTK_DIR})