cmake_minimum_required(VERSION 3.25)

project(TESTS)

set(CMAKE_CXX_STANDARD 20)

add_executable(TESTS test_setBoundaryConditions.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp test_checkNodeSatisfiesBoundaryCondition.cpp)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

link_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(TESTS PRIVATE "${FOO_LIB}")

target_include_directories(TESTS PRIVATE "${FOO_PATH}")

add_test(NAME testSetBoundaryConditions COMMAND TESTS)
add_test(NAME testCheckNodeSatisfiesBoundaryCondition COMMAND TESTS)