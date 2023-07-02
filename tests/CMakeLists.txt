cmake_minimum_required(VERSION 3.25)

project(TESTS)

set(CMAKE_CXX_STANDARD 20)

add_executable(TESTS test_setBoundaryConditions.cpp)
add_test(NAME testSetBoundsryConditions COMMAND TESTS)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

target_link_libraries(TESTS PRIVATE "${FOO_LIB}")

target_include_directories(TESTS PRIVATE "${FOO_PATH}")
