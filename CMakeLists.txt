cmake_minimum_required(VERSION 3.25)
enable_testing()

project(FEM)

add_subdirectory(tests)

set(CMAKE_CXX_STANDARD 20)

set(EXPRTK_DIR /Users/andrei/CLionProjects/FEM/include/exprtk)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

add_executable(FEM src/main.cpp ${EXPRTK_DIR}/exprtk.hpp include/FEM3D.hpp src/FEM3D.cpp)

target_link_libraries(FEM PRIVATE "${FOO_LIB}")

target_include_directories(FEM PRIVATE "${FOO_PATH}")

target_include_directories(FEM PRIVATE ${EXPRTK_DIR})

include_directories(${EXPRTK_DIR})


