cmake_minimum_required(VERSION 3.25)
enable_testing()

project(FEM)

add_subdirectory(tests)

set(CMAKE_CXX_STANDARD 20)

set(EXPRTK_DIR /Users/andrei/CLionProjects/FEM/include/exprtk)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

find_package(OpenMP REQUIRED)

add_executable(FEM src/main.cpp include/FEM3D.hpp src/FEM3D.cpp include/FEM3Dvector.hpp src/FEM3Dvector.cpp include/LinearElasticity3D.hpp src/LinearElasticity3D.cpp include/utils.hpp src/utils.cpp)

target_link_libraries(FEM LINK_PUBLIC "${FOO_LIB}")

target_link_libraries(FEM PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(FEM PUBLIC "${FOO_PATH}")

target_include_directories(FEM PUBLIC ${EXPRTK_DIR})

target_include_directories(FEM PUBLIC include)

include_directories(${EXPRTK_DIR})


