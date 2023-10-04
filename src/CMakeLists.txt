set(CMAKE_CXX_STANDARD 20)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

add_executable(srctarget ../include/LinearElasticity3D.hpp LinearElasticity3D.cpp main.cpp ../include/FEM3D.hpp FEM3D.cpp ../include/FEM3Dvector.hpp FEM3Dvector.cpp ../include/utils.hpp utils.cpp ../include/Mesh.hpp Mesh.cpp)

find_package(OpenMP REQUIRED)

target_link_libraries(srctarget PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(srctarget LINK_PUBLIC "${FOO_LIB}")

target_include_directories(srctarget PUBLIC "${FOO_PATH}")

set(EXPRTK_DIR /Users/andrei/CLionProjects/FEM/include/exprtk)

include_directories(${EXPRTK_DIR})

target_include_directories(srctarget PUBLIC ${EXPRTK_DIR})

target_include_directories(srctarget PUBLIC include)
