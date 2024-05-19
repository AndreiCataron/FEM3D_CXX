set(CMAKE_CXX_STANDARD 20)

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

find_package(OpenMP REQUIRED)

add_executable(srctarget ../include/FEM3D.hpp FEM3D.cpp ../include/FEM3Dvector.hpp FEM3Dvector.cpp ../include/LinearElasticity3D.hpp LinearElasticity3D.cpp main.cpp ../include/utils.hpp utils.cpp ../include/Mesh.hpp Mesh.cpp ../include/params.hpp ../include/AlternantKozlov.hpp AlternantKozlov.cpp ../include/eigen_common.hpp ../include/EigenLogger.hpp)

target_link_libraries(srctarget LINK_PUBLIC "${FOO_LIB}")

target_link_libraries(srctarget PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(srctarget PUBLIC "${FOO_PATH}")

set(EXPRTK_DIR /Users/andrei/CLionProjects/FEM/include/exprtk)

include_directories(${EXPRTK_DIR})

target_include_directories(srctarget PUBLIC ${EXPRTK_DIR})

target_include_directories(srctarget PUBLIC include)
