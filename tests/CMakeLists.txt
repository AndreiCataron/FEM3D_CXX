set(CMAKE_CXX_STANDARD 20)

#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE Release)
#endif()
#
#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_library(FOO_LIB gmsh)

find_path(FOO_PATH gmsh.h)

find_package(OpenMP REQUIRED)

add_executable(sbc test_setBoundaryConditions.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp ../include/Mesh.hpp ../src/Mesh.cpp ../include/params.hpp)
add_executable(cnsbc test_checkNodeSatisfiesBoundaryCondition.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp ../include/Mesh.hpp ../src/Mesh.cpp ../include/params.hpp)
add_executable(cpi test_checkParamsInitialization.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp ../include/Mesh.hpp ../src/Mesh.cpp ../include/params.hpp)
add_executable(lj learn_jabobians.cpp ../include/utils.hpp ../src/utils.cpp)
add_executable(lnbag learn_isNodalBasisAutomaticallyGenerated.cpp)
add_executable(rtst randomCXXTests.cpp)
add_executable(ret randomEigenTests.cpp)
add_executable(stiff test_stiffnessAndLoad.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp ../include/Mesh.hpp ../src/Mesh.cpp ../include/params.hpp)
add_executable(exam testExampleSolutions.cpp ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp ../include/Mesh.hpp ../src/Mesh.cpp ../include/params.hpp)
add_executable(util test_utils.cpp ../include/utils.hpp ../src/utils.cpp)
add_executable(examnonlin ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp ../include/Mesh.hpp ../src/Mesh.cpp ../include/params.hpp testNonlinearSolution.cpp)
add_executable(alter ../include/FEM3D.hpp ../src/FEM3D.cpp ../include/FEM3Dvector.hpp ../src/FEM3Dvector.cpp ../include/LinearElasticity3D.hpp ../src/LinearElasticity3D.cpp ../include/utils.hpp ../src/utils.cpp ../include/Mesh.hpp ../src/Mesh.cpp ../include/params.hpp ../include/Algoritm.hpp ../src/Algoritm.cpp testAlternant.cpp)

target_link_libraries(sbc LINK_PUBLIC "${FOO_LIB}")

target_link_libraries(stiff PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(rtst PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(cpi PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(sbc PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(cnsbc PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(exam PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(examnonlin PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(alter PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(sbc PUBLIC "${FOO_PATH}")

target_link_libraries(cnsbc LINK_PUBLIC "${FOO_LIB}")

target_include_directories(cnsbc PUBLIC "${FOO_PATH}")

target_link_libraries(cpi LINK_PUBLIC "${FOO_LIB}")

target_include_directories(cpi PUBLIC "${FOO_PATH}")

target_link_libraries(lj LINK_PUBLIC "${FOO_LIB}")

target_include_directories(lj PUBLIC "${FOO_PATH}")

target_link_libraries(lnbag LINK_PUBLIC "${FOO_LIB}")

target_include_directories(lnbag PUBLIC "${FOO_PATH}")

target_link_libraries(stiff LINK_PUBLIC "${FOO_LIB}")

target_include_directories(stiff PUBLIC "${FOO_PATH}")

target_link_libraries(exam LINK_PUBLIC "${FOO_LIB}")

target_include_directories(exam PUBLIC "${FOO_PATH}")

target_link_libraries(examnonlin LINK_PUBLIC "${FOO_LIB}")

target_include_directories(examnonlin PUBLIC "${FOO_PATH}")

target_link_libraries(alter LINK_PUBLIC "${FOO_LIB}")

target_include_directories(alter PUBLIC "${FOO_PATH}")

add_test(NAME testSetBoundaryConditions COMMAND sbc)
add_test(NAME testCheckNodeSatisfiesBoundaryCondition COMMAND cnsbc)
add_test(NAME testCheckParamsInitialization COMMAND cpi)
add_test(NAME learnJacobians COMMAND lj)
add_test(NAME learnNodalBasisAutomaticallyGenerated COMMAND lnbag)
add_test(NAME randomTst COMMAND rtst)
add_test(NAME randomEigenTst COMMAND ret)
add_test(NAME testStiffnessAndLoad COMMAND stiff)
add_test(NAME testExamples COMMAND exam)
add_test(NAME testUtils COMMAND util)
add_test(NAME testNonlinear COMMAND examnonlin)
add_test(NAME testAlternant COMMAND alter)