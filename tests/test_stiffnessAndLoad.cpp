#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include "../include/utils.hpp"
#include <gmsh.h>
#include <memory>
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>

int main(int argc, char **argv) {
    return 0;
//    auto par = std::make_shared<ParamsLE>(ParamsLE{
//            0.3, // h,
//            10,
//            1,
//            "0 == 0", // dirichlet BC
//            1, // quadrature precision
//            2, // order of lagrange polynomials
//            { "0", "0", "0" }, // f
//            { "x", "y", "z" }, // g
//            -1, // lambda
//            -1, // mu
//            0.34, // nu
//            12.864 // E
//    });
//
//    utils::checkParamsLE(*par);
//
//    Mesh msh(argc, argv);
//
//    LinearElasticity3D fem(*par, msh);
//
//    gmsh::initialize(argc, argv);
//
//    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
//    gmsh::model::occ::synchronize();
//
//    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.2);
//    gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
//    gmsh::option::setNumber("General.Verbosity", 1);
//
//    fem.setBoundaryConditions();
//
//    auto start = std::chrono::steady_clock::now();
//
//    fem.computeStiffnessMatrixAndLoadVector();
//
//    auto end = std::chrono::steady_clock::now();
//    auto diff = end - start;
//
//    fem.solveDisplacements();
//
//    std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
//
//    Eigen::SparseMatrix<double> sm = fem.getStiffnessMatrix();
//    Eigen::VectorXd lv = fem.getLoadVector();
//
//    std::cout << "Rows: " << sm.rows() << "; Cols: " << sm.cols() << "; Non-zero elements: " << sm.nonZeros() << '\n';
//    std::cout << "Load Vector Size: " << lv.size() << '\n';
//
//    // check stiffness is spd
//
//    std::cout << "Is sym: " << sm.isApprox(sm.transpose()) << '\n';
//
////    Eigen::LLT<Eigen::MatrixXd> lltOfA(sm); // compute the Cholesky decomposition of A
////    if(lltOfA.info() == Eigen::NumericalIssue)
////    {
////        std::cout << "Gresit";
////    }
////    else {
////        std::cout << "OK";
////    }
////    std::cout << '\n';
//
//    // check number of nodes in mesh
//
//    std::vector<std::size_t> nodeTags;
//    std::vector<double> coord, paramCoord;
//    gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord, -1, -1, true, false);
////    for (auto i : nodeTags) std::cout << i << ' ';
////    std::cout << '\n';
//
//    std::sort(nodeTags.begin(), nodeTags.end());
//    nodeTags.erase(std::unique(nodeTags.begin(), nodeTags.end()), nodeTags.end());
//
//    std::unordered_map<std::size_t, std::vector<double> > dir = fem.getDirichletBC();
//
//    std::cout << "No nodes: " << nodeTags.size() << ' ' << nodeTags.size() << '\n';
//    std::cout << "No of indexes: " << fem.getFreeNodes().size() + fem.getConstrainedNodes().size() << '\n';
//    std::cout << "No of constrained nodes: " << fem.getConstrainedNodes().size() << ' ' << dir.size() << '\n';
//
//    std::cout << fem.getNodeCoordinates().size();
//    fem.outputData("/Users/andrei/CLionProjects/FEM/outputs/out.txt");
//
////    for (auto p : fem.getNodeIndexes()) {
////        std::cout << p.first << ' ' << p.second << '\n';
////    }
//
////    std::vector<int> free = fem.getFreeNodes();
////    std::vector<int> constrained = fem.getConstrainedNodes();
////    for (auto i : constrained) std::cout << i << ' ';
////    std::cout << '\n';
////    for (auto i : free) std::cout << i << ' ';
////    std::cout << '\n';
////
////    std::unordered_map<std::size_t, int> ni = fem.getNodeIndexes();
////    for (const auto& n : ni){
////        std::cout << n.first << " " << n.second << '\n';
////    }
////
////    auto it = ni.begin();
////    while (it != ni.end()) {
////        std::cout << it->first << ' ' << it->second << '\n';
////        it++;
////    }
//
////    for (const auto& n : dir) {
////        std::cout << n.first << ' ';
////        for (auto i : n.second) {
////            std::cout << i << ' ';
////        }
////        std::cout << '\n';
////    }
//
//    //Eigen::VectorXd prod = sm * lv;
//
////    std::set<std::string> args(argv, argv + argc);
////    if(!args.count("-nopopup")) gmsh::fltk::run();
//
//
//    gmsh::finalize();
}