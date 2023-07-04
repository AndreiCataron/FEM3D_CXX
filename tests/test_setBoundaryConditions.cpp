#include <iostream>
#include "../include/FEM3D.hpp"
#include <gmsh.h>
#include <unordered_map>

int main(int argc, char **argv) {
//    FEM3DVector::ParamsVector par = {
//            0,
//            "x == 0",
//            "x + y + 2 * z"
//    };
//
//    FEM3D fem(par);
//
//    gmsh::initialize(argc, argv);
//
//    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
//    gmsh::model::occ::synchronize();
//
//    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.5);
//    gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
//    gmsh::option::setNumber("General.Verbosity", 1);
//    gmsh::model::mesh::generate(3);
//
//    fem.setBoundaryConditions();
//
//    std::unordered_map<std::size_t, double> map = fem.getDirichletBC();
//
//    std::cout << "start";
//    for (const auto& n : map){
//        auto tag = n.first;
//        std::vector<double> coord, paramcoord;
//        int dim, t;
//        gmsh::model::mesh::getNode(tag, coord, paramcoord, dim, t);
//        std::cout<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<" "<<n.second<<std::endl;
//    }
//
//
//    gmsh::finalize();

    return 0;
}