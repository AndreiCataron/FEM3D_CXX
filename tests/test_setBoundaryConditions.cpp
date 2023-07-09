#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include <gmsh.h>
#include <unordered_map>
#include <cassert>
#include <typeinfo>

int main(int argc, char **argv) {
    FEM3DVector::ParamsVector par = {
            0,
            "x == 0",
            {"2 * x", "3 * y + 1", "x + z"}
    };

    LinearElasticity3D fem(par);

    gmsh::initialize(argc, argv);

    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.5);
    gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::mesh::generate(3);

    fem.setBoundaryConditions();

    std::unordered_map<std::size_t, std::vector<double> > map = fem.getDirichletBC();

    std::cout << "start";
    for (const auto& n : map){
        auto tag = n.first;
        std::vector<double> coord, paramcoord;
        int dim, t;
        gmsh::model::mesh::getNode(tag, coord, paramcoord, dim, t);
        assert(abs(n.second[0] - 2 * coord[0]) < 1e-8);
        assert(abs(n.second[1] - 3 * coord[1] - 1) < 1e-8);
        assert(abs(n.second[2] - coord[0] - coord[2]) < 1e-8);
    }


    gmsh::finalize();

    return 0;
}