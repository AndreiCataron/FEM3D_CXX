#include <iostream>
#include <cassert>
#include "../include/FEM3D.hpp"
#include <gmsh.h>

int main(int argc, char **argv) {
    FEM3D::Params par = {
            0,
            "x == 0",
            "x + y + 2 * z"
    };

    FEM3D fem(par);

    gmsh::initialize(argc, argv);

    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.5);
    gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::mesh::generate(3);

    gmsh::finalize();
}