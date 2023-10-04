#include "../include/Mesh.hpp"
#include <gmsh.h>

Mesh::Mesh(int argc, char **argv) {
    gmsh::initialize(argc, argv);
}

Mesh::~Mesh() {
    gmsh::finalize();
}