#include "../include/FEM3D.hpp"

#include <gmsh.h>

FEM3D::FEM3D(const FEM3D::Params &params) : params_(params) {}

void FEM3D::setBoundaryConditions() {
    std::vector<std::pair<int, int> > domain_entity;
    gmsh::model::getEntities(domain_entity, 3);

    std::vector<std::pair<int, int> > boundary;
    gmsh::model::getBoundary(domain_entity, boundary, true, false, false);

    for(auto b : boundary){
        std::vector<std::size_t> tags;
        std::vector<double> coord, paramCoords;
        gmsh::model::mesh::getNodes(tags, coord, paramCoords, b.first, b.second, true, false);
    }
}