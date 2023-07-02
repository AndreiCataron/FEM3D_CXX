#include "FEM3D.hpp"

#include <gmsh.h>

FEM3D::FEM3D(const FEM3D::Params &params) : params_(params) {}

void FEM3D::setBoundaryConditions() {
    std::vector<std::pair<int, int> > domain_entity;
    gmsh::model::getEntities(domain_entity, 3);

    std::vector<std::pair<int, int> > boundary;
    gmsh::model::getBoundary(domain_entity, boundary, true, false, false);

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t> > elementTags;
    std::vector<std::vector<std::size_t> > nodeTags;
}