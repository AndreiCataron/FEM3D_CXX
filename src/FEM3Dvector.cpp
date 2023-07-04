#include "../include/FEM3Dvector.hpp"

#include <gmsh.h>

FEM3DVector::FEM3DVector(const FEM3DVector::ParamsVector &params) : params3d_(params), FEM3D(params){}

void FEM3DVector::setBoundaryConditions() {
    std::vector<std::pair<int, int> > domain_entity;
    gmsh::model::getEntities(domain_entity, 3);

    std::vector<std::pair<int, int> > boundary;
    gmsh::model::getBoundary(domain_entity, boundary, true, false, false);

    for(auto b : boundary){
        std::vector<std::size_t> tags;
        std::vector<double> coord, param_coords;
        gmsh::model::mesh::getNodes(tags, coord, param_coords, b.first, b.second, true, false);

        for(int i = 0; i < tags.size(); i++) {
            checkNodeSatisfiesBoundaryEquation(tags[i], coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);
        }
    }
}