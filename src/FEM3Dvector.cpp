#include "../include/FEM3Dvector.hpp"

#include <gmsh.h>
#include <fstream>

FEM3DVector::FEM3DVector(const FEM3DVector::ParamsVector &params) : params3d_(params), FEM3D(params){}

void FEM3DVector::setBoundaryConditions() {
    std::vector<std::pair<int, int> > domain_entity;
    gmsh::model::getEntities(domain_entity, 3);

    std::vector<std::pair<int, int> > boundary;
    gmsh::model::getBoundary(domain_entity, boundary, true, false, false);

    for (auto b : boundary) {
        std::vector<std::size_t> tags;
        std::vector<double> coord, param_coords;
        gmsh::model::mesh::getNodes(tags, coord, param_coords, b.first, b.second, true, false);

        for (int i = 0; i < tags.size(); i++) {
            std::size_t tag = tags[i];
            int bc = checkNodeSatisfiesBoundaryEquation(coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);

            if (bc == 1) {
                std::vector<double> prescribed_condition;
                prescribed_condition.reserve(3);

                for (const auto &cond : params3d_.g) {
                    prescribed_condition.emplace_back(parseExpression(cond, coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]));
                }
                dirichlet_bc[tag] = prescribed_condition;
            }

            // add here other BC

        }
    }

    indexConstrainedNodes();
    indexFreeNodes();
}

void FEM3DVector::indexConstrainedNodes() {
    for (const auto& n : dirichlet_bc) {
        constrainedNodes.emplace_back(nodeIndexes.size());
        nodeIndexes[n.first] = int(nodeIndexes.size());
    }
}

void FEM3DVector::getNodesCoordinates() {
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, true, false);

    for (int i = 0; i < nodeTags.size(); i++) {
        node_coordinates[nodeTags[i]] = {coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]};
    }
}

std::unordered_map<std::size_t, std::vector<double> > FEM3DVector::getDirichletBC() {
    return dirichlet_bc;
}

std::unordered_map<std::size_t, std::tuple<double, double, double> > FEM3DVector::getNodeCoordinates() {
    return node_coordinates;
}

void FEM3DVector::outputData(std::string file) {
    std::ofstream myFile;
    myFile.open(file, std::ios::out | std::ios::trunc);

    // output tag : index pairs
    for (const auto& [tag, idx] : nodeIndexes) {
        myFile << tag << ' ' << idx << ' ';
    }
    myFile << '\n';

    // output node coordinates as tag : coord1 coord2 coord3
    for (const auto& [tag, coord] : node_coordinates) {
        myFile << tag << ' ' << std::get<0>(coord) << ' ' << std::get<1>(coord) << ' ' << std::get<2>(coord) << ' ';
    }
    myFile << '\n';

    // output displacements
    for (auto d : displacements) {
        myFile << d << ' ';
    }

    myFile.close();
}