#include "../include/FEM3Dvector.hpp"

#include <gmsh.h>
#include <fstream>
#include "../include/utils.hpp"
#include <iostream>

FEM3DVector::FEM3DVector(std::shared_ptr<ParamsVector> const &params) : params3d_(params), FEM3D(params) {}

FEM3DVector::FEM3DVector(std::shared_ptr<ParamsVector> const &params, Mesh &msh) : params3d_(params), FEM3D(params, msh){}

void FEM3DVector::setBoundaryConditions() {
    setNeumannBoundaryConditions();

    for (auto b : mesh.elems.boundary) {
        std::vector<std::size_t> tags;
        std::vector<double> coord, param_coords;
        gmsh::model::mesh::getNodes(tags, coord, param_coords, b.first, b.second, true, false);

        for (int i = 0; i < tags.size(); i++) {
            std::size_t tag = tags[i];
            int bc = checkNodeSatisfiesBoundaryEquation(coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);

            if (bc == 1) {
                if (dirichlet_bc.find(tag) == dirichlet_bc.end()) {
                    std::vector<double> prescribed_condition;
                    prescribed_condition.reserve(3);

                    prescribed_condition = params3d_ -> g(coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);

                    dirichlet_bc[tag] = prescribed_condition;
                }
            }
            // add here other BC
        }
    }

    indexConstrainedNodes();
    indexFreeNodes();
}

void FEM3DVector::indexConstrainedNodes() {
    for (const auto& n : dirichlet_bc) {
        // check node is not already in node indexes
        if (nodeIndexes.find(n.first) == nodeIndexes.end()) {
            constrainedNodes.emplace_back(nodeIndexes.size());
            nodeIndexes[n.first] = int(nodeIndexes.size());
        }
    }
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
    for (const auto& [tag, coord] : mesh.elems.node_coordinates) {
        myFile << tag << ' ' << std::get<0>(coord) << ' ' << std::get<1>(coord) << ' ' << std::get<2>(coord) << ' ';
    }
    myFile << '\n';

    // output displacements

    for (const auto& [tag, idx] : nodeIndexes) {
        myFile << tag << ' ' << displacements(3 * idx) << ' ' << displacements(3 * idx + 1) << ' ' << displacements(3 * idx + 2) << ' ';
    }

//    for (auto d : displacements) {
//        myFile << d << ' ';
//    }

    myFile.close();
}

double FEM3DVector::computeL2Error() {
    double result = 0;

    for (int i = 0; i < mesh.elems.elementTags.size(); i++) {
        // get tags of nodes in current element
        std::vector<std::size_t> elementNodeTags = std::vector<std::size_t>(
                mesh.elems.nodeTags.begin() + i * mesh.elems.noNodesPerElement,
                mesh.elems.nodeTags.begin() + (i + 1) * mesh.elems.noNodesPerElement);
        double det = mesh.elems.determinants[i];

        double integral = 0;

        // quadrature
        for (int j = 0; j < mesh.elems.noIntegrationPoints; j++) {
            // the exact solution at the integration point
            std::vector<double> exact;
            exact.reserve(3);

            exact = params3d_ -> exact_solution(mesh.elems.globalCoord[i * mesh.elems.localCoord.size() + 3 * j],
                                                mesh.elems.globalCoord[i * mesh.elems.localCoord.size() + 3 * j + 1],
                                                mesh.elems.globalCoord[i * mesh.elems.localCoord.size() + 3 * j + 2]);

            // the approximate solution at the integration point
            std::vector<double> approxSolution = {0, 0, 0};

            for (int k = 0; k < mesh.elems.noNodesPerElement; k++) {
                std::size_t tag = elementNodeTags[k];
                int index = nodeIndexes[tag];

                approxSolution[0] += displacements(3 * index) *
                                     mesh.elems.basisFunctionsValues[j * mesh.elems.noNodesPerElement + k];
                approxSolution[1] += displacements(3 * index + 1) *
                                     mesh.elems.basisFunctionsValues[j * mesh.elems.noNodesPerElement + k];
                approxSolution[2] += displacements(3 * index + 2) *
                                     mesh.elems.basisFunctionsValues[j * mesh.elems.noNodesPerElement + k];
            }

            for (int t = 0; t < 3; t++) {
                integral += mesh.elems.weights[j] * (exact[t] - approxSolution[t]) * (exact[t] - approxSolution[t]);
            }

        }

        result += det * integral;
    }

    return sqrt(result);
}

std::unordered_map<std::size_t, std::vector<double> > FEM3DVector::getDirichletBC() {
    return dirichlet_bc;
}