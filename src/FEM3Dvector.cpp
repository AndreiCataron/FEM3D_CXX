#include "../include/FEM3Dvector.hpp"

#include <gmsh.h>
#include <fstream>
#include <iostream>

FEM3DVector::FEM3DVector(std::shared_ptr<ParamsVector> const& params) : FEM3D(params), params3d_(params) {}

FEM3DVector::FEM3DVector(std::shared_ptr<ParamsVector> const& params, std::shared_ptr<Mesh> const& msh) : FEM3D(params, msh), params3d_(params){}

void FEM3DVector::setDirichletBoundaryConditions() noexcept {
    for (auto b : mesh -> bdry.boundary) {
        std::vector<std::size_t> tags;
        std::vector<double> coord, param_coords;
        gmsh::model::mesh::getNodes(tags, coord, param_coords, b.first, b.second, true, false);

        for (int i = 0; i < tags.size(); i++) {
            std::size_t tag = tags[i];
            int bc = checkNodeSatisfiesBoundaryEquation(coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);

            if (bc == 1) {
                if (dirichlet_bc.find(tag) == dirichlet_bc.cend()) {
                    std::vector<double> prescribed_condition;
                    prescribed_condition.reserve(3);

                    prescribed_condition = params3d_ -> g(coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);

                    dirichlet_bc[tag] = prescribed_condition;
                }
            }
            else {
                if (std::find(bdryFreeNodes.cbegin(), bdryFreeNodes.cend(), tag) == bdryFreeNodes.cend()) {
                    bdryFreeNodes.emplace_back(tag);
                }
            }
        }
    }
}

void FEM3DVector::setDirichletBoundaryConditions(DirichletMap conditionMap) noexcept {
    dirichlet_bc = conditionMap;
}

void FEM3DVector::indexConstrainedNodes() noexcept {
    for (const auto& n : dirichlet_bc) {
//        std::cout << nodeIndexes.size() << ' ' << n.first << '\n';
        // check node is not already in node indexes
        if (nodeIndexes.find(n.first) == nodeIndexes.cend()) {
            constrainedNodes.emplace_back(nodeIndexes.size());
            nodeIndexes[n.first] = int(nodeIndexes.size());
        }
    }
}

// if boundaryError is True, output a projection on the plane given by equation of the solution on the neumann boundary
void FEM3DVector::outputData(const std::string& file, bool boundaryError, const std::vector<double>& plane) {
    std::ofstream myFile;
    myFile.open(file, std::ios::out | std::ios::trunc);

    // Line 1
    // output tag : index pairs
    for (const auto& [tag, idx] : nodeIndexes) {
        myFile << tag << ' ' << idx << ' ';
    }
    myFile << '\n';

    // Line 2
    // output node coordinates as tag : coord1 coord2 coord3
    for (const auto& [tag, coord] : mesh -> global.node_coordinates) {
        myFile << tag << ' ' << std::get<0>(coord) << ' ' << std::get<1>(coord) << ' ' << std::get<2>(coord) << ' ';
    }
    myFile << '\n';

    // Line 3
    // output displacements

    for (const auto& [tag, idx] : nodeIndexes) {
        myFile << tag << ' ' << displacements(3 * idx) << ' ' << displacements(3 * idx + 1) << ' ' << displacements(3 * idx + 2) << ' ';
    }
    myFile << '\n';

    // Lines 4 and 5
    // neumann boundary error
    std::vector<std::size_t> boundaryNodesUnique = {};
    if (boundaryError) {
        for (auto coef : plane) myFile << coef << ' ';
        myFile << '\n';

        for (int i = 0; i < mesh -> bdry.boundaryFacesTags.size(); i++) {
            std::size_t faceTag = mesh -> bdry.boundaryFacesTags[i];
            if (neumannBoundaryTriangles.find(faceTag) != neumannBoundaryTriangles.end()) {
                std::vector<std::size_t> faceNodeTags = std::vector<std::size_t>(
                        mesh -> bdry.boundaryFacesNodes.begin() + i * mesh -> bdry.noNodesPerTriangle,
                        mesh -> bdry.boundaryFacesNodes.begin() + (i + 1) * mesh -> bdry.noNodesPerTriangle);
                boundaryNodesUnique.insert(boundaryNodesUnique.end(), faceNodeTags.begin(), faceNodeTags.end());
            }
        }
        utils::deleteDuplicatesFromVector(boundaryNodesUnique);

        for (auto tag : boundaryNodesUnique) {
            std::tuple<double, double, double> projection = utils::projectionOnPlane(mesh -> global.node_coordinates[tag], plane);

            myFile << tag << ' ' << std::get<0>(projection) << ' ' << std::get<1>(projection) << ' ' << std::get<2>(projection) << ' ';
        }
    }
    else {
        std::cout << '\n';
    }
    myFile << '\n';

    // output exact solution
    for (const auto& [tag, coord]: mesh -> global.node_coordinates) {
        std::vector<double> exact;
        exact.reserve(3);

        exact = params3d_ -> exact_solution(std::get<0>(coord), std::get<1>(coord), std::get<2>(coord));

        myFile << tag << ' ' << exact[0] << ' ' << exact[1] << ' ' << exact[2] << ' ';
    }
    myFile << '\n';

    // output exact solution on neumann boundary projected on plane
    if (boundaryError) {
        for (auto tag : boundaryNodesUnique) {
            std::tuple<double, double, double> coord = mesh -> global.node_coordinates[tag];

            std::vector<double> exact;
            exact.reserve(3);

            exact = params3d_ -> exact_solution(std::get<0>(coord), std::get<1>(coord), std::get<2>(coord));
            myFile << tag << ' ' << exact[0] << ' ' << exact[1] << ' ' << exact[2] << ' ';
        }
    }

    myFile.close();
}

void FEM3DVector::computeApproximateBoundaryGradients() {
//    for (const auto& [tag, idx] : nodeIndexes) {
//        std::cout << tag << ' ' << idx << '\n';
//    }
    approx_grads.reserve(mesh -> global.trianglesGlobalCoord.size() / 3);
    for (int i = 0; i < mesh -> global.trianglesGlobalCoord.size() / 3; i++) {
        auto tag = mesh -> elems.bdryAdjacentElems[i];
        if (auto it = std::find(mesh -> elems.elementTags.cbegin(), mesh -> elems.elementTags.cend(), tag);
                 it != mesh -> elems.elementTags.cend()) {

            auto idx = it - mesh -> elems.elementTags.cbegin();

            auto elementNodeTags = std::vector<std::size_t>(mesh -> elems.nodeTags.cbegin() + idx * mesh -> elems.noNodesPerElement,
                                                            mesh -> elems.nodeTags.cbegin() + (idx + 1) * mesh -> elems.noNodesPerElement);

            // the approximate solution gradient at the integration point
            Eigen::Matrix3d approxGradient = Eigen::Matrix3d::Zero();

            Eigen::Vector3d elementGrads;
            for (int k = 0; k < mesh -> elems.noNodesPerElement; k++) {
                Eigen::Vector3d referenceGrads;
                referenceGrads << mesh -> bdry.gradientsAtTriangleIntPoints[i * 3 * mesh -> elems.noNodesPerElement + 3 * k],
                                  mesh -> bdry.gradientsAtTriangleIntPoints[i * 3 * mesh -> elems.noNodesPerElement + 3 * k + 1],
                                  mesh -> bdry.gradientsAtTriangleIntPoints[i * 3 * mesh -> elems.noNodesPerElement + 3 * k + 2];

                elementGrads = mesh -> global.inverse_jacobians[idx] * referenceGrads;

                auto nodeTag = elementNodeTags[k];
                int index = nodeIndexes[nodeTag];

//                if (displacements(3 * index) == 0 && displacements(3 * index + 1) == 0 && displacements(3 * index + 2) == 0) {
//                    std::cout << nodeTag << ' ' << index << '\n';
//                }
                approxGradient.row(0) += displacements(3 * index) * elementGrads.transpose();
                approxGradient.row(1) += displacements(3 * index + 1) * elementGrads.transpose();
                approxGradient.row(2) += displacements(3 * index + 2) * elementGrads.transpose();
//                if (tag == 6223) {
//                    std::cout << "NOD: " << nodeTag << '\n' << elementGrads <<'\n'<< approxGradient << '\n';
//                    std::cout << displacements(3 * index) << ' ' << displacements(3 * index + 1) << ' ' << displacements(3 * index + 2) << '\n';
//                }
            }
//            if (approxGradient.isZero(0)) {
//                std::cout << tag << ' ';
//            }
            approx_grads.emplace_back(approxGradient);
        }
    }
}

void FEM3DVector::computeL2Error() {
    double result = 0;

    for (int i = 0; i < mesh -> elems.elementTags.size(); i++) {
        // get tags of nodes in current element
        std::vector<std::size_t> elementNodeTags = std::vector<std::size_t>(
                mesh -> elems.nodeTags.begin() + i * mesh -> elems.noNodesPerElement,
                mesh -> elems.nodeTags.begin() + (i + 1) * mesh -> elems.noNodesPerElement);
        double det = mesh -> global.determinants[i];

        double integral = 0;

        // quadrature
        for (int j = 0; j < mesh -> elems.noIntegrationPoints; j++) {
            // the exact solution at the integration point
            std::vector<double> exact;
            exact.reserve(3);

            exact = params3d_ -> exact_solution(mesh -> global.globalCoord[i * mesh -> elems.localCoord.size() + 3 * j],
                                                mesh -> global.globalCoord[i * mesh -> elems.localCoord.size() + 3 * j + 1],
                                                mesh -> global.globalCoord[i * mesh -> elems.localCoord.size() + 3 * j + 2]);

            // the approximate solution at the integration point
            std::vector<double> approxSolution = {0, 0, 0};

            for (int k = 0; k < mesh -> elems.noNodesPerElement; k++) {
                std::size_t tag = elementNodeTags[k];
                int index = nodeIndexes[tag];

                approxSolution[0] += displacements(3 * index) *
                                     mesh -> elems.basisFunctionsValues[j * mesh -> elems.noNodesPerElement + k];
                approxSolution[1] += displacements(3 * index + 1) *
                                     mesh -> elems.basisFunctionsValues[j * mesh -> elems.noNodesPerElement + k];
                approxSolution[2] += displacements(3 * index + 2) *
                                     mesh -> elems.basisFunctionsValues[j * mesh -> elems.noNodesPerElement + k];
            }

            for (int t = 0; t < 3; t++) {
                integral += mesh -> elems.weights[j] * (exact[t] - approxSolution[t]) * (exact[t] - approxSolution[t]);
            }

        }

        result += det * integral;
    }

    l2_error = sqrt(result);
}

void FEM3DVector::computeH1Error() {
    double result = l2_error * l2_error;

    for (int i = 0; i < mesh -> elems.elementTags.size(); i++) {
        // get tags of nodes in current element
        auto elementNodeTags = std::vector<std::size_t>(
                mesh -> elems.nodeTags.cbegin() + i * mesh -> elems.noNodesPerElement,
                mesh -> elems.nodeTags.cbegin() + (i + 1) * mesh -> elems.noNodesPerElement);
        double det = mesh -> global.determinants[i];

        double integral = 0;

        // quadrature
        for (int j = 0; j < mesh -> elems.noIntegrationPoints; j++) {
            // the exact solution gradient at the integration point
            Eigen::Matrix3d exact = params3d_ -> solution_gradient(mesh -> global.globalCoord[i * mesh -> elems.localCoord.size() + 3 * j],
                                                                   mesh -> global.globalCoord[i * mesh -> elems.localCoord.size() + 3 * j + 1],
                                                                   mesh -> global.globalCoord[i * mesh -> elems.localCoord.size() + 3 * j + 2]);

            // the approximate solution gradient at the integration point
            Eigen::Matrix3d approxGradient = Eigen::Matrix3d::Zero();

            Eigen::Vector3d referenceGrads, elementGrads;
            for (int k = 0; k < mesh -> elems.noNodesPerElement; k++) {
                referenceGrads << mesh -> elems.basisFunctionsGradients[j * 3 * mesh -> elems.noNodesPerElement + 3 * k],
                                  mesh -> elems.basisFunctionsGradients[j * 3 * mesh -> elems.noNodesPerElement + 3 * k + 1],
                                  mesh -> elems.basisFunctionsGradients[j * 3 * mesh -> elems.noNodesPerElement + 3 * k + 2];

                elementGrads = mesh -> global.inverse_jacobians[i] * referenceGrads;

                std::size_t tag = elementNodeTags[k];
                int index = nodeIndexes[tag];

                approxGradient.row(0) += displacements(3 * index) * elementGrads.transpose();
                approxGradient.row(1) += displacements(3 * index + 1) * elementGrads.transpose();
                approxGradient.row(2) += displacements(3 * index + 2) * elementGrads.transpose();
            }

            for (int u = 0; u < 3; u++) {
                for (int v = 0; v < 3; v++) {
                    integral += mesh -> elems.weights[j] * (exact(u, v) - approxGradient(u, v)) * (exact(u, v) - approxGradient(u, v));
                }
            }
        }

        result += det * integral;
    }

    h1_error = sqrt(result);

}

DirichletMap FEM3DVector::getDirichletBC() {
    return dirichlet_bc;
}