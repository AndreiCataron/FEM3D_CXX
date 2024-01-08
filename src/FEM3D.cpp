#include "../include/FEM3D.hpp"

#include <gmsh.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

FEM3D::FEM3D(std::shared_ptr<Params> const& params) : params_(params){
    auto msh = std::make_shared<Mesh>();
    mesh = msh;
}

FEM3D::FEM3D(std::shared_ptr<Params> const& params, std::shared_ptr<Mesh> const& msh) : params_(params), mesh(msh) {
    // used for parsing logical expressions i.e. boundary conditions
    symbol_table.add_variable("x", _x);
    symbol_table.add_variable("y", _y);
    symbol_table.add_variable("z", _z);
    symbol_table.add_constants();

    expression.register_symbol_table(symbol_table);
}

double FEM3D::parseExpression(const std::string& exp, double xx, double yy, double zz) {
    parser.compile(exp, expression);

    _x = xx; _y = yy; _z = zz;

    return expression.value();
}

int FEM3D::checkNodeSatisfiesBoundaryEquation(double nx, double ny, double nz) {

    if (params_ -> dirichlet_bc.empty() != 1) {
        int check_dir = int(parseExpression(params_ -> dirichlet_bc, nx, ny, nz));

        if (check_dir == 1) {
            return 1;
        }
    }
    if (params_ -> neumann_bc.empty() != 1) {
        int check_neu = int(parseExpression(params_ -> neumann_bc, nx, ny, nz));

        if (check_neu == 1) {
            return 2;
        }
    }

    return 0;
}

void FEM3D::setNeumannBoundaryConditions() noexcept {
    for (const auto& b : mesh -> bdry.boundary) {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t> > elementTags, nodeTags;
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, b.first, b.second);

        mesh -> bdry.boundaryFacesTags.insert(mesh -> bdry.boundaryFacesTags.cend(), elementTags[0].cbegin(), elementTags[0].cend());
        mesh -> bdry.boundaryFacesNodes.insert(mesh -> bdry.boundaryFacesNodes.cend(), nodeTags[0].cbegin(), nodeTags[0].cend());

        for (std::size_t i = 0; i < elementTags[0].size(); i++) {
            std::size_t triangleTag = elementTags[0][i];
            std::vector<std::size_t> faceNodeTags = std::vector<std::size_t>(nodeTags[0].cbegin() + int(i) * mesh -> bdry.noNodesPerTriangle,
                                                                             nodeTags[0].cbegin() + int(i + 1) * mesh -> bdry.noNodesPerTriangle);

            bool ok = true;
            for (const auto& tag : faceNodeTags) {
                CoordTuple coord = mesh -> global.node_coordinates[tag];
                int bc = checkNodeSatisfiesBoundaryEquation(std::get<0>(coord), std::get<1>(coord), std::get<2>(coord));

                if (bc != 2) {
                    ok = false;
                    break;
                }
                else continue;
            }

            if (ok) {
                neumannBoundaryTriangles[triangleTag] = b.second;
            }
        }
    }
}

void FEM3D::setBoundaryConditions() noexcept {
    setNeumannBoundaryConditions();
    setDirichletBoundaryConditions();
    indexConstrainedNodes();
    indexFreeNodes();
}

void FEM3D::indexFreeNodes() noexcept {
    // get nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, parametricCoord;

    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, true, false);

    for (const auto& tag : nodeTags) {
        // check node with tag is not constrained i.e. already in nodeIndexes
        if (nodeIndexes.find(tag) == nodeIndexes.cend()) {
            freeNodes.emplace_back(nodeIndexes.size());
            nodeIndexes[tag] = int(nodeIndexes.size());
        }
    }
}

// getters

std::unordered_map<std::size_t, int> FEM3D::getNodeIndexes() {
    return nodeIndexes;
}

std::vector<int> FEM3D::getConstrainedNodes() {
    return constrainedNodes;
}

std::vector<int> FEM3D::getFreeNodes() {
    return freeNodes;
}

Eigen::SparseMatrix<double> FEM3D::getStiffnessMatrix() {
    return stiffness_matrix;
}

Eigen::VectorXd FEM3D::getLoadVector() {
    return load_vector;
}

Eigen::VectorXd FEM3D::getDisplacements() {
    return displacements;
}

double FEM3D::getL2Error() const {
    return l2_error;
}

double FEM3D::getH1Error() const {
    return h1_error;
}
