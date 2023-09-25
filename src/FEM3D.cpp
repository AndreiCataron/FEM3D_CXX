#include "../include/FEM3D.hpp"

#include <gmsh.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <utility>
#include "../include/exprtk/exprtk.hpp"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

FEM3D::FEM3D(FEM3D::Params params) : params_(std::move(params)) {
    symbol_table.add_variable("x", _x);
    symbol_table.add_variable("y", _y);
    symbol_table.add_variable("z", _z);
    symbol_table.add_constants();

    expression.register_symbol_table(symbol_table);
}

double FEM3D::parseExpression(const std::string &exp, double xx, double yy, double zz) {
    parser.compile(exp, expression);

    _x = xx; _y = yy; _z = zz;

    return expression.value();
}

int FEM3D::checkNodeSatisfiesBoundaryEquation(double nx, double ny, double nz) {

    if (params_.dirichlet_bc.empty() != 1) {
        int check_condition = int(parseExpression(params_.dirichlet_bc, nx, ny, nz));

        if (check_condition == 1) {
            return 1;
        }
    }
    return 0;
}

void FEM3D::setupMesh() {
    gmsh::model::mesh::generate(3);

    gmsh::model::mesh::setOrder(int(params_.element_order));

    getNodesCoordinates();
}

void FEM3D::indexFreeNodes() {
    // get nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, parametricCoord;

    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, true, false);

    for (auto tag : nodeTags) {
        // check node with tag is not constrained i.e. already in nodeIndexes
        if (nodeIndexes.find(tag) == nodeIndexes.end()) {
            freeNodes.emplace_back(nodeIndexes.size());
            nodeIndexes[tag] = int(nodeIndexes.size());
        }
    }
}

// getters

FEM3D::Params FEM3D::getParams() {
    return params_;
}

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

