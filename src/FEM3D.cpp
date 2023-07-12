#include "../include/FEM3D.hpp"

#include <gmsh.h>
#include "../include/exprtk/exprtk.hpp"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

FEM3D::FEM3D(const FEM3D::Params &params) : params_(params) {}

double FEM3D::parseExpression(std::string exp, double xx, double yy, double zz) {
    double x, y, z;

    symbol_table_t symbol_table;
    symbol_table.add_variable("x",x);
    symbol_table.add_variable("y", y);
    symbol_table.add_variable("z", z);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    const std::string e = exp;

    parser.compile(e, expression);

    x = xx; y = yy; z = zz;

    return expression.value();
}

int FEM3D::checkNodeSatisfiesBoundaryEquation(const std::size_t tag, double nx, double ny, double nz) {

    if(params_.dirichlet_bc.empty() != 1) {
        int check_condition = int(parseExpression(params_.dirichlet_bc, nx, ny, nz));

        if(check_condition == 1) {
            return 1;
        }
    }
    return 0;
}

void FEM3D::setupMesh() {
    gmsh::model::mesh::setOrder(int(params_.element_order));
}

FEM3D::Params FEM3D::getParams() {
    return params_;
}