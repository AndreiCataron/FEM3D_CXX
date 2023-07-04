#include "../include/FEM3D.hpp"

#include <gmsh.h>
#include "../include/exprtk/exprtk.hpp"

FEM3D::FEM3D(const FEM3D::Params &params) : params_(params) {}

bool FEM3D::checkNodeSatisfiesBoundaryEquation(const std::size_t tag, double nx, double ny, double nz) {
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>   expression_t;
    typedef exprtk::parser<double>       parser_t;

    double x, y, z;

    symbol_table_t symbol_table;
    symbol_table.add_variable("x",x);
    symbol_table.add_variable("y", y);
    symbol_table.add_variable("z", z);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;

    if(params_.dirichlet_bc.empty() != 1) {
        const std::string dirichlet_condition = params_.dirichlet_bc;
        parser.compile(dirichlet_condition,expression);

        x = nx; y = ny; z = nz;

        if(expression.value() == 1) {
            return true;
//            const std::string dirichlet_value = params_.g;
//            parser.compile(dirichlet_value,expression);
//            dirichlet_bc[tag] = expression.value();
        }
        else {
            return false;
        }
    }
    return false;
}

std::unordered_map<std::size_t, double> FEM3D::getDirichletBC() {
    return dirichlet_bc;
}

FEM3D::Params FEM3D::getParams() {
    return params_;
}