#include "../include/FEM3D.hpp"

#include <gmsh.h>
#include "../include/exprtk/exprtk.hpp"

FEM3D::FEM3D(const FEM3D::Params &params) : params_(params) {}

void FEM3D::setBoundaryConditions() {
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

void FEM3D::checkNodeSatisfiesBoundaryEquation(const std::size_t tag, double nx, double ny, double nz) {
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
            const std::string dirichlet_value = params_.g;
            parser.compile(dirichlet_value,expression);
            dirichlet_bc[tag] = expression.value();
        }
    }
}

std::unordered_map<std::size_t, double> FEM3D::getDirichletBC() {
    return dirichlet_bc;
}

FEM3D::Params FEM3D::getParams() {
    return params_;
}