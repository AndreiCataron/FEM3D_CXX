#include "../include/FEM3Dvector.hpp"

#include <iostream>
#include <gmsh.h>
#include "../include/exprtk/exprtk.hpp"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double>   expression_t;
typedef exprtk::parser<double>       parser_t;

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
            std::size_t tag = tags[i];
            int bc = checkNodeSatisfiesBoundaryEquation(tag, coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);

            if(bc == 1){
                std::vector<double> prescribed_condition;
                for(const auto &cond : params3d_.g){
                    prescribed_condition.emplace_back(parseExpression(cond, coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]));
                }
                dirichlet_bc[tag] = prescribed_condition;
            }
        }
    }
}

std::unordered_map<std::size_t, std::vector<double> > FEM3DVector::getDirichletBC() {
    return dirichlet_bc;
}