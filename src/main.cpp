#include <iostream>

#include <set>

#include "exprtk.hpp"

#include <gmsh.h>
#include <set>
#include <cmath>
#include <algorithm>

#include "../include/FEM3D.hpp"

void print_entities(std::vector<std::pair<int, int> > entities) {
    for(auto e : entities){
        int s = e.second;

        std::vector<std::size_t> tags;
        std::vector<double> coord, param;
        gmsh::model::mesh::getNodes(tags, coord, param, 2, s, true);

        for(auto i : tags){
            std::cout << i << " ";
        }
        std::cout << std::endl;

    }
}

int main(int argc, char **argv) {

    FEM3D::Params p = {2};
    FEM3D fem(p);

    auto start = std::chrono::steady_clock::now();

    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>   expression_t;
    typedef exprtk::parser<double>       parser_t;

    const std::string expression_string =
            "x + y == 4";

    double x, y;

    symbol_table_t symbol_table;
    symbol_table.add_variable("x",x);
    symbol_table.add_variable("y", y);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(expression_string,expression);

    x = double(1.00000000000000000001);
    y = double(3);
    double z = expression.value();
    std::cout << z << std::endl;

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;

    std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    gmsh::initialize(argc, argv);

    start = std::chrono::steady_clock::now();

    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();

    std::vector<std::pair<int, int>> faces;
    gmsh::model::getEntities(faces, 2);

    int tag = 1;
    for (auto face : faces) {
        gmsh::model::addPhysicalGroup(2, {face.first}, tag);
        gmsh::model::setPhysicalName(2, tag, "Face_" + std::to_string(tag));
        tag++;
    }

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.5);
    gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::mesh::generate(3);

    end = std::chrono::steady_clock::now();
    diff = end - start;

    std::cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;

    std::vector<std::pair<int, int>> boundary;
    std::vector<std::pair<int, int>> domain;

    domain.emplace_back(3, 1000);

    gmsh::model::getBoundary(domain, boundary, true, false, false);

    for(auto i : boundary) {
        std::cout << i.first << " " << i.second << std::endl;
    }

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t> > elementTags;
    std::vector<std::vector<std::size_t> > nodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);

    //std::cout << nodeTags[0].size() << std::endl;

    std::vector<double> coords;
    std::vector<double> parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags[0], coords, parametricCoords, 2, -1, true, true);

//    for(int i = 0; i < coords.size() / 3; i++){
//        std::cout << nodeTags[0][i] << " " << coords[3 * i] << " " << coords[3 * i + 1] << " " << coords[3 * i + 2] << " " <<
//             parametricCoords[3 * i] << " " << parametricCoords[3 * i + 1] << " " << parametricCoords[3 * i + 2] << std::endl;
//    }

    std::vector<std::pair<int, int> > entities;
    gmsh::model::getEntities(entities, 2);

//    for(auto e : entities){
//        int s = e.second;
//
//        std::vector<std::size_t> tags;
//        std::vector<double> coord, param;
//        gmsh::model::mesh::getNodes(tags, coord, param, 2, s, true);
//
//        for(auto i : coord){
//            std::cout << i << " ";
//        }
//        std::cout << std::endl;
//
//    }

    gmsh::write("/Users/andrei/CLionProjects/FEM/test01.msh");

//    std::set<std::string> args(argv, argv + argc);
//    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
