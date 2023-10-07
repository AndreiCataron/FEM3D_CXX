#include <iostream>
#include "../include/LinearElasticity3D.hpp"
#include <gmsh.h>
#include <unordered_map>

int main(int argc, char **argv) {
    ParamsLE par = {
            0.1, // h
            10,
            1,
            "0 == 0",
            1,
            1,
            {"-1.9831", "-1.9831", "50.016"},
            {"x * z + 1 / 56 * x ^ 2", "y * z + 1 / 56 * y ^ 2", "- z ^ 2 + 1 / 56 * z ^ 2"},
            56,
            36,
            -1,
            -1
    };

    Mesh msh(argc, argv);

    LinearElasticity3D fem(par);

    gmsh::initialize(argc, argv);

    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.1);
    gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::mesh::generate(3);

    fem.setBoundaryConditions();

    std::unordered_map<std::size_t, std::vector<double> > map = fem.getDirichletBC();

    std::cout << "start";
    for (const auto& n : map){
        auto tag = n.first;
        std::vector<double> coord, paramcoord;
        int dim, t;
        gmsh::model::mesh::getNode(tag, coord, paramcoord, dim, t);
        for (auto i : coord) std::cout << i << ' ';
        std::cout << '\n';
//        assert(abs(n.second[0] - 2 * coord[0]) < 1e-8);
//        assert(abs(n.second[1] - 3 * coord[1] - 1) < 1e-8);
//        assert(abs(n.second[2] - coord[0] - coord[2]) < 1e-8);
    }

    // check node indexes

    std::unordered_map<std::size_t, int> ni = fem.getNodeIndexes();
    for (const auto& n : ni){
        std::cout << n.first << " " << n.second << '\n';
    }

    std::vector<int> con, free;
    con = fem.getConstrainedNodes();
    free = fem.getFreeNodes();

    for (const auto& i : con) std::cout << i << " ";

    std::cout << '\n';

    for (const auto& i : free) std::cout << i << " ";

    gmsh::finalize();

    return 0;
}