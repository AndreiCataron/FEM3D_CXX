#include <iostream>
#include <cassert>
#include "../include/FEM3D.hpp"

int main() {
    FEM3D::Params par = {
            0,
            "x == 0",
            "x + y + 2 * z"
    };

    FEM3D fem(par);

    fem.checkNodeSatisfiesBoundaryEquation(1, 0, 0.5, 0.3);

    std::unordered_map<std::size_t, double> bc = fem.getDirichletBC();

    for(auto i : bc){
        std::cout << i.first << " " << i.second << std::endl;
    }

    assert(bc[1] == 1.1);

    return 0;
}