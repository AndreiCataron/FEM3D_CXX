#include <iostream>
#include <cassert>
#include "../include/LinearElasticity3D.hpp"

int main() {
    FEM3DVector::ParamsVector par = {
            0,
            "x == 0"
    };

    LinearElasticity3D fem(par);

    int check = fem.checkNodeSatisfiesBoundaryEquation(1, 0, 0.5, 0.3);

    assert(check == 1);

    check = fem.checkNodeSatisfiesBoundaryEquation(1, 0.2, 0.3, 0);

    assert(check == 0);


    return 0;
}