#include <iostream>
#include <cassert>
#include "../include/LinearElasticity3D.hpp"

int main() {
    LinearElasticity3D::ParamsLE par = {
            0, // h
            "x == 0", // dirichlet BC
            3, // quadrature precision
            2, // order of lagrage polynomials
            {"2 * x", "3 * y + 1", "x + z"}, // g
            {"0", "1", "0"}, // f
            1, // lambda
            2, // mu
    };

    LinearElasticity3D fem(par);

    int check = fem.checkNodeSatisfiesBoundaryEquation(0, 0.5, 0.3);

    assert(check == 1);

    check = fem.checkNodeSatisfiesBoundaryEquation(0.2, 0.3, 0);

    assert(check == 0);


    return 0;
}