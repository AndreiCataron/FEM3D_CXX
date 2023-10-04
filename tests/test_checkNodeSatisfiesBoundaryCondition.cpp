#include <iostream>
#include <cassert>
#include "../include/LinearElasticity3D.hpp"
#include "../include/utils.hpp"

int main() {
    LinearElasticity3D::ParamsLE par2 = {
            0.1, // h
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

    utils::checkParamsLE(par2);

    LinearElasticity3D fem(par2);

    int check = fem.checkNodeSatisfiesBoundaryEquation(0, 0, 0.3);

    assert(check == 1);

    check = fem.checkNodeSatisfiesBoundaryEquation(0.2, 0.3, 0);

    assert(check == 1);


    return 0;
}