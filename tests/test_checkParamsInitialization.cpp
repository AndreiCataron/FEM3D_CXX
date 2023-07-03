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

    FEM3D::Params check = fem.getParams();

    assert(check.h == 0);
    assert(check.dirichlet_bc == "x == 0");
    assert(check.g == "x + y + 2 * z");

    return 0;
}