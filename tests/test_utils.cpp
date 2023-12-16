#include "../include/utils.hpp"
#include <iostream>

int main () {
    std::tuple<double, double, double> punct{-4, -9, -5};
    std::vector<double> plan{5, 11, 4, -23};
    std::tuple<double, double, double> proj = utils::projectionOnPlane(punct, plan);
    std::cout << std::get<0>(proj) << ' ' << std::get<1>(proj) << ' ' << std::get<2>(proj);
    return 0;
}