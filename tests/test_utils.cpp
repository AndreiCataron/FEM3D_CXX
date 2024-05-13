#include "../include/utils.hpp"
#include <iostream>

int main () {
    std::tuple<double, double, double> punct{-4, -9, -5};
    std::vector<double> plan{5, 11, 4, -23};
    std::tuple<double, double, double> proj = utils::projectionOnPlane(punct, plan);
    std::cout << std::get<0>(proj) << ' ' << std::get<1>(proj) << ' ' << std::get<2>(proj);
    int s = 1;
    for (int i = 0; i < 1000000; i++) {
        s *= i;
    }
    return 0;
}