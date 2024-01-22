#include "../include/utils.hpp"
#include <cmath>
#include <iostream>
#include <tuple>

int utils::binomialCoefficient(int n, int k) {
    double coefficient = std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1));
    return int(coefficient);
}

void utils::checkParamsLE(ParamsLE &params) {
    // if E and nu are given, compute lambda and mu
    if (params.lambda < 0 && params.mu < 0 && params.nu >= -1 && params.nu <= 0.5 && params.E > 0) {
        params.lambda = params.E * params.nu / ((1 + params.nu) * (1 - 2 * params.nu));
        params.mu = params.E / (2 * (1 + params.nu));
    }
    else {
        if (params.lambda < 0) {
            // set lambda to aluminium lambda
            params.lambda = 39.2;
        }

        if (params.mu < 0) {
            // set mu to aluminium mu
            params.mu = 25;
        }
    }

    if ((params.nu >= -1 && params.nu <= 0.5) == 0) {
        // set nu to aluminium nu
        params.nu = 0.3;
    }

    if (params.E < 0) {
        // set E to aluminium E
        params.E = 68;
    }
}

void utils::generateCombinations(std::vector<std::vector<int> > &res, int n, int k) {
    std::vector<bool> v(n);
    std::fill(v.end() - k, v.end(), true);

    do {
        std::vector<int> temp;
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                temp.emplace_back(i);
            }
        }
        res.emplace_back(temp);
    } while (std::next_permutation(v.begin(), v.end()));
}

std::vector<double> utils::midpoint(std::vector<double> &p1, std::vector<double> &p2) {
    try {
        if (p1.size() == p2.size()) {
            std::vector<double> mid;
            mid.reserve(p1.size());
            for (int i = 0; i < p1.size(); i++) {
                mid.emplace_back((p1[i] + p2[i]) / 2);
            }
            return mid;
        }
        else {
            throw std::runtime_error("The points must have the same length");
        }
    }
    catch(std::exception const &e) {
        std::cout << "Runtime error: " << e.what() << '\n';
    }
    return {};
}

std::tuple<double, double, double> utils::projectionOnPlane(std::tuple<double, double, double> &point, const std::vector<double> &plane) {
    try {
        if (plane.size() == 4) {

            double coef = (plane[0] * std::get<0>(point) + plane[1] * std::get<1>(point) + plane[2] * std::get<2>(point) + plane[3]) /
                    (plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2]);


            return std::tuple<double, double, double>{std::get<0>(point) - coef * plane[0], std::get<1>(point) - coef * plane[1],
                                                      std::get<2>(point) - coef * plane[2]};

        }
        else {
            throw std::runtime_error("The vector does not correctly discribe the equation of a plane");
        }
    }
    catch(std::exception const &e) {
        std::cout << "Runtime error: " << e.what() << '\n';
    }
    return {};
}