#include "../include/utils.hpp"
#include <cmath>

int utils::binom(int n, int k) {
    double coefficient = std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1));
    return int(coefficient);
}


