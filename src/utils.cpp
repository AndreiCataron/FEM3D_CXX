#include "../include/utils.hpp"
#include <cmath>

int utils::binom(int n, int k) {
    double coefficient = std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1));
    return int(coefficient);
}

void utils::checkParamsLE(ParamsLE &params) {
    // if E and nu are given, compute lambda and mu
    if (params.lambda < 0 && params.mu < 0 && params.nu >= -1 && params.nu <= 0.5 && params.E > 0) {
        params.lambda = params.E / ((1 + params.nu) * (1 - 2 * params.nu));
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
        // set E to aluminimum E
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