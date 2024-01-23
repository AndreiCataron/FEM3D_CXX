#include "../include/Algoritm.hpp"
#include <iostream>

Algoritm::Algoritm(std::shared_ptr<LinearElasticity3D> const& le) : LE(le) {}

Algoritm::Algoritm(std::shared_ptr<LinearElasticity3D> const& le, std::vector<double>& plane) : LE(le), projectionPlane(plane) {}

void Algoritm::computeSolution() {
    LE -> computeStiffnessMatrixAndLoadVector();
    LE -> solveDirectProblem();
    LE -> computeL2Error();
    LE -> computeH1Error();
    LE -> computeApproximateBoundaryGradients();
}

void Algoritm::unaccessibleDirichletPrep() {
    auto displacements = LE -> getDisplacements();
    auto freeIdxs = LE -> getFreeNodes();
    auto constrainedIdxs = LE -> getConstrainedNodes();
    auto tag2idx = LE -> getNodeIndexes();

    DirichletMap nextStepDir = {};
    auto idx2tag = utils::inverseMap(tag2idx);

    for (const auto& idx : freeIdxs) {
        Eigen::VectorXd d = displacements(std::vector<int>{3 * idx, 3 * idx + 1, 3 * idx + 2});
        std::vector<double> bc(d.data(), d.data() + d.size());

        nextStepDir[idx2tag[idx]] = bc;
    }

    auto par = LE -> getParamsPointer();
    par -> swapBC();

    LE -> resetBoundaryConditions();

    LE -> setBoundaryConditions(nextStepDir);
}


void Algoritm::unaccessibleNeumannPrep() {
    LE -> resetBoundaryConditions();

    LE -> computeApproximateBoundaryGradients();
    LE -> computeApproximatedStresses();

    auto par = LE -> getParamsPointer();
    par -> swapBC();

    LE -> setBoundaryConditions();
}

void Algoritm::iteration(int k) {
    if (k == 0) {
        LE -> setBoundaryConditions();
    }
    else if (k % 2 == 1) {
        unaccessibleDirichletPrep();
    }
    else {
        unaccessibleNeumannPrep();
    }
    computeSolution();
}

void Algoritm::iterations(int k, double tol, bool out) {
    for (int i = 0; i < k; i++) {
        auto start = std::chrono::steady_clock::now();
        iteration(i);
        auto end = std::chrono::steady_clock::now();
        auto diff = end - start;

        if (out) {
            LE -> outputData("/Users/andrei/CLionProjects/FEM/outputs/out" + std::to_string(i) + ".txt", true, projectionPlane);
        }
        double l2 = LE -> getL2Error(), h1 = LE -> getH1Error();
        std::cout << "Iteration " << i << '\n' << "L2 error: " << l2 << '\n' << "H1 Error: " << h1 << '\n' << "Time: " <<
                     std::chrono::duration <double, std::milli> (diff).count() << " ms" << '\n';
    }
}