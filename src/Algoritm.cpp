#include "../include/Algoritm.hpp"
#include <iostream>

Algoritm::Algoritm(std::shared_ptr<LinearElasticity3D> const& le) : LE(le) {}

Algoritm::Algoritm(std::shared_ptr<LinearElasticity3D> const& le, std::vector<double>& plane) : LE(le), projectionPlane(plane) {}

void Algoritm::computeSolution() {
    LE -> computeStiffnessMatrixAndLoadVector();
    LE -> solveDirectProblem();
    LE -> computeL2Error();
    LE -> computeH1Error();
}

void Algoritm::unaccessibleDirichletPrep() {
    auto displacements = LE -> getDisplacements();
    auto freeBdryTags = LE -> getBdryFreeNodeTags();
    auto constrainedIdxs = LE -> getConstrainedNodes();
    auto tag2idx = LE -> getNodeIndexes();

    DirichletMap nextStepDir = {};

    // for each free boundary node (i.e. on Gamma_u at the previous step) assign the displacements from the previous step
    // as the next Dirichlet BC
    for (const auto& tag : freeBdryTags) {
        auto idx = tag2idx[tag];
        Eigen::VectorXd d = displacements(std::vector<int>{3 * idx, 3 * idx + 1, 3 * idx + 2});
        std::vector<double> bc(d.data(), d.data() + d.size());

        nextStepDir[tag] = bc;
    }

    // swap Neumann and Dirichlet boundaries
    auto par = LE -> getParamsPointer();
    par -> swapBC();

    // reset stresses and BC
    LE -> resetBoundaryConditions(true, true);
    // set prev-u|_Gamma_u as Dirichlet BC on Gamma_u
    LE -> setBoundaryConditions(nextStepDir);
    // Gamma_a nodes are free
    LE -> setBdryFreeNodeTags(accesibleBdryTags);
    // set initial stresses which will be used as an exact Neumann BC on Gamma_a
    LE -> setStresses(initialStresses);

    // delete stiff, load, displacements
    LE -> resetFEMData();
}


void Algoritm::unaccessibleNeumannPrep() {
    // reset stresses
    LE -> resetBoundaryConditions(true, false);
    // compute aproximated stresses from previous solution
    LE -> computeApproximateBoundaryGradients();
    LE -> computeApproximatedStresses();
    // reset BC
    LE -> resetBoundaryConditions(false, true);

    // swap Neumann and Dirichlet boundaries
    auto par = LE -> getParamsPointer();
    par -> swapBC();

    // set exact Dirichlet BC on Gamma_a
    LE -> setBoundaryConditions(initialDirichletCondition);
    // Gamma_u nodes are free
    LE -> setBdryFreeNodeTags(unaccessibleBdryTags);

    // delete stiff, load, displacements
    LE -> resetFEMData();
}

void Algoritm::iteration(int k) {
    if (k == 0) {
        // set exact Dirichlet BC Gamma_u and make Gamma_a Neumann boundary
        LE -> setBoundaryConditions();
        // get tags of nodes on Gamma_a
        accesibleBdryTags = LE -> getBdryFreeNodeTags();
        // compute boundary stresses from exact gradient
        LE -> computeIntegrationPointsStresses();
        initialStresses = LE -> getStresses();
        // get Dirichlet BC from Gamma_u and modify them to a random guess
        auto dirichletBC = LE -> getDirichletBC();
        for (const auto& [tag, bc] : dirichletBC) {
            dirichletBC[tag] = std::vector<double>{0.1, 0.1, 0.1};
        }
        // reset BC
        LE -> resetBoundaryConditions(false, true);
        // set exact Neumann BC on Gamma_a and guessed Dirichlet BC on Gamma_u
        LE -> setBoundaryConditions(dirichletBC);
        LE -> setBdryFreeNodeTags(accesibleBdryTags);
    }
    else if (k == 1) {
        // reset stresses
        LE -> resetBoundaryConditions(true, false);
        // compute aproximated stresses from previous solution
        LE -> computeApproximateBoundaryGradients();
        LE -> computeApproximatedStresses();
        // reset BC
        LE -> resetBoundaryConditions(false, true);

        // swap Neumann and Dirichlet boundaries
        auto par = LE -> getParamsPointer();
        par -> swapBC();

        // set exact Dirichlet BC on Gamma_a and make Gamma_u Neumann boundary
        LE -> setBoundaryConditions();
        // get tags of nodes on Gamma_u
        unaccessibleBdryTags = LE -> getBdryFreeNodeTags();
        // get exact Dirichlet BC from Gamma_a
        initialDirichletCondition = LE -> getDirichletBC();

        // delete stiff, load, displacements
        LE -> resetFEMData();
    }
    else if (k % 2 == 0) {
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

        auto d = LE -> getDisplacements();
        auto stiff = LE -> getStiffnessMatrix();
        auto load = LE -> getLoadVector();

        if (out) {
            LE -> outputData("/Users/andrei/CLionProjects/FEM/outputs/out" + std::to_string(i) + ".txt", true, projectionPlane);
        }
        double l2 = LE -> getL2Error(), h1 = LE -> getH1Error();
        std::cout << "Iteration " << i << '\n' << "L2 error: " << l2 << '\n' << "H1 Error: " << h1 << '\n' << "Time: " <<
                  std::chrono::duration <double, std::milli> (diff).count() << " ms, Displacements: " << d.size() << ' ' << "; Stiff: " << stiff.rows() << ' '  << stiff.cols() << "; load: " << load.size() << '\n';
    }
}