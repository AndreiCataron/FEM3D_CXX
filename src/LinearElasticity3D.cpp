#include "../include/LinearElasticity3D.hpp"
#include <iostream>
#include <set>
#include <gmsh.h>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <thread>
#include <future>

LinearElasticity3D::LinearElasticity3D(std::shared_ptr<ParamsLE> const& params) : FEM3DVector(params), paramsLE_(params) {
    computeIntegrationPointsStresses();
}

LinearElasticity3D::LinearElasticity3D(std::shared_ptr<ParamsLE> const& params, std::shared_ptr<Mesh> const& msh) : FEM3DVector(params, msh), paramsLE_(params) {
    computeIntegrationPointsStresses();
}

void LinearElasticity3D::resetBoundaryConditions() noexcept {
    dirichlet_bc.clear();
    integrationPointsStresses.clear();
    freeNodes.clear();
    constrainedNodes.clear();
    neumannBoundaryTriangles.clear();
}

void LinearElasticity3D::computeIntegrationPointsStresses() noexcept {
    for (int i = 0; i < mesh -> global.trianglesGlobalCoord.size() / 3; i++) {
        Eigen::Matrix3d strain = 0.5 * (paramsLE_ -> solution_gradient(mesh -> global.trianglesGlobalCoord[3 * i],
                                                                       mesh -> global.trianglesGlobalCoord[3 * i + 1],
                                                                       mesh -> global.trianglesGlobalCoord[3 * i + 2]) +
                                 (paramsLE_ -> solution_gradient(mesh -> global.trianglesGlobalCoord[3 * i],
                                                                 mesh -> global.trianglesGlobalCoord[3 * i + 1],
                                                                 mesh -> global.trianglesGlobalCoord[3 * i + 2])).transpose());

        // compute stress tensor
        Eigen::Matrix3d stress = (paramsLE_ -> E) / (1 + paramsLE_ -> nu) * (strain + (paramsLE_ -> nu) / (1 - 2 * paramsLE_ -> nu) * strain.trace() * Eigen::Matrix3d::Identity());

        integrationPointsStresses.emplace_back(stress);
    }
}

Eigen::Vector3d LinearElasticity3D::h(const int nodeId) {
    return integrationPointsStresses[nodeId] * mesh -> global.normals[nodeId];
}

void LinearElasticity3D::computeApproximatedStresses() noexcept {
    approximatedStresses.reserve(mesh -> global.trianglesGlobalCoord.size() / 3);
    for (const auto& grad : approx_grads) {
        Eigen::Matrix3d strain = 0.5 * (grad + grad.transpose());
        Eigen::Matrix3d stress = (paramsLE_ -> E) / (1 + paramsLE_ -> nu) * (strain + (paramsLE_ -> nu) / (1 - 2 * paramsLE_ -> nu) * strain.trace() * Eigen::Matrix3d::Identity());

        approximatedStresses.emplace_back(stress);
    }
}

void LinearElasticity3D::stiffnessIterations(int start, int end) {
    for (int i = start; i < end; i++) {
        std::vector<std::size_t> elementNodeTags = std::vector<std::size_t>(
                mesh -> elems.nodeTags.cbegin() + i * mesh -> elems.noNodesPerElement,
                mesh -> elems.nodeTags.cbegin() + (i + 1) * mesh -> elems.noNodesPerElement);

        // compute element stiffness matrix
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3 * mesh -> elems.noNodesPerElement, 3 * mesh -> elems.noNodesPerElement);

        for (int j = 0; j < mesh -> elems.noIntegrationPoints; j++) {
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 3 * mesh -> elems.noNodesPerElement);

            for (int k = 0; k < mesh -> elems.noNodesPerElement; k++) {
                Eigen::Vector3d referenceGrads, elementGrads;
                referenceGrads << mesh -> elems.basisFunctionsGradients[j * 3 * mesh -> elems.noNodesPerElement + 3 * k],
                        mesh -> elems.basisFunctionsGradients[j * 3 * mesh -> elems.noNodesPerElement + 3 * k + 1],
                        mesh -> elems.basisFunctionsGradients[j * 3 * mesh -> elems.noNodesPerElement + 3 * k + 2];

                elementGrads = mesh -> global.inverse_jacobians[i] * referenceGrads;

                B(0, 3 * k) = elementGrads(0);
                B(1, 3 * k + 1) = elementGrads(1);
                B(2, 3 * k + 2) = elementGrads(2);
                B(3, 3 * k) = elementGrads(1);
                B(3, 3 * k + 1) = elementGrads(0);
                B(4, 3 * k + 1) = elementGrads(2);
                B(4, 3 * k + 2) = elementGrads(1);
                B(5, 3 * k) = elementGrads(2);
                B(5, 3 * k + 2) = elementGrads(0);
            }

            K = K + mesh -> elems.weights[j] * B.transpose() * D * B;
        }

        double det = mesh -> global.determinants[i];

        // find stiffness and load indexes for the displacements in the element
        std::vector<int> elementNodeIndexes;
        elementNodeIndexes.reserve(3 * mesh -> elems.noNodesPerElement);

        for (const auto &tag: elementNodeTags) {
            int index = nodeIndexes[tag];

            elementNodeIndexes.emplace_back(3 * index);
            elementNodeIndexes.emplace_back(3 * index + 1);
            elementNodeIndexes.emplace_back(3 * index + 2);
        }

        // the vector f^k from Larson
        std::vector<double> fk;
        fk.reserve(3 * mesh -> elems.noNodesPerElement);
        for (const auto &tag: elementNodeTags) {
            // std::tuple<double, double, double> nodeCoord = mesh -> elems.node_coordinates[tag];

            for (auto component: paramsLE_->f(std::get<0>(mesh -> global.node_coordinates[tag]),
                                              std::get<1>(mesh -> global.node_coordinates[tag]),
                                              std::get<2>(mesh -> global.node_coordinates[tag]))) {
                fk.emplace_back(component);
            }
        }

        double *ptr = &fk[0];
        Eigen::Map<Eigen::VectorXd> fkConverted(ptr, int(fk.size()));

        mtx_stiff.lock();
        // add triplets to tripletList
        // repeated pairs of indexes are summed up when initializing the sparse stiffness matrix
        for (int k = 0; k < 3 * mesh -> elems.noNodesPerElement; k++) {
            for (int j = 0; j < 3 * mesh -> elems.noNodesPerElement; j++) {
                tripletList.emplace_back(elementNodeIndexes[k], elementNodeIndexes[j], det * K(k, j));
            }
        }
        load_vector(elementNodeIndexes) = load_vector(elementNodeIndexes) + det * Mk * fkConverted;
        mtx_stiff.unlock();
    }
}

void LinearElasticity3D::neumannIterations(int start, int end) {
    for (int i = start; i < end; i++) {
        std::size_t faceTag = mesh -> bdry.boundaryFacesTags[i];
        if (neumannBoundaryTriangles.find(faceTag) != neumannBoundaryTriangles.cend()) {
            std::vector<std::size_t> faceNodeTags = std::vector<std::size_t>(
                    mesh -> bdry.boundaryFacesNodes.cbegin() + i * mesh -> bdry.noNodesPerTriangle,
                    mesh -> bdry.boundaryFacesNodes.cbegin() + (i + 1) * mesh -> bdry.noNodesPerTriangle);

            // compute integral of h * phi for all basis functions
            for (int k = 0; k < mesh -> bdry.noNodesPerTriangle; k++) {
                Eigen::Vector3d integral = Eigen::Vector3d::Zero();

                for (int j = 0; j < mesh -> bdry.triangleNoIntegrationPoints; j++) {
                      integral += mesh -> bdry.triangleWeights[j] *
                                  mesh -> bdry.triangleBasisFunctionsValues[j * mesh -> bdry.noNodesPerTriangle + k] *
                                  h(i * mesh -> bdry.triangleNoIntegrationPoints + j);
                }

                integral *= mesh -> global.trianglesDeterminants[i];

                int index = nodeIndexes[faceNodeTags[k]];
                std::vector<int> nodeIndexes = {3 * index, 3 * index + 1, 3 * index + 2};

                mtx_neu.lock();
                load_vector(nodeIndexes) = load_vector(nodeIndexes) + integral;
                mtx_neu.unlock();
            }
        }
    }
}

void LinearElasticity3D::computeStiffnessMatrixAndLoadVector() {
    // initialize D matrix
    D = Eigen::MatrixXd(6, 6);
    double lam = paramsLE_ -> lambda;
    double mu = paramsLE_ -> mu;
    D << lam + 2 * mu, lam, lam, 0, 0, 0,
         lam, lam + 2 * mu, lam, 0, 0, 0,
         lam, lam, lam + 2 * mu, 0, 0, 0,
         0, 0, 0, mu, 0, 0,
         0, 0, 0, 0, mu, 0,
         0, 0, 0, 0, 0, mu;

    // number of nodes in the mesh
    // first, delete duplicate nodes (which appear because they are shared by multiple surfaces)
    std::vector<std::size_t> tags;
    std::vector<double> co, pc;
    gmsh::model::mesh::getNodes(tags, co, pc, -1, -1, true, false);
    utils::deleteDuplicatesFromVector(tags);
    int noNodes = int(tags.size());

    // computation of the matrix M^k from Larson page 270

    Mk = Eigen::MatrixXd(3 * mesh -> elems.noNodesPerElement, 3 * mesh -> elems.noNodesPerElement);

    for (int i = 0; i < mesh -> elems.noIntegrationPoints; i++) {
        Eigen::MatrixXd Mktemp = Eigen::MatrixXd::Zero(3 * mesh -> elems.noNodesPerElement, 3);
        for (int j = 0; j < mesh -> elems.noNodesPerElement; j++) {
            Eigen::Matrix3d basisFunctionMatrix = Eigen::Matrix3d::Zero();
            basisFunctionMatrix(0, 0) = mesh -> elems.basisFunctionsValues[i * mesh -> elems.noNodesPerElement + j];
            basisFunctionMatrix(1, 1) = mesh -> elems.basisFunctionsValues[i * mesh -> elems.noNodesPerElement + j];
            basisFunctionMatrix(2, 2) = mesh -> elems.basisFunctionsValues[i * mesh -> elems.noNodesPerElement + j];
            Mktemp.block<3, 3>(3 * j, 0) = basisFunctionMatrix;
        }
        Mk = Mk + mesh -> elems.weights[i] * (Mktemp * Mktemp.transpose());
    }

    // initialize load vector
    load_vector = Eigen::VectorXd::Zero(3 * noNodes);

    // initialize stiffness matrix as Eigen Sparse Matrix
    tripletList.reserve(mesh -> elems.elementTags.size() * 3 * mesh -> elems.noNodesPerElement * 3 * mesh -> elems.noNodesPerElement);
    stiffness_matrix.resize(3 * noNodes, 3 * noNodes);

    // assemble tripletList and load vector by looping through all elements and computing the element stiffness matrix and load vector
    const int num_threads = 4;

    //std::vector<std::jthread> threads;
    std::vector<std::thread> threads;
    int chunk_size = int(mesh -> elems.elementTags.size()) / num_threads;

    for (int t = 0; t < num_threads; t++) {
        int startIdx = t * chunk_size;
        int endIdx = (t + 1) * chunk_size;
        if (t == num_threads - 1) {
            endIdx = int(mesh -> elems.elementTags.size());
        }

        threads.emplace_back([startIdx, endIdx, this]()
                             {stiffnessIterations(startIdx, endIdx);} );
    }

    //neumannIterations(0, int(mesh -> bdry.boundaryFacesTags.size()));
    // auto a = std::async(std::launch::async, [this](){neumannIterations(0, int(mesh -> bdry.boundaryFacesTags.size())); } );

    for (auto &thread : threads) {
        thread.join();
    }

    std::vector<std::thread> neuThreads;
    chunk_size = int(mesh -> bdry.boundaryFacesTags.size()) / num_threads;

    for (int t = 0; t < num_threads; t++) {
        int startIdx = t * chunk_size;
        int endIdx = (t + 1) * chunk_size;
        if (t == num_threads - 1) {
            endIdx = int(mesh -> bdry.boundaryFacesTags.size());
        }

        neuThreads.emplace_back([startIdx, endIdx, this]()
                             { neumannIterations(startIdx, endIdx);} );
    }

    for (auto &thread : neuThreads) {
        thread.join();
    }

    stiffness_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

}

void LinearElasticity3D::solveDirectProblem() {
    displacements = Eigen::VectorXd(3 * nodeIndexes.size());

    // assemble vector of constrained values
    Eigen::VectorXd constrainedValues(3 * constrainedNodes.size());

    std::unordered_map<int, std::size_t> reverseNodeIndexes = utils::inverseMap(nodeIndexes);

    std::vector<double> tempDirBC;

    for (const auto& idx : constrainedNodes) {
        std::size_t tag = reverseNodeIndexes[idx];

        tempDirBC = dirichlet_bc[tag];

        constrainedValues(3 * idx) = tempDirBC[0];
        constrainedValues(3 * idx + 1) = tempDirBC[1];
        constrainedValues(3 * idx + 2) = tempDirBC[2];
    }

//    std::vector<int> freeIndexes, constrainedIndexes;
//    freeIndexes.reserve(3 * freeNodes.size());
//    constrainedIndexes.reserve(3 * constrainedNodes.size());
//    for (const auto& idx : freeNodes) {
//        freeIndexes.emplace_back(3 * idx);
//        freeIndexes.emplace_back(3 * idx + 1);
//        freeIndexes.emplace_back(3 * idx + 2);
//    }
//    for (const auto& idx : constrainedNodes) {
//        constrainedIndexes.emplace_back(3 * idx);
//        constrainedIndexes.emplace_back(3 * idx + 1);
//        constrainedIndexes.emplace_back(3 * idx + 2);
//    }

    displacements.head(3 * constrainedNodes.size()) = constrainedValues;
    //displacements(constrainedIndexes) = constrainedValues;

     load_vector = (load_vector.tail(3 * freeNodes.size()) -
            stiffness_matrix.bottomRows(3 * freeNodes.size()).leftCols(3 * constrainedNodes.size()) * constrainedValues).eval();

    stiffness_matrix = stiffness_matrix.bottomRightCorner(3 * freeNodes.size(), 3 * freeNodes.size()).eval();

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;

    auto start = std::chrono::steady_clock::now();

    displacements.tail(3 * freeNodes.size()) = solver.compute(stiffness_matrix).solve(load_vector);

    auto end = std::chrono::steady_clock::now();

    auto diff = end - start;

    std::cout << "SOLVING TIME:" << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;


    try {
        if (solver.info() != Eigen::Success) {
                // solving failed
                throw 1;
        }
    }
    catch (int err) {
        std::cout << "solving the sparse system failed";
        return;
    }
}