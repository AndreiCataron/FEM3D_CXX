#include "../include/Algoritm.hpp"

Algoritm::Algoritm(std::shared_ptr<LinearElasticity3D> const& le) : LE(le) {}

void Algoritm::iteration() {
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
}
