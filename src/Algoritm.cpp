#include "../include/Algoritm.hpp"

Algoritm::Algoritm(std::shared_ptr<LinearElasticity3D> const& le) : LE(le) {}

void Algoritm::iteration() {
    auto displacements = LE -> getDisplacements();
    auto freeIdx = LE -> getFreeNodes();
    auto constrainedIdx = LE -> getConstrainedNodes();
    auto tag2idx = LE -> getNodeIndexes();
}
