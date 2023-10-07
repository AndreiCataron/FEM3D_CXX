#include "../include/Mesh.hpp"
#include <gmsh.h>
#include "../include/params.hpp"
#include "../include/utils.hpp"
#include <iostream>

Mesh::Mesh() {
    if (gmsh::isInitialized() == 0) {
        gmsh::initialize();
    }
    params = nullptr;
}

Mesh::Mesh(int argc, char **argv) {
    if (gmsh::isInitialized() == 0) {
        gmsh::initialize(argc, argv);
    }
    params = nullptr;
}

Mesh::Mesh(int argc, char **argv, const std::string& meshType) : Mesh(argc, argv) {
    if (meshType.compare("unit_cube") == 1) {
        cubeMesh();
    }
}

Mesh::Mesh(int argc, char **argv, std::shared_ptr<Params> const &par) : Mesh(argc, argv) {
    params = par;
}

//Mesh::Mesh(const Mesh &other) : elems(other.elems) {
//    params = new Params;
//    *params = *other.params;
//}
//
//Mesh& Mesh::operator=(const Mesh &other) {
//    if (this == &other) {
//        return *this;
//    }
//
//    delete params;
//
//    params = new Params;
//    *params = *other.params;
//
//    elems = other.elems;
//
//    return *this;
//
//}

Mesh::~Mesh() {
    if (gmsh::isInitialized() == 1) {
        gmsh::finalize();
    }
}

void Mesh::getNodesCoordinates() {
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, true, false);

    for (int i = 0; i < nodeTags.size(); i++) {
        elems.node_coordinates[nodeTags[i]] = {coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]};
    }
}

void Mesh::cubeMesh() {
    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1000);
    gmsh::model::occ::synchronize();
}

void Mesh::initMesh() {
    if (params == nullptr) {
        gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.1);
        gmsh::option::setNumber("Mesh.MaxNumThreads3D", 10);
        gmsh::option::setNumber("General.Verbosity", 1);

        gmsh::model::mesh::generate(3);

        gmsh::model::mesh::setOrder(1);
    }
    else {
        gmsh::option::setNumber("Mesh.CharacteristicLengthMax", params -> h);
        gmsh::option::setNumber("Mesh.MaxNumThreads3D", params -> max_num_threads);
        gmsh::option::setNumber("General.Verbosity", params -> verbosity);

        gmsh::model::mesh::generate(3);

        gmsh::model::mesh::setOrder(params -> element_order);
    }

    getNodesCoordinates();

    // get elements
    gmsh::model::mesh::getElements(elems.elementTypes, elems.elemTags, elems.nTags, 3);
    // get element type
    elems.elementType = elems.elementTypes[0];
    // element tags
    elems.elementTags = elems.elemTags[0];
    // node tags
    elems.nodeTags = elems.nTags[0];

    // get integration points
    std::string intRule = "Gauss" + std::to_string(params -> quadrature_precision);

    gmsh::model::mesh::getIntegrationPoints(elems.elementType, intRule, elems.localCoord, elems.weights);
    elems.noIntegrationPoints = int(elems.localCoord.size()) / 3;

    // get basis function values at integration points
    std::string functionSpaceType = "Lagrange" + std::to_string(params -> element_order);
    int numOrient, numComp;

    gmsh::model::mesh::getBasisFunctions(elems.elementType, elems.localCoord, functionSpaceType, numComp, elems.basisFunctionsValues, numOrient);

    elems.noNodesPerElement = utils::binom(int(params -> element_order) + 3, 3);

    // get gradients at integration points
    functionSpaceType = "GradLagrange" + std::to_string(params -> element_order);

    gmsh::model::mesh::getBasisFunctions(elems.elementType, elems.localCoord, functionSpaceType, numComp, elems.basisFunctionsGradients, numOrient);


    // get jacobians
    std::vector<double> jacobianCoords = {0.25, 0.25, 0.25};

    gmsh::model::mesh::preallocateJacobians(elems.elementType, 1, false, true, false, elems.jacobians, elems.determinants, elems.coord);
    gmsh::model::mesh::getJacobians(elems.elementType, jacobianCoords, elems.jacobians, elems.determinants, elems.coord);

}
