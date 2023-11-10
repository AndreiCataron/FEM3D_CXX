#ifndef FEM_MESH_HPP
#define FEM_MESH_HPP

#include <string>
#include <memory>
#include <set>
#include <eigen3/Eigen/Dense>
#include "params.hpp"

class Mesh {
public:
    // mesh elements
    struct MeshElements {
        // store for each node [tag : coordinates]
        std::unordered_map<std::size_t, std::tuple<double, double, double> > node_coordinates = {};
        // store for each node on the boundary the tag of the respective boundary surface
        std::unordered_map<std::size_t, std::size_t> boundaryTags = {};

        // elements
        std::vector<int> elementTypes = {};
        std::vector<std::vector<std::size_t> > elemTags = {}, nTags = {};
        int elementType = -1;
        std::vector<std::size_t> elementTags = {};
        std::vector<std::size_t> nodeTags = {};
        // faces
        std::vector<std::size_t> facesVerticesTags = {};
        // neumann boundary nodes
        std::set<std::size_t> neumannBoundaryNodes = {};

        // boundary
        std::vector<std::pair<int, int> > boundary = {};

        // integration points
        std::vector<double> localCoord = {}, weights = {};
        int noIntegrationPoints = -1;

        // basis function values at integration  points
        std::vector<double> basisFunctionsValues = {};
        std::vector<double> basisFunctionsGradients = {};
        int noNodesPerElement = -1;

        // jacobians
        std::vector<double> jacobians = {}, determinants = {}, globalCoord = {};
        // inverse jacobians
        std::vector<Eigen::MatrixXd> inverse_jacobians;
    };

public:
    std::shared_ptr<Params> params;

    MeshElements elems;

public:
    // constructors
    Mesh();
    Mesh(int, char**);
    Mesh(int, char**, const std::string&);
    Mesh(int, char**, std::shared_ptr<Params> const&);
//    // copy constructor
//    Mesh(const Mesh&);
//    // operator=
//    Mesh& operator=(const Mesh&);
    // destructor
    ~Mesh();

    friend class FEM3DVector;
    friend class LinearElasticity3D;

    void getNodesCoordinates();
    //void getBoundaryTriangles();
    static void cubeMesh();
    void computeInverseJacobians();
    void initMesh();
};


#endif
