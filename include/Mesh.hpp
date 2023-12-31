#ifndef FEM_MESH_HPP
#define FEM_MESH_HPP

#include <string>
#include <memory>
#include <set>
#include <eigen3/Eigen/Dense>
#include "params.hpp"
#include "utils.hpp"
#include "MeshInfo.hpp"

using CoordTuple = std::tuple<double, double, double>;

class Mesh {
protected:
    std::shared_ptr<Params> params;

    GlobalData global;
    BoundaryData bdry;
    ElementData elems;

public:
    // constructors
    Mesh();
    Mesh(int, char**);
    Mesh(int, char**, const std::string&);
    Mesh(int, char**, std::shared_ptr<Params> const&);
    ~Mesh();

    friend class FEM3D;
    friend class FEM3DVector;
    friend class LinearElasticity3D;

    void getNodesCoordinates();
    static void cubeMesh();
    void computeInverseJacobians();
    void initMesh();
    static void showMesh(int, char**);
};


#endif
