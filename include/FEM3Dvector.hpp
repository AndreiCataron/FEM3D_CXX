#ifndef FEM_FEM3DVECTOR_HPP
#define FEM_FEM3DVECTOR_HPP

#endif //FEM_FEM3DVECTOR_HPP

#include "FEM3D.hpp"
#include <string>

class FEM3DVector : public FEM3D {
public:
    struct ParamsVector : FEM3D::Params{
        std::string g;
    };
private:
    const ParamsVector params3d_;
public:
    FEM3DVector(const ParamsVector&);
    void setBoundaryConditions();
};