#ifndef FEM_ALGORITM_HPP
#define FEM_ALGORITM_HPP

#include "LinearElasticity3D.hpp"

class Algoritm {
private:
    std::shared_ptr<LinearElasticity3D> LE;
public:
    explicit Algoritm(std::shared_ptr<LinearElasticity3D> const&);

    void iteration();
};


#endif
