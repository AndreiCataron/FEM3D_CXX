#ifndef FEM_FEM3D_HPP
#define FEM_FEM3D_HPP

#endif //FEM_FEM3D_HPP

class FEM3D{
public:
    struct Params{
        double h;
    };

private:
    const Params params_;

public:
    // constructor
    FEM3D(const Params&);

    // methods
    void setBoundaryConditions();
};
