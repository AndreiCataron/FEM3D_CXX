#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Sparse>
#include <iostream>
#include <cassert>

int main() {
    std::vector<Eigen::Triplet<double> > tripletList;
    //tripletList.reserve();
    Eigen::SparseMatrix<double> mat(3 * 1, 3 * 1);
    tripletList.emplace_back(5, 6, 7);
    assert(tripletList.size() == 1);

    tripletList.emplace_back(Eigen::Triplet<double>(8, 9, 10));
    assert(tripletList.size() == 2);

}