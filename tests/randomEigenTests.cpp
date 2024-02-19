#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Dense>
#include </opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/Eigen/Sparse>
#include <iostream>
#include <cassert>

//int main() {
//    std::vector<Eigen::Triplet<double> > tripletList;
//    //tripletList.reserve();
//    Eigen::SparseMatrix<double> mat(3 * 1, 3 * 1);
//    tripletList.emplace_back(5, 6, 7);
//    assert(tripletList.size() == 1);
//
//    tripletList.emplace_back(Eigen::Triplet<double>(8, 9, 10));
//    assert(tripletList.size() == 2);
//
//    Eigen::MatrixXd m(3, 3);
//    m << 1, 0, 0, 0, 1, 0, 0, 0, 1;
//    std::vector<double> vec = {3, 2, 1};
//    //Eigen::Matrix<double, 2, 1> v(vec.data());
//    double *ptr = &vec[0];
//    Eigen::Map<Eigen::VectorXd> v(ptr, 3);
//    Eigen::VectorXd v2(2);
//    v2 << 5, 5;
//
//    std::vector<int> indexes = {0, 2};
//    std::cout << v(indexes) << '\n';
//    v(indexes) += v2;
//    std::cout << v << '\n';
//
//    std::cout << m * v;
//
//}

int main() {
    Eigen::Matrix3d hatz = Eigen::Matrix3d::Identity();

    std::vector<Eigen::Matrix3d> v, vv;
    v.emplace_back(hatz);
    vv = v;
    std::cout << vv[0];

    return 0;
}