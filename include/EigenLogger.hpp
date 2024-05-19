#ifndef FEM_EIGENLOGGER_HPP
#define FEM_EIGENLOGGER_HPP

#include "eigen_common.hpp"
#include <fstream>
#include <utility>

class EigenLogger {
private:
    std::string file;
    std::ofstream myFile;
public:
    explicit EigenLogger(std::string f) : file(std::move(f)) {}

    void outputSparseMatrixtoFile(const Eigen::SparseMatrix<double>&& m) {
        myFile.open(file, std::ios::out | std::ios::trunc);

        for (int k = 0; k < m.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
                myFile << it.row() << ' ' << it.col() << ' ' << it.value() << ' ';
            }
        }

        myFile.close();
    };

    void outputVectortoFile(const Eigen::VectorXd&& v) {
        myFile.open(file, std::ios::out | std::ios::trunc);

        myFile << v.cols() << ' ' << v.rows();
        for (auto it = v.data(); it != v.data() + v.size(); it++) {
            myFile << *it << ' ';
        }

        myFile.close();
    }

    void setFile(std::string f) {
        file = std::move(f);
    }

};

#endif //FEM_EIGENLOGGER_HPP
