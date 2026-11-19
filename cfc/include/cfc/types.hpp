#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace cfc {

using Vector = Eigen::VectorXd;
using Vector3 = Eigen::Vector3d;
using Matrix3 = Eigen::Matrix3d;
// Compressed columns with 32-bit indices, the index width of fclib's files.
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

}  // namespace cfc
