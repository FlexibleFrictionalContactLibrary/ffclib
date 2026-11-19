#pragma once

// The types of CardilloCxx's misc/types.hpp that the copied solver wrappers use (README.md in this folder).
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace cardillo {

using real_t = double;
using VectorXr = Eigen::Matrix<real_t, Eigen::Dynamic, 1>;
using Vector3r = Eigen::Matrix<real_t, 3, 1>;
using CscMatrix = Eigen::SparseMatrix<real_t, Eigen::ColMajor>;

}  // namespace cardillo
