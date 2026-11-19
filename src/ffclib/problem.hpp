#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

// fclib's global frictional contact problem, extended by compliant bilateral rows (FORMAT.md):
//
//   M v = H r + G lambda + f
//   G^T v + b + C lambda = 0                   C = 0: fclib's perfect constraints
//   u = H^T v + w,   K_mu* ∋ u_hat ⊥ r ∈ K_mu   per contact, local frame (N, T1, T2)
//
// with the friction cone K_mu = {r : ||r_T|| <= mu r_N} and Coulomb's law u_hat = u + (mu ||u_T||, 0, 0)
// (De Saxce). Sizes: nv velocities, nc contacts, nb bilateral rows. M is symmetric positive definite, C
// symmetric positive semidefinite (a compliance).
namespace ffclib {

using Vector = Eigen::VectorXd;
// Compressed columns with 32-bit indices, the index width of fclib's files.
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

struct Problem {
    SparseMatrix M;  // nv x nv
    SparseMatrix H;  // nv x 3 nc, the columns (N, T1, T2) of each contact
    SparseMatrix G;  // nv x nb
    SparseMatrix C;  // nb x nb, both triangles stored; no entries: perfect constraints
    Vector f;        // nv
    Vector w;        // 3 nc
    Vector b;        // nb
    Vector mu;       // nc, friction coefficients (0: frictionless)

    int nv() const { return static_cast<int>(M.rows()); }
    int nc() const { return static_cast<int>(mu.size()); }
    int nb() const { return static_cast<int>(b.size()); }
};

struct Solution {
    Vector v;       // nv, velocity
    Vector r;       // 3 nc, contact percussions
    Vector lambda;  // nb, bilateral percussions
};

// Throws std::invalid_argument unless the sizes match, mu >= 0 and C is symmetric. M is checked where it is
// factorized (MassInverse); that C is positive semidefinite is assumed.
void check(const Problem& pb);

// v = M^-1 (f + H r + G lambda), the velocity that the percussions produce.
Vector velocity(const Problem& pb, const Vector& r, const Vector& lambda);

// M^-1 applied to vectors and sparse matrices: a diagonal M directly, otherwise through a sparse LDL^T
// factorization of its lower triangle.
class MassInverse {
   public:
    // Throws std::invalid_argument unless M is positive definite.
    explicit MassInverse(const SparseMatrix& M);
    Vector operator()(const Vector& x) const;
    SparseMatrix operator()(const SparseMatrix& X) const;

   private:
    bool m_diagonal{true};
    Vector m_inverse_diagonal;  // if M is diagonal
    Eigen::SimplicialLDLT<SparseMatrix> m_ldlt;
};

}  // namespace ffclib
