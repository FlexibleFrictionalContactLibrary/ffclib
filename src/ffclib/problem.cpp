#include <ffclib/problem.hpp>

#include <stdexcept>
#include <string>

namespace ffclib {

namespace {

void require(bool ok, const char* what) {
    if (!ok) throw std::invalid_argument(std::string("ffclib problem: ") + what);
}

}  // namespace

void check(const Problem& pb) {
    const Eigen::Index nv = pb.M.rows(), nc = pb.mu.size(), nb = pb.b.size();
    require(pb.M.cols() == nv, "M must be square");
    require(pb.H.rows() == nv && pb.H.cols() == 3 * nc, "H must be nv x 3 nc");
    require(pb.G.rows() == nv && pb.G.cols() == nb, "G must be nv x nb");
    require(pb.C.rows() == nb && pb.C.cols() == nb, "C must be nb x nb");
    require(pb.f.size() == nv, "f must have nv entries");
    require(pb.w.size() == 3 * nc, "w must have 3 nc entries");
    require((pb.mu.array() >= 0).all(), "mu must be nonnegative");
    const SparseMatrix Ct = pb.C.transpose();
    require(SparseMatrix(pb.C - Ct).norm() == 0, "C must be symmetric (both triangles stored)");
}

MassInverse::MassInverse(const SparseMatrix& M) {
    if (M.rows() != M.cols()) throw std::invalid_argument("ffclib: M must be square");
    for (Eigen::Index j = 0; j < M.outerSize() && m_diagonal; ++j)
        for (SparseMatrix::InnerIterator it(M, j); it; ++it)
            if (it.row() != j && it.value() != 0) m_diagonal = false;
    if (m_diagonal) {
        const Vector d = M.diagonal();
        if (!(d.array() > 0).all()) throw std::invalid_argument("ffclib: M must be positive definite");
        m_inverse_diagonal = d.cwiseInverse();
        return;
    }
    m_ldlt.compute(M);
    if (m_ldlt.info() != Eigen::Success || !(m_ldlt.vectorD().array() > 0).all()) throw std::invalid_argument("ffclib: M must be positive definite");
}

Vector MassInverse::operator()(const Vector& x) const {
    if (m_diagonal) return m_inverse_diagonal.cwiseProduct(x);
    return m_ldlt.solve(x);
}

SparseMatrix MassInverse::operator()(const SparseMatrix& X) const {
    if (m_diagonal) return m_inverse_diagonal.asDiagonal() * X;
    return m_ldlt.solve(X);
}

}  // namespace ffclib
