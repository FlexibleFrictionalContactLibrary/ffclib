#include <cfc/residuals.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>

#include <cfc/cone.hpp>

namespace cfc {

namespace {

// NaN-propagating maximum.
double nanMax(double a, double b) {
    return (std::isnan(a) || std::isnan(b)) ? std::numeric_limits<double>::quiet_NaN() : std::max(a, b);
}

double infNorm(const Vector& x) {
    double m = 0;
    for (Eigen::Index i = 0; i < x.size(); ++i) m = nanMax(m, std::abs(x[i]));
    return m;
}

// The norm of the natural map, max(|x_N|, ||x_T||).
double mapNorm(const Vector3& x) {
    return std::max(std::abs(x[0]), std::hypot(x[1], x[2]));
}

// Dot product of column i of A and column j of B (sorted inner indices).
double columnDot(const SparseMatrix& A, Eigen::Index i, const SparseMatrix& B, Eigen::Index j) {
    double sum = 0;
    SparseMatrix::InnerIterator a(A, i), b(B, j);
    while (a && b) {
        if (a.row() < b.row()) {
            ++a;
        } else if (b.row() < a.row()) {
            ++b;
        } else {
            sum += a.value() * b.value();
            ++a;
            ++b;
        }
    }
    return sum;
}

bool isDiagonal(const SparseMatrix& M) {
    for (Eigen::Index j = 0; j < M.outerSize(); ++j)
        for (SparseMatrix::InnerIterator it(M, j); it; ++it)
            if (it.row() != j && it.value() != 0) return false;
    return true;
}

}  // namespace

struct MassInverse::Impl {
    Vector inverse_diagonal;  // used if M is diagonal
    Eigen::SimplicialLDLT<SparseMatrix> ldlt;
    bool diagonal{true};
};

MassInverse::MassInverse(const SparseMatrix& M) : m_impl(std::make_unique<Impl>()) {
    m_impl->diagonal = isDiagonal(M);
    if (m_impl->diagonal) {
        m_impl->inverse_diagonal = Vector(M.diagonal()).cwiseInverse();
        return;
    }
    m_impl->ldlt.compute(M);
    if (m_impl->ldlt.info() != Eigen::Success) throw std::runtime_error("cfc: the mass matrix is not positive definite");
}

MassInverse::~MassInverse() = default;

Vector MassInverse::apply(const Vector& x) const {
    if (m_impl->diagonal) return m_impl->inverse_diagonal.cwiseProduct(x);
    return m_impl->ldlt.solve(x);
}

SparseMatrix MassInverse::apply(const SparseMatrix& X) const {
    if (m_impl->diagonal) return m_impl->inverse_diagonal.asDiagonal() * X;
    SparseMatrix Y = m_impl->ldlt.solve(X);
    Y.prune(double(0));
    return Y;
}

double Residuals::convexifiedMax() const {
    return nanMax(nanMax(balance, bilateral), natural_map);
}
double Residuals::coulombMax() const {
    return nanMax(nanMax(balance, bilateral), natural_map_coulomb);
}

Residuals residuals(const Problem& pb, const Solution& sol, const Vector& shift) {
    const int n = pb.n(), nc = pb.nc(), p = pb.p();
    if (sol.v.size() != n || sol.r.size() != 3 * nc || sol.lambda.size() != p || shift.size() != nc) throw std::invalid_argument("cfc residuals: the solution does not match the problem");
    const MassInverse Minv(pb.M);
    Residuals res;

    // Velocity scale: the largest contact velocity of the motion without forces.
    Vector ufree = pb.w;
    ufree.noalias() += pb.H.transpose() * Minv.apply(pb.f);
    double scale = 0;
    for (int c = 0; c < nc; ++c) scale = nanMax(scale, ufree.segment<3>(3 * c).norm());
    res.velocity_scale = scale > 0 ? scale : double(1);

    // Momentum balance.
    Vector force = pb.f;
    force.noalias() += pb.H * sol.r;
    if (p > 0) force.noalias() += pb.G * sol.lambda;
    res.balance = infNorm(sol.v - Minv.apply(force));

    // Bilateral rows.
    if (p > 0) {
        Vector row = pb.b;
        row.noalias() += pb.G.transpose() * sol.v;
        if (pb.hasCompliance()) row.noalias() += pb.C * sol.lambda;
        res.bilateral = infNorm(row);
    }

    // Contacts. u = H^T v + w; the Delassus block H_c^T M^-1 H_c of contact c.
    Vector u = pb.w;
    u.noalias() += pb.H.transpose() * sol.v;
    const SparseMatrix MinvH = Minv.apply(pb.H);
    for (int c = 0; c < nc; ++c) {
        Matrix3 D;
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) D(i, j) = D(j, i) = columnDot(pb.H, 3 * c + i, MinvH, 3 * c + j);
        Eigen::SelfAdjointEigenSolver<Matrix3> eig;
        eig.computeDirect(D, Eigen::EigenvaluesOnly);
        const double lambdaMax = eig.eigenvalues()[2];
        if (!(lambdaMax > 0)) continue;
        const double rho = 1 / lambdaMax;

        const double mu = pb.mu[c];
        const Vector3 r = sol.r.segment<3>(3 * c);
        const Vector3 uc = u.segment<3>(3 * c);
        Vector3 uhat = uc;
        uhat[0] += shift[c];
        Vector3 uhatCoulomb = uc;
        uhatCoulomb[0] += deSaxce(uc, mu);

        res.natural_map = nanMax(res.natural_map, mapNorm(r - projectCoulombCone(r - rho * uhat, mu)) * lambdaMax);
        res.natural_map_coulomb = nanMax(res.natural_map_coulomb, mapNorm(r - projectCoulombCone(r - rho * uhatCoulomb, mu)) * lambdaMax);
        res.cone_primal = nanMax(res.cone_primal, (r - projectCoulombCone(r, mu)).norm() * lambdaMax);
        // dist(x, K*) = ||proj_K(-x)|| (Moreau decomposition, the polar of K* is -K).
        res.cone_dual = nanMax(res.cone_dual, projectCoulombCone(-uhat, mu).norm());
        res.complementarity = nanMax(res.complementarity, std::sqrt(std::abs(r.dot(uhat)) * lambdaMax));
    }
    return res;
}

}  // namespace cfc
