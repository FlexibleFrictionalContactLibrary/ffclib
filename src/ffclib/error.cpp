#include <ffclib/error.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <Eigen/Eigenvalues>

#include <ffclib/cone.hpp>

namespace ffclib {

namespace {

// Maximum that keeps NaN (std::max may drop it).
double nanMax(double a, double b) {
    return std::isnan(a) || std::isnan(b) ? std::numeric_limits<double>::quiet_NaN() : std::max(a, b);
}

double maxAbs(const Vector& x) {
    double m = 0;
    for (const double xi : x) m = nanMax(m, std::abs(xi));
    return m;
}

// Dot product of column i of A and column j of B (both with ascending row indices).
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

}  // namespace

Errors errors(const Problem& pb, const Vector& r, const Vector& lambda) {
    if (r.size() != 3 * pb.nc() || lambda.size() != pb.nb()) throw std::invalid_argument("ffclib errors: r or lambda does not match the problem");
    const MassInverse Minv(pb.M);

    // The velocity that the percussions produce.
    Vector force = pb.f;
    force.noalias() += pb.H * r;
    force.noalias() += pb.G * lambda;
    const Vector v = Minv(force);

    Vector rows = pb.b;
    rows.noalias() += pb.G.transpose() * v;
    rows.noalias() += pb.C * lambda;
    const double bilateral = maxAbs(rows);

    Vector u = pb.w;
    u.noalias() += pb.H.transpose() * v;
    const SparseMatrix MinvH = Minv(pb.H);
    double coulomb = 0, convexified = 0;
    for (int c = 0; c < pb.nc(); ++c) {
        Eigen::Matrix3d D;  // H_c^T M^-1 H_c
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) D(i, j) = D(j, i) = columnDot(pb.H, 3 * c + i, MinvH, 3 * c + j);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig;
        eig.computeDirect(D, Eigen::EigenvaluesOnly);
        const double lambdaMax = eig.eigenvalues()[2];
        if (!(lambdaMax > 0)) continue;
        const double rho = 1 / lambdaMax;

        const double mu = pb.mu[c];
        const Eigen::Vector3d rc = r.segment<3>(3 * c), uc = u.segment<3>(3 * c);
        const Eigen::Vector3d shifted(uc[0] + mu * std::hypot(uc[1], uc[2]), uc[1], uc[2]);
        coulomb = nanMax(coulomb, (rc - projectOntoCone(rc - rho * shifted, mu)).norm() * lambdaMax);
        convexified = nanMax(convexified, (rc - projectOntoCone(rc - rho * uc, mu)).norm() * lambdaMax);
    }
    return {nanMax(bilateral, coulomb), nanMax(bilateral, convexified)};
}

}  // namespace ffclib
