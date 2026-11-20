#include <ffclib/solver/gauss_seidel.hpp>

#include <algorithm>
#include <chrono>
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

}  // namespace

Result gaussSeidel(const Problem& pb, const GaussSeidelSettings& st, const Solution* guess) {
    const auto start = std::chrono::steady_clock::now();
    const int nc = pb.nc(), nb = pb.nb();
    if (!(st.relaxation > 0 && st.relaxation < 2) || !(st.atol > 0) || !(st.rtol >= 0) || st.max_sweeps < 0) throw std::invalid_argument("ffclib Gauss-Seidel: invalid settings");
    if (guess && (guess->r.size() != 3 * nc || guess->lambda.size() != nb)) throw std::invalid_argument("ffclib Gauss-Seidel: the guess does not match the problem");

    const MassInverse Minv(pb.M);
    const SparseMatrix MinvH = Minv(pb.H), MinvG = Minv(pb.G);
    // Diagonal entries of the Delassus matrix of the bilateral rows, G_k^T M^-1 G_k + C_kk.
    Vector rowDiagonal(nb);
    for (int k = 0; k < nb; ++k) rowDiagonal[k] = pb.G.col(k).dot(MinvG.col(k)) + pb.C.coeff(k, k);
    // Step of each contact, rho_c = 1 / lambda_max(H_c^T M^-1 H_c); 0 if the contact acts on no velocity.
    Vector rho = Vector::Zero(nc);
    for (int c = 0; c < nc; ++c) {
        Eigen::Matrix3d W;
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) W(i, j) = W(j, i) = pb.H.col(3 * c + i).dot(MinvH.col(3 * c + j));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig;
        eig.computeDirect(W, Eigen::EigenvaluesOnly);
        if (eig.eigenvalues()[2] > 0) rho[c] = 1 / eig.eigenvalues()[2];
    }

    Result res;
    Solution& s = res.solution;
    s.r = guess ? guess->r : Vector::Zero(3 * nc);
    s.lambda = guess ? guess->lambda : Vector::Zero(nb);
    Vector v = velocity(pb, s.r, s.lambda);

    // Largest error of the sweep, relative to its tolerance (<= 1: converged).
    auto relative = [&](double error, double scale) { return error / (st.atol + st.rtol * scale); };
    int sweep = 0;
    res.status = "sweep limit";
    while (sweep < st.max_sweeps) {
        double worst = 0;
        for (int k = 0; k < nb; ++k) {
            if (!(rowDiagonal[k] > 0)) continue;
            const double velocityTerms = pb.G.col(k).dot(v) + pb.b[k];
            const double residual = velocityTerms + pb.C.col(k).dot(s.lambda);
            worst = nanMax(worst, relative(std::abs(residual), std::abs(velocityTerms)));
            const double delta = -st.relaxation * residual / rowDiagonal[k];
            s.lambda[k] += delta;
            v += delta * MinvG.col(k);
        }
        for (int c = 0; c < nc; ++c) {
            if (rho[c] == 0) continue;
            Eigen::Vector3d u;
            for (int j = 0; j < 3; ++j) u[j] = pb.H.col(3 * c + j).dot(v) + pb.w[3 * c + j];
            Eigen::Vector3d uHat = u;
            if (st.law == Law::Coulomb) uHat[0] += pb.mu[c] * std::hypot(u[1], u[2]);
            const double step = st.relaxation * rho[c];
            const Eigen::Vector3d rc = s.r.segment<3>(3 * c);
            const Eigen::Vector3d delta = projectOntoCone(rc - step * uHat, pb.mu[c]) - rc;
            worst = nanMax(worst, relative(std::abs(delta[0]) / step, std::abs(u[0])));
            worst = nanMax(worst, relative(std::hypot(delta[1], delta[2]) / step, std::hypot(u[1], u[2])));
            s.r.segment<3>(3 * c) += delta;
            for (int j = 0; j < 3; ++j) v += delta[j] * MinvH.col(3 * c + j);
        }
        ++sweep;
        if (std::isnan(worst)) {
            res.status = "diverged (NaN)";
            break;
        }
        if (worst <= 1) {
            res.converged = true;
            res.status = "converged";
            break;
        }
    }
    res.iterations = sweep;
    s.v = velocity(pb, s.r, s.lambda);  // without the round-off accumulated by the updates
    res.errors = errors(pb, s);
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return res;
}

}  // namespace ffclib
