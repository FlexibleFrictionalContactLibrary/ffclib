// The error measures:
//   1. closed-form solutions: a sticking and a sliding point mass, a slider on a compliant spring row;
//   2. NaN in a candidate gives NaN;
//   3. a random problem against a dense evaluation of the definition, for diagonal and non-diagonal M.

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <ffclib/cone.hpp>
#include <ffclib/error.hpp>

#include "common.hpp"

namespace {

using test::expect;

// One point mass (dofs 0-2, m = 2) on the ground with the contact frame (N, T1, T2) = (e_z, e_x, e_y),
// mu = 0.5, and a slider (dof 3, m = 2) on a compliant spring row: G = e_3, C = 0.5, b = -0.2, so that
// lambda = 0.2 and v_3 = 0.1. f is the momentum without forces. Sticking: f = (0.3, 0, -1, 0), r = (1, -0.3, 0),
// v = 0. Sliding: f = (2, 0, -1, 0), r = (1, -0.5, 0), v = (0.75, 0, 0).
void testClosedForm() {
    ffclib::Problem pb;
    pb.M.resize(4, 4);
    for (int i = 0; i < 4; ++i) pb.M.insert(i, i) = 2;
    pb.H.resize(4, 3);
    pb.H.insert(2, 0) = 1;
    pb.H.insert(0, 1) = 1;
    pb.H.insert(1, 2) = 1;
    pb.G.resize(4, 1);
    pb.G.insert(3, 0) = 1;
    pb.C.resize(1, 1);
    pb.C.insert(0, 0) = 0.5;
    pb.w = ffclib::Vector::Zero(3);
    pb.mu = ffclib::Vector::Constant(1, 0.5);
    pb.b = ffclib::Vector::Constant(1, -0.2);
    const ffclib::Vector lambda = ffclib::Vector::Constant(1, 0.2);
    const double tol = 1e-15;

    pb.f = Eigen::Vector4d(0.3, 0, -1, 0);
    ffclib::Errors e = ffclib::errors(pb, Eigen::Vector3d(1, -0.3, 0), lambda);
    expect(e.coulomb < tol && e.convexified < tol, "sticking contact: both laws hold");

    // Sliding: u = (0, 0.75, 0). Coulomb's law holds; the convexified law needs u_N = mu ||u_T|| = 0.375. With
    // rho = 2 (H_c^T M^-1 H_c = I / 2), r - rho u = (1, -2, 0) projects onto (1.6, -0.8, 0), so its error is
    // ||(-0.6, 0.3, 0)|| / 2.
    pb.f = Eigen::Vector4d(2, 0, -1, 0);
    const Eigen::Vector3d slide(1, -0.5, 0);
    e = ffclib::errors(pb, slide, lambda);
    expect(e.coulomb < tol && std::abs(e.convexified - 0.5 * std::sqrt(0.45)) < tol, "sliding contact: Coulomb's law holds, the convexified law not");

    // lambda + 1e-3 moves the slider by 0.5e-3 and adds 0.5e-3 through C: the spring row is off by 1e-3.
    e = ffclib::errors(pb, slide, (lambda.array() + 1e-3).matrix());
    expect(std::abs(e.coulomb - 1e-3) < tol, "perturbed bilateral percussion");

    Eigen::Vector3d nan = slide;
    nan[1] = std::numeric_limits<double>::quiet_NaN();
    e = ffclib::errors(pb, nan, lambda);
    expect(std::isnan(e.coulomb) && std::isnan(e.convexified), "NaN gives NaN");
}

// The definition in error.hpp, evaluated with dense matrices.
ffclib::Errors denseErrors(const ffclib::Problem& pb, const ffclib::Vector& r, const ffclib::Vector& lambda) {
    const Eigen::MatrixXd M = pb.M.toDense(), H = pb.H.toDense(), G = pb.G.toDense(), C = pb.C.toDense();
    const Eigen::LLT<Eigen::MatrixXd> llt(M);
    const ffclib::Vector v = llt.solve(pb.f + H * r + G * lambda);
    const ffclib::Vector u = H.transpose() * v + pb.w;
    double bilateral = pb.nb() > 0 ? (G.transpose() * v + pb.b + C * lambda).cwiseAbs().maxCoeff() : 0;
    ffclib::Errors e{bilateral, bilateral};
    for (int c = 0; c < pb.nc(); ++c) {
        const Eigen::MatrixXd Hc = H.middleCols(3 * c, 3);
        const Eigen::Matrix3d D = Hc.transpose() * llt.solve(Hc);
        const double lambdaMax = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(D).eigenvalues()[2];
        const Eigen::Vector3d rc = r.segment<3>(3 * c), uc = u.segment<3>(3 * c);
        const double mu = pb.mu[c];
        const Eigen::Vector3d coulomb(uc[0] + mu * uc.tail<2>().norm(), uc[1], uc[2]);
        e.coulomb = std::max(e.coulomb, (rc - ffclib::projectOntoCone(rc - coulomb / lambdaMax, mu)).norm() * lambdaMax);
        e.convexified = std::max(e.convexified, (rc - ffclib::projectOntoCone(rc - uc / lambdaMax, mu)).norm() * lambdaMax);
    }
    return e;
}

void testAgainstDense() {
    for (const bool coupledMass : {false, true}) {
        const ffclib::Record rec = test::randomRecord(3, 7, 6, true, coupledMass, 3);
        const ffclib::Errors e = ffclib::errors(rec.problem, rec.solution->r, rec.solution->lambda);
        const ffclib::Errors d = denseErrors(rec.problem, rec.solution->r, rec.solution->lambda);
        const double tol = 1e-12 * std::max(1.0, d.coulomb);
        expect(std::abs(e.coulomb - d.coulomb) < tol && std::abs(e.convexified - d.convexified) < tol,
               std::string("random problem against the dense definition") + (coupledMass ? ", non-diagonal M" : ", diagonal M"));
    }
}

}  // namespace

int main() {
    testClosedForm();
    testAgainstDense();
    return test::finish("test_error");
}
