// The error measures:
//   1. closed-form solutions: a sticking and a sliding point mass, a slider on a compliant spring row, each
//      condition perturbed alone;
//   2. NaN in a candidate gives NaN in the errors that use it;
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

// The point mass and slider of test::pointMass().
void testClosedForm() {
    ffclib::Problem pb = test::pointMass(test::kSticking);
    const ffclib::Vector lambda = ffclib::Vector::Constant(1, 0.2);
    const double tol = 1e-15;

    ffclib::Errors e = ffclib::errors(pb, {Eigen::Vector4d(0, 0, 0, 0.1), Eigen::Vector3d(1, -0.3, 0), lambda});
    expect(e.balance < tol && e.bilateral < tol && e.coulomb < tol && e.convexified < tol, "sticking contact: every condition holds");

    // Sliding: u = (0, 0.75, 0). Coulomb's law holds; the convexified law needs u_N = mu ||u_T|| = 0.375. With
    // rho = 2 (H_c^T M^-1 H_c = I / 2), r - rho u = (1, -2, 0) projects onto (1.6, -0.8, 0), so its error is
    // ||(-0.6, 0.3, 0)|| / 2.
    pb.f = test::kSliding;
    const ffclib::Solution slide{Eigen::Vector4d(0.75, 0, 0, 0.1), Eigen::Vector3d(1, -0.5, 0), lambda};
    e = ffclib::errors(pb, slide);
    expect(e.balance < tol && e.bilateral < tol && e.coulomb < tol && std::abs(e.convexified - 0.5 * std::sqrt(0.45)) < tol,
           "sliding contact: Coulomb's law holds, the convexified law not");

    // lambda + 1e-3 at the same velocity: the slider's balance is off by G 1e-3, its spring row by C 1e-3.
    ffclib::Solution perturbed = slide;
    perturbed.lambda[0] += 1e-3;
    e = ffclib::errors(pb, perturbed);
    expect(std::abs(e.balance - 1e-3) < tol && std::abs(e.bilateral - 0.5e-3) < tol && e.coulomb < tol, "perturbed bilateral percussion");

    // v_3 + 1e-3 at the same percussions: the balance is off by M 1e-3, the spring row by G^T 1e-3.
    perturbed = slide;
    perturbed.v[3] += 1e-3;
    e = ffclib::errors(pb, perturbed);
    expect(std::abs(e.balance - 2e-3) < tol && std::abs(e.bilateral - 1e-3) < tol && e.coulomb < tol, "perturbed velocity of the slider");

    perturbed = slide;
    perturbed.r[1] = std::numeric_limits<double>::quiet_NaN();
    e = ffclib::errors(pb, perturbed);
    expect(std::isnan(e.balance) && std::isnan(e.coulomb) && std::isnan(e.convexified) && e.bilateral < tol, "NaN in r gives NaN in the errors that use r");
}

// The definition in error.hpp, evaluated with dense matrices.
ffclib::Errors denseErrors(const ffclib::Problem& pb, const ffclib::Solution& s) {
    const Eigen::MatrixXd M = pb.M.toDense(), H = pb.H.toDense(), G = pb.G.toDense(), C = pb.C.toDense();
    const Eigen::LLT<Eigen::MatrixXd> llt(M);
    const ffclib::Vector u = H.transpose() * s.v + pb.w;
    ffclib::Errors e;
    e.balance = (M * s.v - pb.f - H * s.r - G * s.lambda).cwiseAbs().maxCoeff();
    e.bilateral = pb.nb() > 0 ? (G.transpose() * s.v + pb.b + C * s.lambda).cwiseAbs().maxCoeff() : 0;
    for (int c = 0; c < pb.nc(); ++c) {
        const Eigen::MatrixXd Hc = H.middleCols(3 * c, 3);
        const Eigen::Matrix3d D = Hc.transpose() * llt.solve(Hc);
        const double lambdaMax = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(D).eigenvalues()[2];
        const Eigen::Vector3d rc = s.r.segment<3>(3 * c), uc = u.segment<3>(3 * c);
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
        const ffclib::Errors e = ffclib::errors(rec.problem, *rec.solution);
        const ffclib::Errors d = denseErrors(rec.problem, *rec.solution);
        auto near = [](double a, double b) { return std::abs(a - b) < 1e-12 * std::max(1.0, std::abs(b)); };
        expect(near(e.balance, d.balance) && near(e.bilateral, d.bilateral) && near(e.coulomb, d.coulomb) && near(e.convexified, d.convexified),
               std::string("random problem against the dense definition") + (coupledMass ? ", non-diagonal M" : ", diagonal M"));
    }
}

}  // namespace

int main() {
    testClosedForm();
    testAgainstDense();
    return test::finish("test_error");
}
