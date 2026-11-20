// The conic solvers:
//   1. the cone program: at a candidate (v, r, lambda) its stationarity and slack rows are the problem's
//      equations, and solutionFromConeProgram() recovers the candidate;
//   2. every built solver (Clarabel, SCS, QOCO, ConicXX, cuOpt) on the closed-form point mass (convexified law) and on a random
//      feasible problem. A solver without a solution because of its platform (cuOpt without GPU) is skipped.

#include <string>
#include <utility>
#include <vector>

#include <ffclib/solver/conic/cone_program.hpp>
#include <ffclib/solver/solver.hpp>
#ifdef FFCLIB_WITH_CLARABEL
#include <ffclib/solver/conic/clarabel.hpp>
#endif
#ifdef FFCLIB_WITH_SCS
#include <ffclib/solver/conic/scs.hpp>
#endif
#ifdef FFCLIB_WITH_QOCO
#include <ffclib/solver/conic/qoco.hpp>
#endif
#ifdef FFCLIB_WITH_CONICXX
#include <ffclib/solver/conic/conicxx.hpp>
#endif
#ifdef FFCLIB_WITH_CUOPT
#include <ffclib/solver/conic/cuopt.hpp>
#endif

#include "common.hpp"

namespace {

using test::expect;

double maxDiff(const ffclib::Vector& a, const ffclib::Vector& b) {
    return a.size() == b.size() ? (a - b).lpNorm<Eigen::Infinity>() : 1;
}

void testConeProgram() {
    const ffclib::Record rec = test::randomRecord(3, 7, 6, true, true, 3);
    const ffclib::Problem& pb = rec.problem;
    const int nv = pb.nv(), nb = pb.nb(), npi = nb / 2;  // the first half of the rows is compliant
    ffclib::Vector v = rec.solution->v, r = rec.solution->r, lambda = rec.solution->lambda;
    std::vector<int> frictionless, frictional;
    for (int c = 0; c < pb.nc(); ++c) {
        (pb.mu[c] > 0 ? frictional : frictionless).push_back(c);
        if (pb.mu[c] == 0) r.segment<2>(3 * c + 1).setZero();
    }

    const ffclib::ConeProgram cp = ffclib::coneProgram(pb);
    const int rows = nb + static_cast<int>(frictionless.size() + 3 * frictional.size());
    expect(cp.P.cols() == nv + npi && cp.A.rows() == rows && cp.dims.zero == nb && cp.dims.nonneg == static_cast<int>(frictionless.size()) && cp.dims.soc == std::vector<int>(frictional.size(), 3),
           "cone program: sizes");

    // x = (v, pi = lambda of the compliant rows), z = (-lambda, S_c^-1 r_c).
    ffclib::Vector x(nv + npi), z(rows);
    x << v, lambda.head(npi);
    z.head(nb) = -lambda;
    int row = nb;
    for (int c : frictionless) z[row++] = r[3 * c];
    for (int c : frictional) {
        z.segment<3>(row) << r[3 * c], r[3 * c + 1] / pb.mu[c], r[3 * c + 2] / pb.mu[c];
        row += 3;
    }
    const ffclib::SparseMatrix P = cp.P.selfadjointView<Eigen::Upper>();
    const ffclib::Vector stationarity = P * x + cp.q + cp.A.transpose() * z;
    const ffclib::Vector balance = pb.M * v - pb.f - pb.H * r - pb.G * lambda;
    const double tol = 1e-12;
    expect(maxDiff(stationarity.head(nv), balance) < tol && stationarity.tail(npi).lpNorm<Eigen::Infinity>() < tol, "cone program: stationarity is the balance");

    const ffclib::Vector slack = cp.b - cp.A * x;
    const ffclib::Vector u = pb.H.transpose() * v + pb.w;
    ffclib::Vector expected(rows);
    expected.head(nb) = -(pb.G.transpose() * v + pb.b + pb.C * lambda);
    row = nb;
    for (int c : frictionless) expected[row++] = u[3 * c];
    for (int c : frictional) {
        expected.segment<3>(row) << u[3 * c], pb.mu[c] * u[3 * c + 1], pb.mu[c] * u[3 * c + 2];
        row += 3;
    }
    expect(maxDiff(slack, expected) < tol, "cone program: slacks are the bilateral rows and S_c u_c");

    const ffclib::Solution s = ffclib::solutionFromConeProgram(pb, x, z);
    expect(maxDiff(s.v, v) < tol && maxDiff(s.r, r) < tol && maxDiff(s.lambda, lambda) < tol, "cone program: the solution from x and z");
}

using SolverFunction = ffclib::Result (*)(const ffclib::Problem&, const ffclib::ConicSettings&);

void testSolver(const std::string& name, SolverFunction solve) {
    ffclib::ConicSettings settings;
    settings.atol = settings.rtol = 1e-9;
    // With data of order one the balance and the bilateral rows (the stationarity and the zero rows of the cone
    // program) meet the native tolerance; the contact error falls only with the square root of the gap: on the
    // random problem the interior point methods reach 3e-6 to 4e-6 at native 1e-9 (SCS 1e-14), on the point mass
    // the percussions agree to 1e-6.
    const double feasibility = settings.atol, pointMassAccuracy = 1e-6, randomAccuracy = 1e-5;

    // Convexified law of the sliding point mass: r = (1.6, -0.8, 0), lambda = 0.2, v = (0.6, 0, 0.3, 0.1).
    ffclib::Result res = solve(test::pointMass(test::kSliding), settings);
    if (res.status.rfind("no solution", 0) == 0) {
        std::cout << name << ": skipped (" << res.status << ")\n";
        return;
    }
    expect(res.converged && maxDiff(res.solution.r, Eigen::Vector3d(1.6, -0.8, 0)) < pointMassAccuracy && maxDiff(res.solution.lambda, ffclib::Vector::Constant(1, 0.2)) < pointMassAccuracy &&
               maxDiff(res.solution.v, Eigen::Vector4d(0.6, 0, 0.3, 0.1)) < pointMassAccuracy,
           name + ": point mass, " + res.status);

    res = solve(test::feasibleProblem(), settings);
    expect(res.converged && res.errors.balance < feasibility && res.errors.bilateral < feasibility && res.errors.convexified < randomAccuracy, name + ": random problem, " + res.status);
}

}  // namespace

int main() {
    testConeProgram();
    const std::vector<std::pair<std::string, SolverFunction>> solvers = {
#ifdef FFCLIB_WITH_CLARABEL
        {"clarabel", &ffclib::clarabel},
#endif
#ifdef FFCLIB_WITH_SCS
        {"scs", &ffclib::scs},
#endif
#ifdef FFCLIB_WITH_QOCO
        {"qoco", &ffclib::qoco},
#endif
#ifdef FFCLIB_WITH_CONICXX
        {"conicxx", &ffclib::conicxx},
#endif
#ifdef FFCLIB_WITH_CUOPT
        {"cuopt", &ffclib::cuopt},
#endif
    };
    for (const auto& [name, solve] : solvers) testSolver(name, solve);
    return test::finish("test_conic");
}
