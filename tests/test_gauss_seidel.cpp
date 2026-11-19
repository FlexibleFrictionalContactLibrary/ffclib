// The Gauss-Seidel solver:
//   1. the closed-form point mass and slider (test::pointMass) for both laws;
//   2. a random feasible problem: the convexified law converges for omega = 1 and 1.5;
//   3. a converged guess needs one sweep, the sweep limit stops, invalid settings throw.

#include <cmath>
#include <stdexcept>
#include <string>

#include <ffclib/solver/gauss_seidel.hpp>

#include "common.hpp"

namespace {

using test::expect;

bool near(const ffclib::Vector& a, const ffclib::Vector& b, double tol) {
    return a.size() == b.size() && (a - b).lpNorm<Eigen::Infinity>() <= tol;
}

void testClosedForm() {
    ffclib::GaussSeidelSettings st;
    st.atol = 1e-12;
    const double tol = 1e-10;  // on r, lambda and v: the errors bound them through M and C of order one

    ffclib::Result res = ffclib::gaussSeidel(test::pointMass(test::kSticking), st);
    expect(res.converged && near(res.solution.r, Eigen::Vector3d(1, -0.3, 0), tol) && near(res.solution.lambda, ffclib::Vector::Constant(1, 0.2), tol) &&
               near(res.solution.v, Eigen::Vector4d(0, 0, 0, 0.1), tol),
           "sticking contact, Coulomb's law");

    res = ffclib::gaussSeidel(test::pointMass(test::kSliding), st);
    expect(res.converged && near(res.solution.r, Eigen::Vector3d(1, -0.5, 0), tol) && near(res.solution.v, Eigen::Vector4d(0.75, 0, 0, 0.1), tol), "sliding contact, Coulomb's law");

    // The convexified law: u_N = mu ||u_T|| and r on the cone's boundary give r = (1.6, -0.8, 0), v = (0.6, 0, 0.3, 0.1).
    st.law = ffclib::Law::Convexified;
    res = ffclib::gaussSeidel(test::pointMass(test::kSliding), st);
    expect(res.converged && near(res.solution.r, Eigen::Vector3d(1.6, -0.8, 0), tol) && near(res.solution.v, Eigen::Vector4d(0.6, 0, 0.3, 0.1), tol) && res.errors.coulomb > 0.1,
           "sliding contact, convexified law");
}

void testRandom() {
    const ffclib::Problem pb = test::feasibleProblem();
    for (const double omega : {1.0, 1.5}) {
        ffclib::GaussSeidelSettings st;
        st.law = ffclib::Law::Convexified;
        st.atol = 1e-10;
        st.relaxation = omega;
        st.max_sweeps = 100000;
        const ffclib::Result res = ffclib::gaussSeidel(pb, st);
        // The stopping test and errors() measure the same natural map with other steps and norms.
        expect(res.converged && res.errors.convexified <= 10 * st.atol, "random problem, convexified law, omega = " + std::to_string(omega));

        // From the converged solution one sweep confirms it.
        const ffclib::Result again = ffclib::gaussSeidel(pb, st, &res.solution);
        expect(again.converged && again.iterations == 1, "a converged guess needs one sweep");

        // A relative tolerance alone converges as well.
        st.atol = 1e-14;
        st.rtol = 1e-8;
        expect(ffclib::gaussSeidel(pb, st).converged, "random problem, relative tolerance, omega = " + std::to_string(omega));
    }
}

void testLimits() {
    ffclib::GaussSeidelSettings st;
    st.atol = 1e-14;
    st.max_sweeps = 25;
    const ffclib::Result res = ffclib::gaussSeidel(test::feasibleProblem(), st);
    expect(!res.converged && res.iterations == st.max_sweeps, "the sweep limit stops");

    st.relaxation = 2;
    bool threw = false;
    try {
        ffclib::gaussSeidel(test::feasibleProblem(), st);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    expect(threw, "relaxation 2 is rejected");

    st.relaxation = 1;
    st.atol = 0;
    threw = false;
    try {
        ffclib::gaussSeidel(test::feasibleProblem(), st);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    expect(threw, "atol 0 is rejected");
}

}  // namespace

int main() {
    testClosedForm();
    testRandom();
    testLimits();
    return test::finish("test_gauss_seidel");
}
