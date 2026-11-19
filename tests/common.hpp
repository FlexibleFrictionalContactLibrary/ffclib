#pragma once

// Helpers shared by the tests: a failure counter and a random problem.

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <ffclib/io.hpp>

namespace test {

inline int failures = 0;

inline void expect(bool ok, const std::string& what) {
    if (!ok) {
        ++failures;
        std::cerr << "FAILED: " << what << '\n';
    }
}

// One point mass (dofs 0-2, m = 2) on the ground with the contact frame (N, T1, T2) = (e_z, e_x, e_y),
// mu = 0.5, and a slider (dof 3, m = 2) on a compliant spring row: G = e_3, C = 0.5, b = -0.2, so that
// lambda = 0.2 and v_3 = 0.1. f is the momentum without forces:
//   kSticking: r = (1, -0.3, 0), v = (0, 0, 0, 0.1)
//   kSliding:  r = (1, -0.5, 0), v = (0.75, 0, 0, 0.1) (Coulomb's law)
inline const Eigen::Vector4d kSticking(0.3, 0, -1, 0);
inline const Eigen::Vector4d kSliding(2, 0, -1, 0);

inline ffclib::Problem pointMass(const Eigen::Vector4d& f) {
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
    pb.f = f;
    pb.w = ffclib::Vector::Zero(3);
    pb.mu = ffclib::Vector::Constant(1, 0.5);
    pb.b = ffclib::Vector::Constant(1, -0.2);
    return pb;
}

inline int finish(const char* name) {
    std::cout << name << (failures ? ": FAILED " + std::to_string(failures) : std::string(": all passed")) << '\n';
    return failures ? 1 : 0;
}

// `bodies` rigid bodies (6 dofs each), nc contacts (every third frictionless), nb bilateral rows. With `coupledMass`,
// M couples two rotational dofs per body (not diagonal, still positive definite). With `compliant`, the first
// half of the rows is compliant: a 2 x 2 block of rows 0 and 1, the others diagonal. Also a random solution,
// two guesses and info strings.
inline ffclib::Record randomRecord(int bodies, int nc, int nb, bool compliant, bool coupledMass, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> u(-1, 1), pos(0.5, 2);
    auto random = [&](Eigen::Index size) { return ffclib::Vector(ffclib::Vector::NullaryExpr(size, [&] { return u(gen); })); };
    const int nv = 6 * bodies;
    ffclib::Record rec;
    ffclib::Problem& pb = rec.problem;
    std::vector<Eigen::Triplet<double>> t;
    for (int i = 0; i < nv; ++i) t.emplace_back(i, i, pos(gen));
    if (coupledMass)
        for (int b = 0; b < bodies; ++b) {
            const double c = 0.1 * u(gen);
            t.emplace_back(6 * b + 3, 6 * b + 4, c);
            t.emplace_back(6 * b + 4, 6 * b + 3, c);
        }
    pb.M.resize(nv, nv);
    pb.M.setFromTriplets(t.begin(), t.end());
    t.clear();
    for (int c = 0; c < nc; ++c)
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 6; ++i) t.emplace_back(6 * (c % bodies) + i, 3 * c + j, u(gen));
    pb.H.resize(nv, 3 * nc);
    pb.H.setFromTriplets(t.begin(), t.end());
    t.clear();
    for (int k = 0; k < nb; ++k)
        for (int i = 0; i < 6; ++i) t.emplace_back(6 * (k % bodies) + i, k, u(gen));
    pb.G.resize(nv, nb);
    pb.G.setFromTriplets(t.begin(), t.end());
    t.clear();
    if (compliant) {
        for (int k = 0; k < nb / 2; ++k) t.emplace_back(k, k, 1 + pos(gen));
        if (nb >= 2) {
            t.emplace_back(0, 1, 0.2);
            t.emplace_back(1, 0, 0.2);
        }
    }
    pb.C.resize(nb, nb);
    pb.C.setFromTriplets(t.begin(), t.end());
    pb.f = random(nv);
    pb.w = random(3 * nc);
    pb.b = random(nb);
    pb.mu = ffclib::Vector::NullaryExpr(nc, [&] { return 0.3 + 0.2 * u(gen); });
    for (int c = 2; c < nc; c += 3) pb.mu[c] = 0;

    rec.info = {"random problem", "test::randomRecord", compliant ? "compliant" : ""};
    rec.solution = ffclib::Solution{random(nv), random(3 * nc), random(nb)};
    rec.guesses = {ffclib::Solution{random(nv), random(3 * nc), random(nb)}, ffclib::Solution{random(nv), random(3 * nc), random(nb)}};
    return rec;
}

// randomRecord(4, 10, 6, ...) made feasible for v = 0: w_c inside K_mu* and b = 0 on the rows without
// compliance. The convexified law then has a solution.
inline ffclib::Problem feasibleProblem() {
    ffclib::Problem pb = randomRecord(4, 10, 6, true, true, 13).problem;
    for (int c = 0; c < pb.nc(); ++c) pb.w[3 * c] = std::abs(pb.w[3 * c]) + pb.mu[c] * std::hypot(pb.w[3 * c + 1], pb.w[3 * c + 2]);
    for (int k = 0; k < pb.nb(); ++k)
        if (pb.C.coeff(k, k) == 0) pb.b[k] = 0;
    return pb;
}

}  // namespace test
