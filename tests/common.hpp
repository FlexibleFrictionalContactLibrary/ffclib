#pragma once

// Helpers shared by the tests: a failure counter and a random problem.

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

}  // namespace test
