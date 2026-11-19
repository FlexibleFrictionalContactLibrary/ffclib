#pragma once

#include <vector>

#include <ffclib/problem.hpp>

// The conic solvers solve the convexified law (no De Saxce shift). Clarabel and SCS take it as the cone program
//
//   min 1/2 x'Px + q'x   s.t.   Ax + s = b,  s in K = {0}^zero x R+^nonneg x SOC(3)^k
//
// (KKT: Px + q + A'z = 0, z in K*, s'z = 0), here in the velocity form with the cones on the rows:
//
//   x = (v, pi),  P = blkdiag(M, C_KK),  q = (-f, 0)    pi: one entry per bilateral row in K, the rows with a
//                                                      nonzero column of C
//   zero rows     G^T v + C_:K pi + s = -b              the bilateral rows
//   nonneg rows   -H_N,c^T v + s_c = w_N,c              contacts without friction (mu_c = 0)
//   SOC rows      -S_c H_c^T v + s_c = S_c w_c          contacts with friction, S_c = diag(1, mu_c, mu_c)
//
// The percussions are the duals of the rows: lambda = -z_zero, r_c = S_c z_c (r_c = (z_c, 0, 0) without
// friction). Contacts that act on no velocity (H_c = 0) are left out and get r_c = 0.
namespace ffclib {

// Settings of the conic solvers. They go to each solver's own stopping test, which is not the error of
// error.hpp; 0 keeps the solver's default. atol and rtol set the solver's absolute and relative tolerances
// (Clarabel, ConicXX: gap_abs = atol, gap_rel = feas = rtol; SCS: eps_abs, eps_rel; QOCO: abstol, reltol;
// cuOpt: the absolute and relative primal, dual and gap tolerances).
struct ConicSettings {
    double atol{0};
    double rtol{0};
    int max_iterations{0};
    double time_limit{0};  // seconds
    bool verbose{false};   // the solver's own output
};

struct ConeDims {
    int zero{0};
    int nonneg{0};
    std::vector<int> soc;  // dimension of each second-order cone {(t, x) : ||x|| <= t}
};

struct ConeProgram {
    SparseMatrix P;  // upper triangle
    Vector q;
    SparseMatrix A;
    Vector b;
    ConeDims dims;
};

ConeProgram coneProgram(const Problem& pb);

// Whether contact c acts on a velocity (H_c != 0); the conic solvers leave out the others (r_c = 0).
bool actsOnVelocity(const Problem& pb, int c);

// The solution from the duals z of the rows of coneProgram(pb), with v = M^-1 (f + H r + G lambda).
Solution solutionFromDuals(const Problem& pb, const Vector& z);

}  // namespace ffclib
