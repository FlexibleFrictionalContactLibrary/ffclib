#pragma once

#include <vector>

#include "problem.hpp"

namespace cfc {

// Cone composition of the rows of A, in this order: zero cone, nonnegative orthant, then one
// second-order cone {(t, x) : ||x|| <= t} per entry of soc (Clarabel's and SCS's order).
struct ConeDims {
    int zero{0};
    int nonneg{0};
    std::vector<int> soc;
    int rows() const;
};

// min 1/2 x'Px + q'x  s.t.  Ax + s = b, s in K (KKT: Px + q + A'z = 0, z in K*, s'z = 0), with the upper
// triangle of P stored.
struct ConeProgram {
    SparseMatrix P;
    Vector q;
    SparseMatrix A;
    Vector b;
    ConeDims dims;
};

// The convexified problem (frictional contacts with the frozen shift `shift`) as a cone program.
// Contacts with mu = 0 are nonnegative rows on the normal component, the others second-order cones of
// the scaled variables with S_c = diag(1, mu, mu), r_c = S_c z_c. K: the bilateral rows with a nonzero
// column in C.
//
// Form V (velocity form, cones on rows; CardilloCxx's conic backends solve it):
//   x = (v, pi_K), P = blkdiag(M, C_KK), q = (-f, 0),
//   zero rows   G^T v + C_:K pi_K + s = -b,
//   contact     -S_c H_c^T v + s_c = S_c w_c + shift_c e_N,  s_c in K.
//   Solution: v = x_v, lambda = -z_zero, r_c = S_c z_c.
//
// Form F (force form, the Wolfe dual of form V with the forces as variables, cones on variables):
//   y = (v, lambda, rt), P = blkdiag(M, C, 0), q = (0, b, S_c w_c + shift_c e_N),
//   zero rows   M v - G lambda - H S rt + s = f,
//   contact     -rt_c + s_c = 0,  s_c in K.
//   Solution: v = y_v, lambda = y_lambda, r_c = S_c rt_c.
enum class Form { V, F };
const char* toString(Form form);

// Indices of the contacts in the order of the cones: frictionless (mu = 0) contacts, then frictional ones.
struct ContactOrder {
    std::vector<int> frictionless;
    std::vector<int> frictional;
};
ContactOrder contactOrder(const Problem& pb);

ConeProgram buildConeProgram(const Problem& pb, const Vector& shift, Form form);

// (v, r, lambda) from a primal-dual solution (x, z) of the cone program of `form`.
Solution recoverSolution(const Problem& pb, Form form, const Vector& x, const Vector& z);

}  // namespace cfc
