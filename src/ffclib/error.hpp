#pragma once

#include "problem.hpp"

// The two error measures of a candidate (r, lambda), in velocity units. Both are evaluated at the velocity
// that the percussions produce,
//
//   v = M^-1 (f + H r + G lambda),
//
// so that the balance holds by construction and a solver needs to return only (r, lambda):
//
//   error = max( ||G^T v + b + C lambda||_inf,  max_c ||r_c - proj_K(r_c - rho_c u_hat_c)|| / rho_c )
//
// with u = H^T v + w, the Euclidean norm, rho_c = 1 / lambda_max(H_c^T M^-1 H_c) per contact (contacts with
// H_c = 0 have no velocity and are skipped), and
//
//   Coulomb's law     u_hat_c = u_c + (mu_c ||u_T,c||, 0, 0)    (De Saxce's shift)
//   convexified law   u_hat_c = u_c                             (what conic solvers solve)
//
// Both vanish exactly at the solutions of their law. NaN in (r, lambda) gives NaN.
namespace ffclib {

struct Errors {
    double coulomb{0};
    double convexified{0};
};

// `pb` must pass check(); throws std::invalid_argument if the sizes of r or lambda do not match it.
Errors errors(const Problem& pb, const Vector& r, const Vector& lambda);

}  // namespace ffclib
