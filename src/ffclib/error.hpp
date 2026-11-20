#pragma once

#include "problem.hpp"

// The errors of a candidate solution (v, r, lambda), one per condition of the problem, each evaluated at the
// candidate itself:
//
//   balance      ||M v - (f + H r + G lambda)||_inf                              momentum
//   bilateral    ||G^T v + b + C lambda||_inf                                    velocity
//   contact      max_c ||r_c - proj_K(r_c - rho_c u_hat_c)|| / rho_c             velocity
//
// with u = H^T v + w, the Euclidean norm, rho_c = 1 / lambda_max(H_c^T M^-1 H_c) per contact (contacts with
// H_c = 0 have no velocity and are skipped), and the contact error for both laws:
//
//   Coulomb's law     u_hat_c = u_c + (mu_c ||u_T,c||, 0, 0)    (De Saxce's shift)
//   convexified law   u_hat_c = u_c                             (what conic solvers solve)
//
// Each vanishes exactly where its condition holds. NaN in the candidate gives NaN in the errors that use it.
namespace ffclib {

struct Errors {
    double balance{0};      // ||M v - (f + H r + G lambda)||_inf
    double bilateral{0};    // ||G^T v + b + C lambda||_inf
    double coulomb{0};      // contacts, Coulomb's law
    double convexified{0};  // contacts, convexified law
};

// `pb` must pass check(); throws std::invalid_argument if the sizes of v, r or lambda do not match it.
Errors errors(const Problem& pb, const Solution& candidate);

}  // namespace ffclib
