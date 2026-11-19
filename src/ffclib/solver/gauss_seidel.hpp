#pragma once

#include <ffclib/solver/solver.hpp>

// Naive projected Gauss-Seidel / SOR on the global problem. It keeps v = M^-1 (f + H r + G lambda) and, in
// every sweep, visits each bilateral row k and then each contact c once, updating v after every change:
//
//   row k:      lambda_k <- lambda_k - omega (G_k^T v + b_k + (C lambda)_k) / (G_k^T M^-1 G_k + C_kk)
//   contact c:  r_c      <- proj_K(r_c - omega rho_c u_hat_c),   rho_c = 1 / lambda_max(H_c^T M^-1 H_c)
//
// with u_hat_c = H_c^T v + w_c, plus De Saxce's shift (mu_c ||u_T,c||, 0, 0) of the current velocity for
// Coulomb's law. For the convexified law every update decreases a convex function of (r, lambda) for
// 0 < omega < 2; for Coulomb's law convergence is not guaranteed. Rows and contacts that act on no velocity
// and have no compliance keep their starting value.
//
// Stopping test (as CardilloCxx's projection solvers): a sweep converges if every update in it met
//
//   row k:      |G_k^T v + b_k + (C lambda)_k| <= atol + rtol |G_k^T v + b_k|
//   contact c:  |d_N| <= atol + rtol |u_N|,  ||d_T|| <= atol + rtol ||u_T||,   d = (r_c - r_c^new) / (omega rho_c)
//
// at the state of the update, u = H_c^T v + w_c the relative velocity (without the shift). d is the natural map
// of the law at that state with the step omega rho_c, in velocity units; normal and tangential parts are scaled
// separately, so that a normal error is not hidden behind a large slip. Result::errors are evaluated afterwards.
namespace ffclib {

struct GaussSeidelSettings {
    Law law{Law::Coulomb};
    double atol{1e-6};  // absolute tolerance, velocity units; > 0
    double rtol{0};     // relative tolerance; >= 0
    int max_sweeps{10000};
    double relaxation{1};  // omega: 1 Gauss-Seidel, between 1 and 2 over-relaxation (SOR)
};

// Starts from the percussions of `guess` (its velocity is not used) or from zero. Throws
// std::invalid_argument for invalid settings or a guess that does not match the problem.
Result gaussSeidel(const Problem& pb, const GaussSeidelSettings& settings, const Solution* guess = nullptr);

}  // namespace ffclib
