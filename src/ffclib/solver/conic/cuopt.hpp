#pragma once

#include <ffclib/solver/conic/cone_program.hpp>
#include <ffclib/solver/solver.hpp>

namespace ffclib {

// The convexified law with NVIDIA cuOpt's barrier method (GPU). cuOpt returns no duals for problems with
// quadratic constraints, so it gets the force form, the dual of coneProgram(pb), with the percussions as
// variables, r_c = S_c t_c, S_c = diag(1, mu_c, mu_c):
//
//   min 1/2 v'Mv + 1/2 lambda'C lambda + b'lambda + sum_c (S_c w_c)'t_c
//   s.t. M v - G lambda - sum_c H_c S_c t_c = f,   t_c in SOC(3) (without friction: t_c = r_N,c >= 0)
//
// Contacts that act on no velocity are left out (r_c = 0). Built with FFCLIB_WITH_CUOPT (cuOpt >= 26.10.1).
// Result::converged: cuOpt's termination status OPTIMAL. Without a solution (e.g. no usable NVIDIA GPU),
// Result::status says why and the percussions are zero.
Result cuopt(const Problem& pb, const ConicSettings& settings = {});

}  // namespace ffclib
