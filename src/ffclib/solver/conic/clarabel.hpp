#pragma once

#include <ffclib/solver/conic/cone_program.hpp>
#include <ffclib/solver/solver.hpp>

namespace ffclib {

// The convexified law with Clarabel (interior point method, CPU) on coneProgram(pb). Built with
// FFCLIB_WITH_CLARABEL. Result::converged: Clarabel's status Solved.
Result clarabel(const Problem& pb, const ConicSettings& settings = {});

}  // namespace ffclib
