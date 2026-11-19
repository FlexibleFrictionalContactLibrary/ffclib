#pragma once

#include <ffclib/solver/conic/cone_program.hpp>
#include <ffclib/solver/solver.hpp>

namespace ffclib {

// The convexified law with SCS (operator splitting, CPU, direct linear-system solver) on coneProgram(pb).
// Built with FFCLIB_WITH_SCS. Result::converged: SCS's status solved.
Result scs(const Problem& pb, const ConicSettings& settings = {});

}  // namespace ffclib
