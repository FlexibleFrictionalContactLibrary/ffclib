#pragma once

#include <ffclib/solver/conic/cone_program.hpp>
#include <ffclib/solver/solver.hpp>

namespace ffclib {

// The convexified law with ConicXX (interior point method, CPU) on coneProgram(pb). Built with
// FFCLIB_WITH_CONICXX. Result::converged: ConicXX's status Solved.
Result conicxx(const Problem& pb, const ConicSettings& settings = {});

}  // namespace ffclib
