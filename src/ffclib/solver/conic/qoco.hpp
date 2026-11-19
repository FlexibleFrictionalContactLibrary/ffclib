#pragma once

#include <ffclib/solver/conic/cone_program.hpp>
#include <ffclib/solver/solver.hpp>

namespace ffclib {

// The convexified law with QOCO (interior point method, CPU) on coneProgram(pb). Built with FFCLIB_WITH_QOCO.
// QOCO has no time limit (ConicSettings::time_limit is not used). Result::converged: QOCO's status solved.
Result qoco(const Problem& pb, const ConicSettings& settings = {});

}  // namespace ffclib
