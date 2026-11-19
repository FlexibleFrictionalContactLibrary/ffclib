#include <ffclib/solver/conic/conicxx.hpp>

#include <chrono>

#include <conicxx/solver.h>

namespace ffclib {

Result conicxx(const Problem& pb, const ConicSettings& st) {
    const auto start = std::chrono::steady_clock::now();
    const ConeProgram cp = coneProgram(pb);
    conicxx::Settings settings;
    settings.verbose = st.verbose ? 1 : 0;
    if (st.atol > 0) settings.tol_gap_abs = st.atol;
    if (st.rtol > 0) settings.tol_gap_rel = settings.tol_feas = st.rtol;
    if (st.max_iterations > 0) settings.max_iter = st.max_iterations;
    if (st.time_limit > 0) settings.time_limit = st.time_limit;
    conicxx::ConeSpec spec;
    spec.zero_dim = cp.dims.zero;
    spec.nonneg_dim = cp.dims.nonneg;
    spec.soc_dims.assign(cp.dims.soc.begin(), cp.dims.soc.end());

    Result res;
    Vector z = Vector::Zero(cp.A.rows());
    conicxx::Solver solver;
    if (!solver.setup(cp.P, cp.q, cp.A, cp.b, spec, settings)) {
        res.status = "setup failed";
    } else {
        const conicxx::Solution& sol = solver.solve();
        res.converged = sol.status == conicxx::Status::Solved;
        res.status = conicxx::toString(sol.status);
        res.iterations = sol.info.iterations;
        z = sol.z;
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    res.solution = solutionFromDuals(pb, z);
    res.errors = errors(pb, res.solution.r, res.solution.lambda);
    return res;
}

}  // namespace ffclib
