#include <ffclib/solver/conic/scs.hpp>

#include <chrono>
#include <type_traits>

#include <scs.h>

namespace ffclib {

// The matrices are handed to SCS without conversion: SCS is built with its default types (DLONG and SFLOAT off).
static_assert(std::is_same_v<scs_int, SparseMatrix::StorageIndex>, "SCS's index type must match SparseMatrix");
static_assert(std::is_same_v<scs_float, double>, "SCS's float type must be double");

Result scs(const Problem& pb, const ConicSettings& st) {
    const auto start = std::chrono::steady_clock::now();
    ConeProgram cp = coneProgram(pb);
    const scs_int n = static_cast<scs_int>(cp.A.cols()), m = static_cast<scs_int>(cp.A.rows());
    ScsMatrix A{cp.A.valuePtr(), cp.A.innerIndexPtr(), cp.A.outerIndexPtr(), m, n};
    ScsMatrix P{cp.P.valuePtr(), cp.P.innerIndexPtr(), cp.P.outerIndexPtr(), n, n};
    ScsData data{m, n, &A, &P, cp.b.data(), cp.q.data()};
    std::vector<scs_int> soc(cp.dims.soc.begin(), cp.dims.soc.end());
    ScsCone cone{};  // all other cone types empty
    cone.z = cp.dims.zero;
    cone.l = cp.dims.nonneg;
    cone.q = soc.data();
    cone.qsize = static_cast<scs_int>(soc.size());

    ScsSettings settings;
    scs_set_default_settings(&settings);
    settings.verbose = st.verbose ? 1 : 0;
    if (st.atol > 0) settings.eps_abs = st.atol;
    if (st.rtol > 0) settings.eps_rel = st.rtol;
    if (st.max_iterations > 0) settings.max_iters = st.max_iterations;
    if (st.time_limit > 0) settings.time_limit_secs = st.time_limit;

    Vector x = Vector::Zero(n), y = Vector::Zero(m), s = Vector::Zero(m);
    ScsSolution sol{x.data(), y.data(), s.data()};
    ScsInfo info{};
    const scs_int flag = ::scs(&data, &cone, &settings, &sol, &info);  // SCS's C function, not this one

    Result res;
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    res.converged = flag == SCS_SOLVED;
    res.status = info.status;
    res.iterations = static_cast<int>(info.iter);
    res.solution = solutionFromConeProgram(pb, x, y);
    res.errors = errors(pb, res.solution);
    return res;
}

}  // namespace ffclib
