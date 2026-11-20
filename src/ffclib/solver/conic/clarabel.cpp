#include <ffclib/solver/conic/clarabel.hpp>

#include <chrono>
#include <cstdint>

#include <clarabel.hpp>

namespace ffclib {

namespace {

const char* statusName(clarabel::SolverStatus status) {
    switch (status) {
        case clarabel::SolverStatus::Solved:
            return "Solved";
        case clarabel::SolverStatus::AlmostSolved:
            return "AlmostSolved";
        case clarabel::SolverStatus::PrimalInfeasible:
            return "PrimalInfeasible";
        case clarabel::SolverStatus::DualInfeasible:
            return "DualInfeasible";
        case clarabel::SolverStatus::AlmostPrimalInfeasible:
            return "AlmostPrimalInfeasible";
        case clarabel::SolverStatus::AlmostDualInfeasible:
            return "AlmostDualInfeasible";
        case clarabel::SolverStatus::MaxIterations:
            return "MaxIterations";
        case clarabel::SolverStatus::MaxTime:
            return "MaxTime";
        case clarabel::SolverStatus::NumericalError:
            return "NumericalError";
        case clarabel::SolverStatus::InsufficientProgress:
            return "InsufficientProgress";
        default:
            return "Unsolved";
    }
}

}  // namespace

Result clarabel(const Problem& pb, const ConicSettings& st) {
    const auto start = std::chrono::steady_clock::now();
    const ConeProgram cp = coneProgram(pb);
    std::vector<clarabel::SupportedConeT<double>> cones;
    if (cp.dims.zero > 0) cones.push_back(clarabel::ZeroConeT<double>(static_cast<std::uintptr_t>(cp.dims.zero)));
    if (cp.dims.nonneg > 0) cones.push_back(clarabel::NonnegativeConeT<double>(static_cast<std::uintptr_t>(cp.dims.nonneg)));
    for (const int d : cp.dims.soc) cones.push_back(clarabel::SecondOrderConeT<double>(static_cast<std::uintptr_t>(d)));

    auto settings = clarabel::DefaultSettings<double>::default_settings();
    settings.verbose = st.verbose;
    if (st.atol > 0) settings.tol_gap_abs = st.atol;
    if (st.rtol > 0) settings.tol_gap_rel = settings.tol_feas = st.rtol;
    if (st.max_iterations > 0) settings.max_iter = static_cast<std::uint32_t>(st.max_iterations);
    if (st.time_limit > 0) settings.time_limit = st.time_limit;

    Vector q = cp.q, b = cp.b;  // Clarabel takes non-const references
    clarabel::DefaultSolver<double> solver(cp.P, q, cp.A, b, cones, settings);
    solver.solve();
    const auto sol = solver.solution();

    Result res;
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    res.converged = sol.status == clarabel::SolverStatus::Solved;
    res.status = statusName(sol.status);
    res.iterations = static_cast<int>(sol.iterations);
    res.solution = solutionFromConeProgram(pb, sol.x, sol.z);
    res.errors = errors(pb, res.solution);
    return res;
}

}  // namespace ffclib
