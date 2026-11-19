#include "backends.hpp"

#include <chrono>
#include <stdexcept>

#ifdef CFCBENCH_HAVE_CLARABEL
#include <clarabel.hpp>
#endif
#ifdef CFCBENCH_HAVE_QOCO
extern "C" {
#include <qoco.h>
}
#endif
#ifdef CFCBENCH_HAVE_CONICXX
#include <conicxx/solver.h>
#endif
#ifdef CFCBENCH_HAVE_SCS
#include "vendor/cardillo/scs_backend.hpp"
#endif
#ifdef CFCBENCH_HAVE_CUOPT
#include "vendor/cardillo/cuopt_backend.hpp"
#endif

namespace cfcbench {

namespace {

using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

#ifdef CFCBENCH_HAVE_CLARABEL
const char* clarabelStatus(clarabel::SolverStatus s) {
    switch (s) {
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

Result solveClarabel(const ConeProgram& cp, const Settings& st) {
    std::vector<clarabel::SupportedConeT<double>> cones;
    if (cp.dims.zero > 0) cones.push_back(clarabel::ZeroConeT<double>((uintptr_t)cp.dims.zero));
    if (cp.dims.nonneg > 0) cones.push_back(clarabel::NonnegativeConeT<double>((uintptr_t)cp.dims.nonneg));
    for (int d : cp.dims.soc) cones.push_back(clarabel::SecondOrderConeT<double>((uintptr_t)d));
    auto settings = clarabel::DefaultSettings<double>::default_settings();
    settings.verbose = st.verbose;
    settings.max_iter = static_cast<uint32_t>(st.max_iter);
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = st.tol;
    settings.tol_infeas_abs = settings.tol_infeas_rel = st.tol;
    if (st.time_limit > 0) settings.time_limit = st.time_limit;
    VectorXr q = cp.q, b = cp.b;  // Clarabel takes non-const references

    Result r;
    const auto t0 = Clock::now();
    clarabel::DefaultSolver<double> solver(cp.P, q, cp.A, b, cones, settings);
    const auto t1 = Clock::now();
    solver.solve();
    const auto t2 = Clock::now();
    const auto sol = solver.solution();
    r.t_setup = seconds(t0, t1);
    r.t_solve = seconds(t1, t2);
    r.t_solver = sol.solve_time;
    r.iterations = static_cast<int>(sol.iterations);
    r.status = clarabelStatus(sol.status);
    r.solved = sol.status == clarabel::SolverStatus::Solved;
    r.ok = r.solved || sol.status == clarabel::SolverStatus::AlmostSolved || sol.status == clarabel::SolverStatus::MaxIterations || sol.status == clarabel::SolverStatus::MaxTime ||
           sol.status == clarabel::SolverStatus::InsufficientProgress;
    r.objective = sol.obj_val;
    r.x = sol.x;
    r.z = sol.z;
    return r;
}
#endif

#ifdef CFCBENCH_HAVE_QOCO
// QOCO: min 1/2 x'Px + c'x s.t. Ax = b (the zero-cone rows), Gx <=_K h (the other rows: h - Gx = s in K).
// KKT Px + c + A'y + G'z = 0, so the duals of the rows are (y, z).
Result solveQoco(const ConeProgram& cp, const Settings& st) {
    const int n = static_cast<int>(cp.P.cols()), p = cp.dims.zero, m = static_cast<int>(cp.A.rows()) - p;
    cfc::SparseMatrix P = cp.P, Aeq = cp.A.topRows(p), G = cp.A.bottomRows(m);
    P.makeCompressed();
    Aeq.makeCompressed();
    G.makeCompressed();
    VectorXr c = cp.q, beq = cp.b.head(p), h = cp.b.tail(m);
    std::vector<QOCOInt> qdims(cp.dims.soc.begin(), cp.dims.soc.end());
    QOCOCscMatrix Pq{}, Aq{}, Gq{};
    qoco_set_csc(&Pq, n, n, static_cast<QOCOInt>(P.nonZeros()), P.valuePtr(), P.outerIndexPtr(), P.innerIndexPtr());
    qoco_set_csc(&Aq, p, n, static_cast<QOCOInt>(Aeq.nonZeros()), Aeq.valuePtr(), Aeq.outerIndexPtr(), Aeq.innerIndexPtr());
    qoco_set_csc(&Gq, m, n, static_cast<QOCOInt>(G.nonZeros()), G.valuePtr(), G.outerIndexPtr(), G.innerIndexPtr());

    QOCOSettings settings;
    set_default_settings(&settings);
    settings.verbose = st.verbose ? 1 : 0;
    settings.max_iters = st.max_iter;
    settings.abstol = settings.reltol = st.tol;
    // QOCO's "inaccurate" level: one order of magnitude looser than the target.
    settings.abstol_inacc = settings.reltol_inacc = 10 * st.tol;

    Result r;
    QOCOSolver* solver = static_cast<QOCOSolver*>(std::malloc(sizeof(QOCOSolver)));
    const auto t0 = Clock::now();
    const QOCOInt setup = qoco_setup(solver, n, m, p, &Pq, c.data(), p > 0 ? &Aq : nullptr, p > 0 ? beq.data() : nullptr, m > 0 ? &Gq : nullptr, m > 0 ? h.data() : nullptr, cp.dims.nonneg,
                                     static_cast<QOCOInt>(qdims.size()), qdims.empty() ? nullptr : qdims.data(), &settings);
    const auto t1 = Clock::now();
    r.t_setup = seconds(t0, t1);
    if (setup != QOCO_NO_ERROR) {
        r.status = "SetupError";
        r.message = "qoco_setup returned " + std::to_string(setup);
        std::free(solver);
        return r;
    }
    const QOCOInt exit = qoco_solve(solver);
    const auto t2 = Clock::now();
    r.t_solve = seconds(t1, t2);
    const QOCOSolution* sol = solver->sol;
    r.iterations = sol->iters;
    r.t_solver = sol->solve_time_sec;
    r.objective = sol->obj;
    r.solved = exit == QOCO_SOLVED;
    r.ok = r.solved || exit == QOCO_SOLVED_INACCURATE || exit == QOCO_MAX_ITER;
    r.status = exit == QOCO_SOLVED              ? "Solved"
               : exit == QOCO_SOLVED_INACCURATE ? "SolvedInaccurate"
               : exit == QOCO_MAX_ITER          ? "MaxIterations"
               : exit == QOCO_NUMERICAL_ERROR   ? "NumericalError"
                                                : "Unsolved";
    r.x = Eigen::Map<const VectorXr>(sol->x, n);
    r.z.resize(p + m);
    if (p > 0) r.z.head(p) = Eigen::Map<const VectorXr>(sol->y, p);
    if (m > 0) r.z.tail(m) = Eigen::Map<const VectorXr>(sol->z, m);
    qoco_cleanup(solver);
    return r;
}
#endif

#ifdef CFCBENCH_HAVE_CONICXX
Result solveConicxx(const ConeProgram& cp, const Settings& st) {
    conicxx::Settings settings;
    settings.verbose = st.verbose ? 1 : 0;
    settings.max_iter = st.max_iter;
    settings.tol_feas = settings.tol_gap_abs = settings.tol_gap_rel = st.tol;
    if (st.time_limit > 0) settings.time_limit = st.time_limit;
    settings.warm_start = false;
    conicxx::ConeSpec spec;
    spec.zero_dim = cp.dims.zero;
    spec.nonneg_dim = cp.dims.nonneg;
    spec.soc_dims.assign(cp.dims.soc.begin(), cp.dims.soc.end());

    Result r;
    conicxx::Solver solver;
    const auto t0 = Clock::now();
    const bool setup = solver.setup(cp.P, cp.q, cp.A, cp.b, spec, settings);
    const auto t1 = Clock::now();
    r.t_setup = seconds(t0, t1);
    if (!setup) {
        r.status = "SetupError";
        return r;
    }
    const conicxx::Solution& sol = solver.solve();
    const auto t2 = Clock::now();
    r.t_solve = seconds(t1, t2);
    r.t_solver = sol.info.solve_time_s;
    r.iterations = sol.info.iterations;
    r.status = conicxx::toString(sol.status);
    r.solved = sol.status == conicxx::Status::Solved;
    r.ok = r.solved || sol.status == conicxx::Status::AlmostSolved || sol.status == conicxx::Status::MaxIterations || sol.status == conicxx::Status::MaxTime ||
           sol.status == conicxx::Status::InsufficientProgress;
    r.objective = sol.objective;
    r.x = sol.x;
    r.z = sol.z;
    return r;
}
#endif

#ifdef CFCBENCH_HAVE_SCS
Result solveScs(const ConeProgram& cp, const Settings& st) {
    namespace scs = cardillo::solver::scs;
    scs::Settings settings;
    settings.eps_abs = settings.eps_rel = st.tol;
    settings.max_iters = 100 * st.max_iter;
    if (st.time_limit > 0) settings.time_limit_secs = st.time_limit;
    settings.verbose = st.verbose;
    scs::ConeDims dims;
    dims.zero = cp.dims.zero;
    dims.nonneg = cp.dims.nonneg;
    dims.soc = cp.dims.soc;

    Result r;
    scs::Solver solver;
    const auto t0 = Clock::now();
    const bool setup = solver.setup(cp.P, cp.q, cp.A, cp.b, dims, settings);
    const auto t1 = Clock::now();
    r.t_setup = seconds(t0, t1);
    if (!setup) {
        r.status = "SetupError";
        return r;
    }
    const scs::Status status = solver.solve(false);
    const auto t2 = Clock::now();
    r.t_solve = seconds(t1, t2);
    r.t_solver = solver.info().solve_ms * 1e-3;
    r.iterations = solver.info().iterations;
    r.status = scs::toString(status);
    r.solved = status == scs::Status::Solved;
    r.ok = scs::isAcceptable(status) || status == scs::Status::Unknown;
    r.objective = solver.info().pobj;
    r.x = solver.x();
    r.z = solver.y();
    return r;
}
#endif

#ifdef CFCBENCH_HAVE_CUOPT
Result solveCuopt(const ConeProgram& cp, const Settings& st) {
    namespace cu = cardillo::solver::cuopt;
    cu::Settings settings;
    settings.tol_abs = settings.tol_rel = st.tol;
    settings.iteration_limit = st.max_iter;
    if (st.time_limit > 0) settings.time_limit = st.time_limit;
    settings.verbose = st.verbose;
    cu::ConeDims dims;
    dims.zero = cp.dims.zero;
    dims.nonneg = cp.dims.nonneg;
    dims.soc = cp.dims.soc;

    Result r;
    cu::Solver solver(settings);
    const auto t0 = Clock::now();
    solver.setup(cp.P, cp.q, cp.A, cp.b, dims);  // Wolfe dual in cuOpt's form (CPU) and cuOpt's problem
    const auto t1 = Clock::now();
    const cu::Status status = solver.solve();  // transfer, presolve, barrier on the GPU
    const auto t2 = Clock::now();
    r.t_setup = seconds(t0, t1);
    r.t_solve = seconds(t1, t2);
    const cu::Info& info = solver.info();
    r.t_solver = info.solve_time;
    r.iterations = info.iterations;
    r.status = info.termination >= 0 ? std::string(cu::toString(status)) + " (termination " + std::to_string(info.termination) + ")" : cu::toString(status);
    r.solved = status == cu::Status::Optimal;
    r.ok = status != cu::Status::Failed;
    r.message = info.message;
    r.objective = -info.objective;  // the Wolfe dual's optimum is minus the cone program's
    if (r.ok) {
        solver.recover(cp.A, cp.b);
        r.x = solver.x();
        r.z = solver.z();
    }
    return r;
}
#endif

}  // namespace

std::vector<std::string> available() {
    std::vector<std::string> b;
#ifdef CFCBENCH_HAVE_CLARABEL
    b.push_back("clarabel");
#endif
#ifdef CFCBENCH_HAVE_QOCO
    b.push_back("qoco");
#endif
#ifdef CFCBENCH_HAVE_CONICXX
    b.push_back("conicxx");
#endif
#ifdef CFCBENCH_HAVE_SCS
    b.push_back("scs");
#endif
#ifdef CFCBENCH_HAVE_CUOPT
    b.push_back("cuopt");
#endif
    return b;
}

std::string version(const std::string& backend) {
#ifdef CFCBENCH_HAVE_CUOPT
    if (backend == "cuopt") return cardillo::solver::cuopt::Solver::version();
#endif
    if (backend == "clarabel") return CFCBENCH_CLARABEL_VERSION;
    if (backend == "qoco") return CFCBENCH_QOCO_VERSION;
    if (backend == "scs") return CFCBENCH_SCS_VERSION;
    if (backend == "conicxx") return CFCBENCH_CONICXX_VERSION;
    return "unknown";
}

bool acceptsForm(const std::string& backend, cfc::Form form) {
    return backend != "cuopt" || form == cfc::Form::V;
}

Result solve(const std::string& backend, const ConeProgram& cp, const Settings& settings) {
#ifdef CFCBENCH_HAVE_CLARABEL
    if (backend == "clarabel") return solveClarabel(cp, settings);
#endif
#ifdef CFCBENCH_HAVE_QOCO
    if (backend == "qoco") return solveQoco(cp, settings);
#endif
#ifdef CFCBENCH_HAVE_CONICXX
    if (backend == "conicxx") return solveConicxx(cp, settings);
#endif
#ifdef CFCBENCH_HAVE_SCS
    if (backend == "scs") return solveScs(cp, settings);
#endif
#ifdef CFCBENCH_HAVE_CUOPT
    if (backend == "cuopt") return solveCuopt(cp, settings);
#endif
    throw std::invalid_argument("unknown or not built backend '" + backend + "'");
}

double cudaWarmup(std::string* reason) {
#ifdef CFCBENCH_HAVE_CUOPT
    const auto t0 = Clock::now();
    if (!cardillo::solver::cuopt::deviceAvailable(reason)) return -1;
    return seconds(t0, Clock::now());
#else
    if (reason) *reason = "cfc_solve was built without cuOpt";
    return -1;
#endif
}

}  // namespace cfcbench
