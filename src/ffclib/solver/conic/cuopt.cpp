#include <ffclib/solver/conic/cuopt.hpp>

#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <cuopt/mathematical_optimization/cuopt_c.h>

namespace ffclib {

// The arrays are handed to cuOpt without conversion: cuOpt's default build has 32-bit indices and doubles.
static_assert(std::is_same_v<cuopt_int_t, SparseMatrix::StorageIndex>, "cuOpt's index type must match SparseMatrix");
static_assert(std::is_same_v<cuopt_float_t, double>, "cuOpt's float type must be double");

namespace {

// Owns cuOpt's handles; its destroy functions take the address of a handle.
struct Handles {
    cuOptSolverSettings settings{nullptr};
    cuOptOptimizationProblem problem{nullptr};
    cuOptSolution solution{nullptr};

    Handles() = default;
    Handles(const Handles&) = delete;
    Handles& operator=(const Handles&) = delete;
    ~Handles() {
        if (solution) cuOptDestroySolution(&solution);
        if (problem) cuOptDestroyProblem(&problem);
        if (settings) cuOptDestroySolverSettings(&settings);
    }
};

void require(cuopt_int_t code, const char* call) {
    if (code != CUOPT_SUCCESS) throw std::runtime_error(std::string("ffclib cuOpt: ") + call + " returned " + std::to_string(code));
}

const char* terminationName(cuopt_int_t status) {
    switch (status) {
        case CUOPT_TERMINATION_STATUS_NO_TERMINATION:
            return "NO_TERMINATION";
        case CUOPT_TERMINATION_STATUS_OPTIMAL:
            return "OPTIMAL";
        case CUOPT_TERMINATION_STATUS_INFEASIBLE:
            return "INFEASIBLE";
        case CUOPT_TERMINATION_STATUS_UNBOUNDED:
            return "UNBOUNDED";
        case CUOPT_TERMINATION_STATUS_ITERATION_LIMIT:
            return "ITERATION_LIMIT";
        case CUOPT_TERMINATION_STATUS_TIME_LIMIT:
            return "TIME_LIMIT";
        case CUOPT_TERMINATION_STATUS_NUMERICAL_ERROR:
            return "NUMERICAL_ERROR";
        case CUOPT_TERMINATION_STATUS_PRIMAL_FEASIBLE:
            return "PRIMAL_FEASIBLE";
        case CUOPT_TERMINATION_STATUS_FEASIBLE_FOUND:
            return "FEASIBLE_FOUND";
        case CUOPT_TERMINATION_STATUS_CONCURRENT_LIMIT:
            return "CONCURRENT_LIMIT";
        case CUOPT_TERMINATION_STATUS_WORK_LIMIT:
            return "WORK_LIMIT";
        case CUOPT_TERMINATION_STATUS_UNBOUNDED_OR_INFEASIBLE:
            return "UNBOUNDED_OR_INFEASIBLE";
        default:
            return "unknown termination status";
    }
}

}  // namespace

Result cuopt(const Problem& pb, const ConicSettings& st) {
    const auto start = std::chrono::steady_clock::now();
    const int nv = pb.nv(), nb = pb.nb(), nc = pb.nc();

    // Variables (v, lambda, t). first[c]: the first variable of t_c (one without friction, three with), -1 if
    // the contact is left out.
    std::vector<int> first(static_cast<std::size_t>(nc), -1);
    int nx = nv + nb;
    for (int c = 0; c < nc; ++c)
        if (actsOnVelocity(pb, c)) {
            first[c] = nx;
            nx += pb.mu[c] > 0 ? 3 : 1;
        }
    auto components = [&](int c) { return pb.mu[c] > 0 ? 3 : 1; };
    auto scale = [&](int c, int j) { return j == 0 ? 1.0 : pb.mu[c]; };

    // Rows M v - G lambda - H_c S_c t_c = f, in compressed rows.
    std::vector<Eigen::Triplet<double>> t;
    for (int j = 0; j < nv; ++j)
        for (SparseMatrix::InnerIterator it(pb.M, j); it; ++it) t.emplace_back(static_cast<int>(it.row()), j, it.value());
    for (int k = 0; k < nb; ++k)
        for (SparseMatrix::InnerIterator it(pb.G, k); it; ++it) t.emplace_back(static_cast<int>(it.row()), nv + k, -it.value());
    for (int c = 0; c < nc; ++c)
        for (int j = 0; first[c] >= 0 && j < components(c); ++j)
            for (SparseMatrix::InnerIterator it(pb.H, 3 * c + j); it; ++it) t.emplace_back(static_cast<int>(it.row()), first[c] + j, -scale(c, j) * it.value());
    Eigen::SparseMatrix<double, Eigen::RowMajor, int> rows(nv, nx);
    rows.setFromTriplets(t.begin(), t.end());
    rows.makeCompressed();

    // Objective x'Qx + c'x (cuOpt's form has no factor 1/2): Q = blkdiag(M, C) / 2 over both triangles,
    // c = (0, b, S_c w_c). Bounds: t_N,c >= 0, everything else free.
    std::vector<int> qRows, qCols;
    std::vector<double> qValues;
    auto quadratic = [&](const SparseMatrix& A, int offset) {
        for (int j = 0; j < A.outerSize(); ++j)
            for (SparseMatrix::InnerIterator it(A, j); it; ++it)
                if (it.value() != 0) {
                    qRows.push_back(offset + static_cast<int>(it.row()));
                    qCols.push_back(offset + j);
                    qValues.push_back(0.5 * it.value());
                }
    };
    quadratic(pb.M, 0);
    quadratic(pb.C, nv);
    std::vector<double> cost(static_cast<std::size_t>(nx), 0.0);
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> lower(static_cast<std::size_t>(nx), -inf), upper(static_cast<std::size_t>(nx), inf);
    for (int k = 0; k < nb; ++k) cost[nv + k] = pb.b[k];
    for (int c = 0; c < nc; ++c) {
        if (first[c] < 0) continue;
        for (int j = 0; j < components(c); ++j) cost[first[c] + j] = scale(c, j) * pb.w[3 * c + j];
        lower[first[c]] = 0;
    }
    std::vector<double> rhs(pb.f.data(), pb.f.data() + nv);
    std::vector<char> senses(static_cast<std::size_t>(nv), CUOPT_EQUAL), types(static_cast<std::size_t>(nx), CUOPT_CONTINUOUS);

    Handles h;
    require(cuOptCreateSolverSettings(&h.settings), "cuOptCreateSolverSettings");
    auto setInt = [&](const char* name, cuopt_int_t value) { require(cuOptSetIntegerParameter(h.settings, name, value), name); };
    auto setFloat = [&](const char* name, cuopt_float_t value) { require(cuOptSetFloatParameter(h.settings, name, value), name); };
    setInt(CUOPT_METHOD, CUOPT_METHOD_BARRIER);  // the method for quadratic constraints
    setInt(CUOPT_LOG_TO_CONSOLE, st.verbose ? 1 : 0);
    if (st.atol > 0)
        for (const char* name : {CUOPT_ABSOLUTE_PRIMAL_TOLERANCE, CUOPT_ABSOLUTE_DUAL_TOLERANCE, CUOPT_ABSOLUTE_GAP_TOLERANCE}) setFloat(name, st.atol);
    if (st.rtol > 0)
        for (const char* name : {CUOPT_RELATIVE_PRIMAL_TOLERANCE, CUOPT_RELATIVE_DUAL_TOLERANCE, CUOPT_RELATIVE_GAP_TOLERANCE}) setFloat(name, st.rtol);
    if (st.max_iterations > 0) setInt(CUOPT_ITERATION_LIMIT, st.max_iterations);
    if (st.time_limit > 0) setFloat(CUOPT_TIME_LIMIT, st.time_limit);

    require(cuOptCreateProblem(nv, nx, CUOPT_MINIMIZE, 0.0, cost.data(), rows.outerIndexPtr(), rows.innerIndexPtr(), rows.valuePtr(), senses.data(), rhs.data(), lower.data(), upper.data(),
                               types.data(), &h.problem),
            "cuOptCreateProblem");
    require(cuOptSetQuadraticObjective(h.problem, static_cast<cuopt_int_t>(qValues.size()), qRows.data(), qCols.data(), qValues.data()), "cuOptSetQuadraticObjective");
    // t_T1^2 + t_T2^2 - t_N^2 <= 0 with t_N >= 0: t_c in SOC(3).
    for (int c = 0; c < nc; ++c) {
        if (first[c] < 0 || components(c) == 1) continue;
        const cuopt_int_t vars[3] = {first[c], first[c] + 1, first[c] + 2};
        const cuopt_float_t coefficients[3] = {-1, 1, 1};
        require(cuOptAddQuadraticConstraint(h.problem, 3, vars, vars, coefficients, 0, nullptr, nullptr, CUOPT_LESS_THAN, 0.0), "cuOptAddQuadraticConstraint");
    }

    Result res;
    Vector x = Vector::Zero(nx);
    const cuopt_int_t rc = cuOptSolve(h.problem, h.settings, &h.solution);
    if (!h.solution || rc != CUOPT_SUCCESS) {
        res.status = "no solution: cuOptSolve returned " + std::to_string(rc);
        if (rc == CUOPT_VALIDATION_ERROR) res.status += " (cuOpt reports this also without a usable NVIDIA GPU and driver)";
    } else {
        cuopt_int_t termination = -1;
        require(cuOptGetTerminationStatus(h.solution, &termination), "cuOptGetTerminationStatus");
        res.status = terminationName(termination);
        res.converged = termination == CUOPT_TERMINATION_STATUS_OPTIMAL;
        cuopt_int_t iterations = 0;
        if (cuOptGetSolutionIntAttribute(h.solution, CUOPT_SOLUTION_ATTR_LP_NUM_ITERATIONS, &iterations) == CUOPT_SUCCESS) res.iterations = iterations;
        require(cuOptGetPrimalSolution(h.solution, x.data()), "cuOptGetPrimalSolution");
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    const Vector lambda = x.segment(nv, nb);
    Vector r = Vector::Zero(3 * nc);
    for (int c = 0; c < nc; ++c)
        for (int j = 0; first[c] >= 0 && j < components(c); ++j) r[3 * c + j] = scale(c, j) * x[first[c] + j];
    res.solution = {velocity(pb, r, lambda), r, lambda};
    res.errors = errors(pb, r, lambda);
    return res;
}

}  // namespace ffclib
