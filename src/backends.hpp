#pragma once

#include <limits>
#include <string>
#include <vector>

#include <cfc/cone_program.hpp>

// The conic solvers of the benchmark, applied to one cone program
//   min 1/2 x'Px + q'x  s.t.  Ax + s = b,  s in K = {0}^z x R+^l x SOC(3)^k
// (cfc::ConeProgram, upper triangle of P). Every backend gets the same tolerance for its
// own stopping test; the results are judged by the common residuals (cfc_residuals), not by the
// solvers' own measures.
namespace cfcbench {

using VectorXr = cfc::Vector;
using cfc::ConeProgram;

struct Settings {
    double tol{1e-8};      // absolute and relative feasibility/gap tolerances of the solver
    int max_iter{500};     // IPM iterations (SCS: ADMM iterations x 100)
    double time_limit{0};  // seconds per solve, 0: none
    bool verbose{false};
};

struct Result {
    std::string status;  // the solver's own status
    bool ok{false};      // a solution was returned (solved, or inaccurate at a limit)
    bool solved{false};  // the solver reports convergence to its tolerances
    int iterations{-1};
    double t_setup{0};    // wall time of the setup (data copies, equilibration, symbolic analysis) [s]
    double t_solve{0};    // wall time of the solve [s]
    double t_solver{-1};  // solve time reported by the solver [s], -1 if none
    double objective{std::numeric_limits<double>::quiet_NaN()};
    VectorXr x, z;  // primal solution and duals of the rows (KKT: Px + q + A'z = 0)
    std::string message;
};

// Backends built into this binary: clarabel, qoco, conicxx, scs, cuopt.
std::vector<std::string> available();
std::string version(const std::string& backend);
// Whether the backend takes a program with cones on rows only (cuOpt: its Wolfe dual is built
// internally, so it is given form V).
bool acceptsForm(const std::string& backend, cfc::Form form);

// One solve, cold start. Throws std::invalid_argument for an unknown backend.
Result solve(const std::string& backend, const ConeProgram& cp, const Settings& settings);

// cuOpt only: creates the CUDA context and cuDSS handles with a tiny problem; returns the time [s]
// (excluded from all timings) or a negative value with `reason` if no GPU is usable.
double cudaWarmup(std::string* reason);

}  // namespace cfcbench
