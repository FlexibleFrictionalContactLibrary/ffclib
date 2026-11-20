#pragma once

#include <string>

#include <ffclib/error.hpp>
#include <ffclib/problem.hpp>

// What all solvers share: the contact law they solve and their result.
namespace ffclib {

// Coulomb's law (with De Saxce's shift mu ||u_T||) or the convexified law (without), as in error.hpp.
enum class Law { Coulomb, Convexified };

struct Result {
    Solution solution;      // the solver's own v (Gauss-Seidel: M^-1 (f + H r + G lambda)), r, lambda
    bool converged{false};  // the solver met its stopping criterion (atol and rtol of its settings;
                            // conic solvers: in their own measures)
    std::string status;     // the solver's own description of the outcome
    int iterations{0};      // the solver's iterations (Gauss-Seidel: sweeps)
    Errors errors;          // errors(problem, solution)
    double seconds{0};      // wall time of the solve, setup included
};

}  // namespace ffclib
