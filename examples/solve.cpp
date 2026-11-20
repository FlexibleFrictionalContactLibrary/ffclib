// Reads a problem, solves it and prints the errors of the solution (and of the solution stored in the file, if
// any): balance (momentum units), bilateral rows and contacts with both laws (velocity units), error.hpp.
//
//   ffclib_solve FILE [--solver NAME] [--atol A] [--rtol R] [--max-iterations N] [--out FILE]
//                     [--law coulomb|convexified] [--relaxation W] [--cold]
//
// NAME: gauss-seidel (default) or a conic solver that was built (clarabel, scs, qoco, conicxx, cuopt); an
// unknown name lists the built ones. --atol, --rtol and --max-iterations go to the solver's own stopping test
// (Gauss-Seidel: the natural map of its law and the sweeps). Gauss-Seidel only: --law, --relaxation, and --cold to
// start from zero instead of the file's first guess. The conic solvers solve the convexified law. --out writes
// the problem with the new solution.

#include <cstdio>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <ffclib/error.hpp>
#include <ffclib/io.hpp>
#include <ffclib/solver/conic/cone_program.hpp>
#include <ffclib/solver/gauss_seidel.hpp>
#ifdef FFCLIB_WITH_CLARABEL
#include <ffclib/solver/conic/clarabel.hpp>
#endif
#ifdef FFCLIB_WITH_SCS
#include <ffclib/solver/conic/scs.hpp>
#endif
#ifdef FFCLIB_WITH_QOCO
#include <ffclib/solver/conic/qoco.hpp>
#endif
#ifdef FFCLIB_WITH_CONICXX
#include <ffclib/solver/conic/conicxx.hpp>
#endif
#ifdef FFCLIB_WITH_CUOPT
#include <ffclib/solver/conic/cuopt.hpp>
#endif

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s FILE [--solver NAME] [--atol A] [--rtol R] [--max-iterations N] [--out FILE] [--law coulomb|convexified] [--relaxation W] [--cold]\n", argv[0]);
        return 2;
    }
    std::string name = "gauss-seidel", out;
    ffclib::GaussSeidelSettings gs;
    ffclib::ConicSettings conic;
    bool cold = false;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        const char* value = i + 1 < argc ? argv[i + 1] : "";
        if (arg == "--cold") {
            cold = true;
            continue;
        }
        if (arg == "--solver")
            name = value;
        else if (arg == "--atol")
            gs.atol = conic.atol = std::stod(value);
        else if (arg == "--rtol")
            gs.rtol = conic.rtol = std::stod(value);
        else if (arg == "--max-iterations")
            gs.max_sweeps = conic.max_iterations = std::stoi(value);
        else if (arg == "--out")
            out = value;
        else if (arg == "--law")
            gs.law = std::string(value) == "convexified" ? ffclib::Law::Convexified : ffclib::Law::Coulomb;
        else if (arg == "--relaxation")
            gs.relaxation = std::stod(value);
        else {
            std::fprintf(stderr, "unknown option %s\n", arg.c_str());
            return 2;
        }
        ++i;
    }

    // The solvers of this build: a problem and a starting point (used by Gauss-Seidel only) to a result.
    using Solve = std::function<ffclib::Result(const ffclib::Problem&, const ffclib::Solution*)>;
    const std::vector<std::pair<std::string, Solve>> solvers = {
        {"gauss-seidel", [&](const ffclib::Problem& pb, const ffclib::Solution* guess) { return ffclib::gaussSeidel(pb, gs, guess); }},
#ifdef FFCLIB_WITH_CLARABEL
        {"clarabel", [&](const ffclib::Problem& pb, const ffclib::Solution*) { return ffclib::clarabel(pb, conic); }},
#endif
#ifdef FFCLIB_WITH_SCS
        {"scs", [&](const ffclib::Problem& pb, const ffclib::Solution*) { return ffclib::scs(pb, conic); }},
#endif
#ifdef FFCLIB_WITH_QOCO
        {"qoco", [&](const ffclib::Problem& pb, const ffclib::Solution*) { return ffclib::qoco(pb, conic); }},
#endif
#ifdef FFCLIB_WITH_CONICXX
        {"conicxx", [&](const ffclib::Problem& pb, const ffclib::Solution*) { return ffclib::conicxx(pb, conic); }},
#endif
#ifdef FFCLIB_WITH_CUOPT
        {"cuopt", [&](const ffclib::Problem& pb, const ffclib::Solution*) { return ffclib::cuopt(pb, conic); }},
#endif
    };
    const Solve* solve = nullptr;
    for (const auto& [solverName, function] : solvers)
        if (solverName == name) solve = &function;
    if (!solve) {
        std::fprintf(stderr, "unknown solver %s; built:", name.c_str());
        for (const auto& solver : solvers) std::fprintf(stderr, " %s", solver.first.c_str());
        std::fprintf(stderr, "\n");
        return 2;
    }

    auto print = [](const char* what, const ffclib::Errors& e) {
        std::printf("%-16s balance %.3e, bilateral %.3e, contacts Coulomb %.3e, convexified %.3e\n", what, e.balance, e.bilateral, e.coulomb, e.convexified);
    };
    const ffclib::Record rec = ffclib::read(argv[1]);
    const ffclib::Problem& pb = rec.problem;
    std::printf("%s: nv %d, nc %d, nb %d\n", argv[1], pb.nv(), pb.nc(), pb.nb());
    if (rec.solution) {
        print("stored solution", ffclib::errors(pb, *rec.solution));
    }

    const ffclib::Result res = (*solve)(pb, !cold && !rec.guesses.empty() ? &rec.guesses.front() : nullptr);
    std::printf("%s: %s, %d iterations, %.3f s\n", name.c_str(), res.status.c_str(), res.iterations, res.seconds);
    print("solution", res.errors);

    if (!out.empty()) {
        ffclib::Record solved = rec;
        solved.solution = res.solution;
        ffclib::write(out, solved);
    }
    return res.converged ? 0 : 1;
}
