#include <ffclib/solver/conic/qoco.hpp>

#include <chrono>
#include <cstdlib>
#include <string>
#include <type_traits>

extern "C" {
#include <qoco.h>
}

namespace ffclib {

// The matrices are handed to QOCO without conversion: QOCO is built with its default types.
static_assert(std::is_same_v<QOCOInt, SparseMatrix::StorageIndex>, "QOCO's index type must match SparseMatrix");
static_assert(std::is_same_v<QOCOFloat, double>, "QOCO's float type must be double");

namespace {

const char* statusName(QOCOInt exit) {
    switch (exit) {
        case QOCO_SOLVED:
            return "Solved";
        case QOCO_SOLVED_INACCURATE:
            return "SolvedInaccurate";
        case QOCO_NUMERICAL_ERROR:
            return "NumericalError";
        case QOCO_MAX_ITER:
            return "MaxIterations";
        default:
            return "Unsolved";
    }
}

}  // namespace

// QOCO solves min 1/2 x'Px + c'x s.t. Ax = b (the zero rows of the cone program), Gx <=_K h (the other rows,
// h - Gx = s in K), with the duals (y, z) of these rows: Px + c + A'y + G'z = 0, the duals of the cone program.
Result qoco(const Problem& pb, const ConicSettings& st) {
    const auto start = std::chrono::steady_clock::now();
    const ConeProgram cp = coneProgram(pb);
    const QOCOInt n = static_cast<QOCOInt>(cp.A.cols()), p = cp.dims.zero, m = static_cast<QOCOInt>(cp.A.rows()) - p;
    SparseMatrix P = cp.P, A = cp.A.topRows(p), G = cp.A.bottomRows(m);
    P.makeCompressed();
    A.makeCompressed();
    G.makeCompressed();
    Vector c = cp.q, b = cp.b.head(p), h = cp.b.tail(m);
    std::vector<QOCOInt> soc(cp.dims.soc.begin(), cp.dims.soc.end());
    QOCOCscMatrix Pq{}, Aq{}, Gq{};
    qoco_set_csc(&Pq, n, n, static_cast<QOCOInt>(P.nonZeros()), P.valuePtr(), P.outerIndexPtr(), P.innerIndexPtr());
    qoco_set_csc(&Aq, p, n, static_cast<QOCOInt>(A.nonZeros()), A.valuePtr(), A.outerIndexPtr(), A.innerIndexPtr());
    qoco_set_csc(&Gq, m, n, static_cast<QOCOInt>(G.nonZeros()), G.valuePtr(), G.outerIndexPtr(), G.innerIndexPtr());

    QOCOSettings settings;
    set_default_settings(&settings);
    settings.verbose = st.verbose ? 1 : 0;
    if (st.atol > 0) settings.abstol = st.atol;
    if (st.rtol > 0) settings.reltol = st.rtol;
    if (st.max_iterations > 0) settings.max_iters = st.max_iterations;

    Result res;
    Vector z = Vector::Zero(p + m);
    auto* solver = static_cast<QOCOSolver*>(std::malloc(sizeof(QOCOSolver)));  // freed by qoco_cleanup
    const QOCOInt setup = qoco_setup(solver, n, m, p, &Pq, c.data(), p > 0 ? &Aq : nullptr, p > 0 ? b.data() : nullptr, m > 0 ? &Gq : nullptr, m > 0 ? h.data() : nullptr, cp.dims.nonneg,
                                     static_cast<QOCOInt>(soc.size()), soc.empty() ? nullptr : soc.data(), &settings);
    if (setup != QOCO_NO_ERROR) {
        std::free(solver);
        res.status = "setup failed (QOCO error " + std::to_string(setup) + ")";
    } else {
        const QOCOInt exit = qoco_solve(solver);
        res.converged = exit == QOCO_SOLVED;
        res.status = statusName(exit);
        res.iterations = solver->sol->iters;
        if (p > 0) z.head(p) = Eigen::Map<const Vector>(solver->sol->y, p);
        if (m > 0) z.tail(m) = Eigen::Map<const Vector>(solver->sol->z, m);
        qoco_cleanup(solver);
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    res.solution = solutionFromDuals(pb, z);
    res.errors = errors(pb, res.solution.r, res.solution.lambda);
    return res;
}

}  // namespace ffclib
