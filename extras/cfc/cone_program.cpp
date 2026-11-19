#include <cfc/cone_program.hpp>

#include <numeric>
#include <stdexcept>

namespace cfc {

namespace {

using Triplets = std::vector<Eigen::Triplet<double>>;

SparseMatrix fromTriplets(Eigen::Index rows, Eigen::Index cols, const Triplets& t) {
    SparseMatrix A(rows, cols);
    A.setFromTriplets(t.begin(), t.end());
    A.makeCompressed();
    return A;
}

// Upper triangle of the symmetric matrix M at offset (k, k).
void addUpper(Triplets& t, const SparseMatrix& M, int k) {
    for (Eigen::Index j = 0; j < M.outerSize(); ++j)
        for (SparseMatrix::InnerIterator it(M, j); it; ++it)
            if (it.row() <= j) t.emplace_back(k + static_cast<int>(it.row()), k + static_cast<int>(j), it.value());
}

// Scale S = (1, mu, mu) of the components of contact c.
double componentScale(const Problem& pb, int c, int j) {
    return j == 0 ? 1.0 : pb.mu[c];
}

// Number of cone rows (and scaled force variables) of the contacts: 1 per frictionless, 3 per frictional.
int contactRowCount(const ContactOrder& order) {
    return static_cast<int>(order.frictionless.size() + 3 * order.frictional.size());
}

// Components (contact, j) of the cone rows in their order.
std::vector<std::pair<int, int>> coneComponents(const ContactOrder& order) {
    std::vector<std::pair<int, int>> comp;
    comp.reserve(static_cast<std::size_t>(contactRowCount(order)));
    for (int c : order.frictionless) comp.emplace_back(c, 0);
    for (int c : order.frictional)
        for (int j = 0; j < 3; ++j) comp.emplace_back(c, j);
    return comp;
}

ConeDims coneDims(int zero, const ContactOrder& order) {
    ConeDims dims;
    dims.zero = zero;
    dims.nonneg = static_cast<int>(order.frictionless.size());
    dims.soc.assign(order.frictional.size(), 3);
    return dims;
}

}  // namespace

int ConeDims::rows() const {
    return zero + nonneg + std::accumulate(soc.begin(), soc.end(), 0);
}

const char* toString(Form form) {
    return form == Form::V ? "V" : "F";
}

ContactOrder contactOrder(const Problem& pb) {
    ContactOrder order;
    for (int c = 0; c < pb.nc(); ++c) (pb.mu[c] > 0 ? order.frictional : order.frictionless).push_back(c);
    return order;
}

ConeProgram buildConeProgram(const Problem& pb, const Vector& shift, Form form) {
    const int n = pb.n(), p = pb.p();
    if (shift.size() != pb.nc()) throw std::invalid_argument("cfc cone program: the shift does not match the contacts");
    validate(pb);
    const ContactOrder order = contactOrder(pb);
    const auto comp = coneComponents(order);
    const int mc = static_cast<int>(comp.size());
    const SparseMatrix Gt = pb.G.transpose();

    ConeProgram cp;
    Triplets P, A;
    if (form == Form::V) {
        // Variables (v, pi_K): one auxiliary variable per bilateral row with a nonzero column in C.
        std::vector<int> aux(static_cast<std::size_t>(p), -1);
        int na = 0;
        for (int k = 0; k < p; ++k)
            for (SparseMatrix::InnerIterator it(pb.C, k); it; ++it)
                if (it.value() != 0) {
                    aux[(size_t)k] = na++;
                    break;
                }
        const int nx = n + na;
        addUpper(P, pb.M, 0);
        for (int k = 0; k < pb.C.outerSize(); ++k)
            for (SparseMatrix::InnerIterator it(pb.C, k); it; ++it)
                if (it.value() != 0 && it.row() <= k) P.emplace_back(n + aux[(size_t)it.row()], n + aux[(size_t)k], it.value());
        cp.q = Vector::Zero(nx);
        cp.q.head(n) = -pb.f;

        cp.b.resize(p + mc);
        for (int k = 0; k < Gt.outerSize(); ++k)
            for (SparseMatrix::InnerIterator it(Gt, k); it; ++it) A.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
        for (int k = 0; k < pb.C.outerSize(); ++k)
            for (SparseMatrix::InnerIterator it(pb.C, k); it; ++it)
                if (it.value() != 0) A.emplace_back(static_cast<int>(it.row()), n + aux[(size_t)k], it.value());
        cp.b.head(p) = -pb.b;
        for (int row = 0; row < mc; ++row) {
            const auto [c, j] = comp[(size_t)row];
            const double s = componentScale(pb, c, j);
            for (SparseMatrix::InnerIterator it(pb.H, 3 * c + j); it; ++it) A.emplace_back(p + row, static_cast<int>(it.row()), -s * it.value());
            cp.b[p + row] = s * pb.w[3 * c + j] + (j == 0 ? shift[c] : 0.0);
        }
        cp.P = fromTriplets(nx, nx, P);
        cp.A = fromTriplets(p + mc, nx, A);
        cp.dims = coneDims(p, order);
        return cp;
    }

    // Form F: variables (v, lambda, rt).
    const int ny = n + p + mc;
    addUpper(P, pb.M, 0);
    for (int k = 0; k < pb.C.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(pb.C, k); it; ++it)
            if (it.value() != 0 && it.row() <= k) P.emplace_back(n + static_cast<int>(it.row()), n + k, it.value());
    cp.q = Vector::Zero(ny);
    cp.q.segment(n, p) = pb.b;
    for (Eigen::Index j = 0; j < pb.M.outerSize(); ++j)
        for (SparseMatrix::InnerIterator it(pb.M, j); it; ++it) A.emplace_back(static_cast<int>(it.row()), static_cast<int>(j), it.value());
    for (int k = 0; k < p; ++k)
        for (SparseMatrix::InnerIterator it(pb.G, k); it; ++it) A.emplace_back(static_cast<int>(it.row()), n + k, -it.value());
    for (int col = 0; col < mc; ++col) {
        const auto [c, j] = comp[(size_t)col];
        const double s = componentScale(pb, c, j);
        for (SparseMatrix::InnerIterator it(pb.H, 3 * c + j); it; ++it) A.emplace_back(static_cast<int>(it.row()), n + p + col, -s * it.value());
        cp.q[n + p + col] = s * pb.w[3 * c + j] + (j == 0 ? shift[c] : 0.0);
        A.emplace_back(n + col, n + p + col, -1.0);
    }
    cp.b = Vector::Zero(n + mc);
    cp.b.head(n) = pb.f;
    cp.P = fromTriplets(ny, ny, P);
    cp.A = fromTriplets(n + mc, ny, A);
    cp.dims = coneDims(n, order);
    return cp;
}

Solution recoverSolution(const Problem& pb, Form form, const Vector& x, const Vector& z) {
    const int n = pb.n(), p = pb.p();
    const auto comp = coneComponents(contactOrder(pb));
    Solution sol;
    sol.v = x.head(n);
    sol.r = Vector::Zero(3 * pb.nc());
    if (form == Form::V) {
        sol.lambda = -z.head(p);
        for (std::size_t row = 0; row < comp.size(); ++row) {
            const auto [c, j] = comp[row];
            sol.r[3 * c + j] = componentScale(pb, c, j) * z[p + static_cast<Eigen::Index>(row)];
        }
    } else {
        sol.lambda = x.segment(n, p);
        for (std::size_t col = 0; col < comp.size(); ++col) {
            const auto [c, j] = comp[col];
            sol.r[3 * c + j] = componentScale(pb, c, j) * x[n + p + static_cast<Eigen::Index>(col)];
        }
    }
    return sol;
}

}  // namespace cfc
