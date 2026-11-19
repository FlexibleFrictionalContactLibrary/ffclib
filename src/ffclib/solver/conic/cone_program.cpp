#include <ffclib/solver/conic/cone_program.hpp>

#include <stdexcept>

namespace ffclib {

namespace {

using Triplets = std::vector<Eigen::Triplet<double>>;

// The contacts in the order of their rows: those without friction, then those with friction. Contacts that
// act on no velocity are left out.
struct ContactRows {
    std::vector<int> frictionless;
    std::vector<int> frictional;
};

ContactRows contactRows(const Problem& pb) {
    ContactRows rows;
    for (int c = 0; c < pb.nc(); ++c)
        if (actsOnVelocity(pb, c)) (pb.mu[c] > 0 ? rows.frictional : rows.frictionless).push_back(c);
    return rows;
}

}  // namespace

bool actsOnVelocity(const Problem& pb, int c) {
    for (int j = 0; j < 3; ++j)
        for (SparseMatrix::InnerIterator it(pb.H, 3 * c + j); it; ++it)
            if (it.value() != 0) return true;
    return false;
}

ConeProgram coneProgram(const Problem& pb) {
    const int nv = pb.nv(), nb = pb.nb();
    const ContactRows contacts = contactRows(pb);

    // pi_k for every bilateral row with a nonzero column of C; -1 otherwise.
    std::vector<int> pi(static_cast<std::size_t>(nb), -1);
    int npi = 0;
    for (int k = 0; k < nb; ++k)
        for (SparseMatrix::InnerIterator it(pb.C, k); it; ++it)
            if (it.value() != 0) {
                pi[k] = npi++;
                break;
            }

    Triplets P, A;
    for (int j = 0; j < nv; ++j)
        for (SparseMatrix::InnerIterator it(pb.M, j); it; ++it)
            if (it.row() <= j) P.emplace_back(static_cast<int>(it.row()), j, it.value());
    for (int k = 0; k < nb; ++k)
        for (SparseMatrix::InnerIterator it(pb.C, k); it; ++it)
            if (it.row() <= k && it.value() != 0) P.emplace_back(nv + pi[it.row()], nv + pi[k], it.value());

    ConeProgram cp;
    const int rows = nb + static_cast<int>(contacts.frictionless.size() + 3 * contacts.frictional.size());
    cp.b.resize(rows);
    // Zero rows: row k is column k of G and the entries C_kl of the rows l with pi_l.
    for (int k = 0; k < nb; ++k)
        for (SparseMatrix::InnerIterator it(pb.G, k); it; ++it) A.emplace_back(k, static_cast<int>(it.row()), it.value());
    for (int l = 0; l < nb; ++l)
        for (SparseMatrix::InnerIterator it(pb.C, l); it; ++it)
            if (it.value() != 0) A.emplace_back(static_cast<int>(it.row()), nv + pi[l], it.value());
    cp.b.head(nb) = -pb.b;
    // Contact rows: -s H_cj^T v + s_row = s w_cj with the scale s of the component.
    int row = nb;
    auto contactRow = [&](int c, int j, double scale) {
        for (SparseMatrix::InnerIterator it(pb.H, 3 * c + j); it; ++it) A.emplace_back(row, static_cast<int>(it.row()), -scale * it.value());
        cp.b[row++] = scale * pb.w[3 * c + j];
    };
    for (int c : contacts.frictionless) contactRow(c, 0, 1);
    for (int c : contacts.frictional) {
        contactRow(c, 0, 1);
        contactRow(c, 1, pb.mu[c]);
        contactRow(c, 2, pb.mu[c]);
    }

    cp.P.resize(nv + npi, nv + npi);
    cp.P.setFromTriplets(P.begin(), P.end());
    cp.A.resize(rows, nv + npi);
    cp.A.setFromTriplets(A.begin(), A.end());
    cp.q = Vector::Zero(nv + npi);
    cp.q.head(nv) = -pb.f;
    cp.dims.zero = nb;
    cp.dims.nonneg = static_cast<int>(contacts.frictionless.size());
    cp.dims.soc.assign(contacts.frictional.size(), 3);
    return cp;
}

Solution solutionFromDuals(const Problem& pb, const Vector& z) {
    const int nb = pb.nb();
    const ContactRows contacts = contactRows(pb);
    if (z.size() != nb + static_cast<Eigen::Index>(contacts.frictionless.size() + 3 * contacts.frictional.size())) throw std::invalid_argument("ffclib: the duals do not match the cone program");
    const Vector lambda = -z.head(nb);
    Vector r = Vector::Zero(3 * pb.nc());
    int row = nb;
    for (int c : contacts.frictionless) r[3 * c] = z[row++];
    for (int c : contacts.frictional) {
        r[3 * c] = z[row];
        r[3 * c + 1] = pb.mu[c] * z[row + 1];
        r[3 * c + 2] = pb.mu[c] * z[row + 2];
        row += 3;
    }
    return {velocity(pb, r, lambda), r, lambda};
}

}  // namespace ffclib
