// Tests of the cfc library:
//   1. round trip write -> read, bit-exact, with and without compliance and metadata;
//   2. interoperability with fclib: fclib reads our files, we read files that fclib wrote (no extension);
//   3. residuals against closed-form solutions (a sticking and a sliding point mass, a compliant spring);
//   4. the cone programs: their KKT and feasibility residuals at (v, r, lambda) are the problem's.
// Files of the FCLIB collection given as arguments are read and checked for consistency.
//
//   test_cfc [FCLIB file ...]

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>

#include <unistd.h>

#include <cfc/cone.hpp>
#include <cfc/cone_program.hpp>
#include <cfc/file.hpp>
#include <cfc/residuals.hpp>

extern "C" {
#include <fclib.h>
}

namespace {

int failures = 0;

void expect(bool ok, const std::string& what) {
    if (!ok) {
        ++failures;
        std::cerr << "FAILED: " << what << '\n';
    }
}

bool same(const cfc::SparseMatrix& A, const cfc::SparseMatrix& B) {
    return A.rows() == B.rows() && A.cols() == B.cols() && (cfc::SparseMatrix(A - B)).norm() == 0;
}
bool same(const cfc::Vector& a, const cfc::Vector& b) {
    return a.size() == b.size() && (a.size() == 0 || (a - b).lpNorm<Eigen::Infinity>() == 0);
}

cfc::SparseMatrix fromFclib(const fclib_matrix& m) {
    std::vector<Eigen::Triplet<double>> t;
    for (int j = 0; j < m.n; ++j)
        for (int k = m.p[j]; k < m.p[j + 1]; ++k) t.emplace_back(m.i[k], j, m.x[k]);
    cfc::SparseMatrix A(m.m, m.n);
    A.setFromTriplets(t.begin(), t.end());
    return A;
}

cfc::Vector vec(const double* x, int n) {
    return x ? cfc::Vector(Eigen::Map<const cfc::Vector>(x, n)) : cfc::Vector::Zero(n);
}

// A random problem: nb rigid bodies, nc contacts (every third frictionless), p bilateral rows, of which
// the first half compliant (with a coupling between rows 0 and 1), and a random solution.
cfc::Record randomRecord(int nb, int nc, int p, bool compliant, bool diagonalMass, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> u(-1, 1), pos(0.5, 2);
    const int n = 6 * nb;
    cfc::Record rec;
    cfc::Problem& pb = rec.problem;
    std::vector<Eigen::Triplet<double>> t;
    for (int i = 0; i < n; ++i) t.emplace_back(i, i, pos(gen));
    if (!diagonalMass)
        for (int b = 0; b < nb; ++b) {
            const double c = 0.1 * u(gen);
            t.emplace_back(6 * b + 3, 6 * b + 4, c);
            t.emplace_back(6 * b + 4, 6 * b + 3, c);
        }
    pb.M.resize(n, n);
    pb.M.setFromTriplets(t.begin(), t.end());
    t.clear();
    for (int c = 0; c < nc; ++c) {
        const int b = c % nb;
        for (int j = 0; j < 3; ++j) {
            if (c % 3 == 2 && j > 0) continue;  // frictionless: zero tangential columns
            for (int i = 0; i < 6; ++i) t.emplace_back(6 * b + i, 3 * c + j, u(gen));
        }
    }
    pb.H.resize(n, 3 * nc);
    pb.H.setFromTriplets(t.begin(), t.end());
    t.clear();
    for (int k = 0; k < p; ++k)
        for (int i = 0; i < 6; ++i) t.emplace_back(6 * (k % nb) + i, k, u(gen));
    pb.G.resize(n, p);
    pb.G.setFromTriplets(t.begin(), t.end());
    t.clear();
    if (compliant) {
        for (int k = 0; k < p / 2; ++k) t.emplace_back(k, k, 0.5 + pos(gen));
        if (p >= 2) {
            t.emplace_back(0, 1, 0.2);
            t.emplace_back(1, 0, 0.2);
        }
    }
    pb.C.resize(p, p);
    pb.C.setFromTriplets(t.begin(), t.end());
    pb.f = cfc::Vector::NullaryExpr(n, [&] { return u(gen); });
    pb.w = cfc::Vector::NullaryExpr(3 * nc, [&] { return u(gen); });
    pb.b = cfc::Vector::NullaryExpr(p, [&] { return u(gen); });
    pb.mu = cfc::Vector::NullaryExpr(nc, [&] { return 0.3 + 0.2 * u(gen); });
    for (int c = 2; c < nc; c += 3) {
        pb.mu[c] = 0;
        pb.w.segment<2>(3 * c + 1).setZero();
    }
    pb.shift = cfc::Vector::NullaryExpr(nc, [&] { return std::abs(u(gen)); });
    for (cfc::SparseMatrix* A : {&pb.M, &pb.H, &pb.G, &pb.C}) A->makeCompressed();

    cfc::Solution s;
    s.v = cfc::Vector::NullaryExpr(n, [&] { return u(gen); });
    s.r = cfc::Vector::NullaryExpr(3 * nc, [&] { return u(gen); });
    s.lambda = cfc::Vector::NullaryExpr(p, [&] { return u(gen); });
    for (int c = 2; c < nc; c += 3) s.r.segment<2>(3 * c + 1).setZero();
    rec.reference = s;
    s.v *= 0.5;
    rec.guess = s;
    return rec;
}

void fillMetadata(cfc::Record& rec) {
    cfc::Metadata& md = rec.meta;
    const int nc = rec.problem.nc(), p = rec.problem.p(), nb = rec.problem.n() / 6;
    md.present = true;
    md.producer = "test_cfc";
    md.producer_version = "v0";
    md.scene = "random";
    md.step = 7;
    md.time = 0.25;
    md.h = 1e-3;
    md.theta = 0.6;
    md.units = "SI";
    for (int c = 0; c < nc; ++c) {
        md.contact_body_a.push_back(c % nb);
        md.contact_body_b.push_back(-1);
        md.contact_point.emplace_back(c, 2 * c, 3 * c);
        md.contact_frame.push_back(cfc::Matrix3::Identity() * (c + 1));
    }
    md.contact_gap = cfc::Vector::LinSpaced(nc, -1e-3, 0);
    md.restitution_normal = cfc::Vector::Constant(nc, 0.5);
    md.restitution_tangential = cfc::Vector::Zero(nc);
    for (int b = 0; b < nb; ++b) {
        md.body_dof_offset.push_back(6 * b);
        md.body_dof_count.push_back(6);
        md.body_type.push_back(0);
    }
    for (int k = 0; k < p; ++k) {
        md.bilateral_kind.push_back(k % 2);
        md.bilateral_constraint.push_back(k / 3);
        md.bilateral_constraint_row.push_back(k % 3);
        md.bilateral_body_a.push_back(k % nb);
        md.bilateral_body_b.push_back(-1);
    }
    md.bilateral_law_compliance = cfc::Vector::LinSpaced(p, 0, 1);
    md.bilateral_law_scale = cfc::Vector::Constant(p, 4);
    md.solver = "solver";
    md.law = cfc::Law::Convexified;
    md.iterations = 42;
    md.solve_seconds = 0.125;
    md.producer_natural_map = 1e-9;
    md.producer_natural_map_coulomb = 2e-9;
    md.producer_bilateral = 3e-9;
    md.previous_file = "prev.cfc.h5";
}

void expectSameProblem(const cfc::Problem& a, const cfc::Problem& b, const std::string& what) {
    expect(same(a.M, b.M) && same(a.H, b.H) && same(a.G, b.G) && same(a.C, b.C), what + ": matrices");
    expect(same(a.f, b.f) && same(a.w, b.w) && same(a.b, b.b) && same(a.mu, b.mu) && same(a.shift, b.shift), what + ": vectors");
}
void expectSameSolution(const cfc::Solution& a, const cfc::Solution& b, const std::string& what) {
    expect(same(a.v, b.v) && same(a.r, b.r) && same(a.lambda, b.lambda), what);
}

void testRoundTrip(const std::string& dir) {
    for (const bool compliant : {false, true})
        for (const bool meta : {false, true}) {
            cfc::Record rec = randomRecord(3, 8, 4, compliant, compliant, 11);
            if (meta) fillMetadata(rec);
            const std::string what = std::string("round trip") + (compliant ? " compliant" : "") + (meta ? " metadata" : "");
            const std::string path = dir + "/roundtrip.cfc.h5";
            cfc::write(path, rec);
            const cfc::Record back = cfc::read(path);
            expectSameProblem(rec.problem, back.problem, what);
            expect(back.reference && back.guess, what + ": solutions present");
            if (back.reference) expectSameSolution(*rec.reference, *back.reference, what + ": reference");
            if (back.guess) expectSameSolution(*rec.guess, *back.guess, what + ": guess");
            expect(back.meta.present == meta && back.mass_asymmetry == 0, what + ": flags");
            if (meta) {
                const cfc::Metadata &a = rec.meta, &b = back.meta;
                expect(a.producer == b.producer && a.producer_version == b.producer_version && a.scene == b.scene && a.step == b.step && a.time == b.time && a.h == b.h && a.theta == b.theta &&
                           a.units == b.units,
                       what + ": provenance");
                bool frames = a.contact_frame.size() == b.contact_frame.size() && a.contact_point.size() == b.contact_point.size();
                for (std::size_t c = 0; frames && c < a.contact_frame.size(); ++c) frames = a.contact_frame[c] == b.contact_frame[c] && a.contact_point[c] == b.contact_point[c];
                expect(frames && a.contact_body_a == b.contact_body_a && a.contact_body_b == b.contact_body_b && same(a.contact_gap, b.contact_gap) &&
                           same(a.restitution_normal, b.restitution_normal) && same(a.restitution_tangential, b.restitution_tangential),
                       what + ": contacts");
                expect(a.body_dof_offset == b.body_dof_offset && a.body_dof_count == b.body_dof_count && a.body_type == b.body_type, what + ": bodies");
                expect(a.bilateral_kind == b.bilateral_kind && a.bilateral_constraint == b.bilateral_constraint && a.bilateral_constraint_row == b.bilateral_constraint_row &&
                           a.bilateral_body_a == b.bilateral_body_a && a.bilateral_body_b == b.bilateral_body_b && same(a.bilateral_law_compliance, b.bilateral_law_compliance) &&
                           same(a.bilateral_law_scale, b.bilateral_law_scale),
                       what + ": bilateral rows");
                expect(a.solver == b.solver && a.law == b.law && a.iterations == b.iterations && a.solve_seconds == b.solve_seconds && a.producer_natural_map == b.producer_natural_map &&
                           a.producer_natural_map_coulomb == b.producer_natural_map_coulomb && a.producer_bilateral == b.producer_bilateral && a.previous_file == b.previous_file,
                       what + ": reference");
            }
        }
}

void testFclib(const std::string& dir) {
    // fclib reads our file: the same problem without C and the shift, the same solutions.
    cfc::Record rec = randomRecord(4, 9, 5, true, true, 5);
    fillMetadata(rec);
    const std::string ours = dir + "/ours.cfc.h5";
    cfc::write(ours, rec);
    const cfc::Problem& pb = rec.problem;
    fclib_global* g = fclib_read_global(ours.c_str());
    expect(g != nullptr, "fclib reads our file");
    if (g) {
        expect(g->spacedim == 3 && same(fromFclib(*g->M), pb.M) && same(fromFclib(*g->H), pb.H) && g->G && same(fromFclib(*g->G), pb.G), "fclib: matrices");
        expect(same(vec(g->f, pb.n()), pb.f) && same(vec(g->w, 3 * pb.nc()), pb.w) && same(vec(g->mu, pb.nc()), pb.mu) && same(vec(g->b, pb.p()), pb.b), "fclib: vectors");
        fclib_delete_global(g);
        std::free(g);
    }
    if (fclib_solution* s = fclib_read_solution(ours.c_str())) {
        const cfc::Vector u = pb.H.transpose() * rec.reference->v + pb.w;
        expect(same(vec(s->v, pb.n()), rec.reference->v) && same(vec(s->r, 3 * pb.nc()), rec.reference->r) && same(vec(s->l, pb.p()), rec.reference->lambda) && same(vec(s->u, 3 * pb.nc()), u),
               "fclib: solution");
        fclib_delete_solutions(s, 1);
    } else {
        expect(false, "fclib reads our solution");
    }

    // We read a file that fclib wrote: no extension, so C = 0 and shift = 0.
    const std::string theirs = dir + "/fclib.hdf5";
    std::filesystem::remove(theirs);
    auto toFclib = [](const cfc::SparseMatrix& A, std::vector<int>& p, std::vector<int>& i, std::vector<double>& x) {
        p.assign(A.outerIndexPtr(), A.outerIndexPtr() + A.cols() + 1);
        i.assign(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros());
        x.assign(A.valuePtr(), A.valuePtr() + A.nonZeros());
        fclib_matrix m{};
        m.nzmax = static_cast<int>(A.nonZeros());
        m.m = static_cast<int>(A.rows());
        m.n = static_cast<int>(A.cols());
        m.p = p.data();
        m.i = i.data();
        m.x = x.data();
        m.nz = -1;
        return m;
    };
    std::vector<int> Mp, Mi, Hp, Hi, Gp, Gi;
    std::vector<double> Mx, Hx, Gx;
    fclib_matrix M = toFclib(pb.M, Mp, Mi, Mx), H = toFclib(pb.H, Hp, Hi, Hx), G = toFclib(pb.G, Gp, Gi, Gx);
    cfc::Vector f = pb.f, w = pb.w, b = pb.b, mu = pb.mu;
    fclib_global problem{};
    problem.M = &M;
    problem.H = &H;
    problem.G = &G;
    problem.mu = mu.data();
    problem.f = f.data();
    problem.b = b.data();
    problem.w = w.data();
    problem.spacedim = 3;
    expect(fclib_write_global(&problem, theirs.c_str()) == 1, "fclib writes");
    cfc::Vector v = rec.reference->v, r = rec.reference->r, l = rec.reference->lambda, u = pb.H.transpose() * v + pb.w;
    fclib_solution sol{v.data(), u.data(), r.data(), l.data()};
    expect(fclib_write_solution(&sol, theirs.c_str()) == 1, "fclib writes a solution");
    const cfc::Record back = cfc::read(theirs);
    cfc::Problem rigid = pb;
    rigid.C = cfc::SparseMatrix(pb.p(), pb.p());
    rigid.shift = cfc::Vector::Zero(pb.nc());
    expectSameProblem(rigid, back.problem, "reading fclib's file");
    expect(back.reference && !back.guess && !back.meta.present, "reading fclib's file: solution, no guess, no metadata");
    if (back.reference) expectSameSolution(*rec.reference, *back.reference, "reading fclib's file: solution");
}

// One point mass (m = 2) on the ground, contact frame (N, T1, T2) = (e_z, e_x, e_y), mu = 0.5, and a
// slider (dof 3, m = 2) on a compliant spring row G = e_3, C = 0.5, b = -0.2: lambda = 0.2, v_3 = 0.1.
// f is the momentum without forces. Sticking: f = (0.3, 0, -1) -> v = 0, r = (1, -0.3, 0). Sliding:
// f = (2, 0, -1) -> r = (1, -0.5, 0), v = (0.75, 0, 0), De Saxce shift mu ||u_T|| = 0.375.
void testResiduals() {
    cfc::Problem pb;
    pb.M.resize(4, 4);
    for (int i = 0; i < 4; ++i) pb.M.insert(i, i) = 2;
    pb.H.resize(4, 3);
    pb.H.insert(2, 0) = 1;
    pb.H.insert(0, 1) = 1;
    pb.H.insert(1, 2) = 1;
    pb.G.resize(4, 1);
    pb.G.insert(3, 0) = 1;
    pb.C.resize(1, 1);
    pb.C.insert(0, 0) = 0.5;
    pb.w = cfc::Vector::Zero(3);
    pb.mu = cfc::Vector::Constant(1, 0.5);
    pb.b = cfc::Vector::Constant(1, -0.2);
    const double tol = 1e-15;

    pb.f = Eigen::Vector4d(0.3, 0, -1, 0);
    pb.shift = cfc::Vector::Zero(1);
    const cfc::Vector lambda = cfc::Vector::Constant(1, 0.2);
    const cfc::Solution stick{Eigen::Vector4d(0, 0, 0, 0.1), cfc::Vector3(1, -0.3, 0), lambda};
    cfc::Residuals r = cfc::residuals(pb, stick, pb.shift);
    expect(r.balance < tol && r.bilateral < tol && r.natural_map < tol && r.natural_map_coulomb < tol && r.cone_primal < tol && r.cone_dual < tol && r.complementarity < tol,
           "residuals: sticking contact solved exactly");
    expect(std::abs(r.velocity_scale - std::sqrt(0.15 * 0.15 + 0.5 * 0.5)) < tol, "residuals: velocity scale");

    pb.f = Eigen::Vector4d(2, 0, -1, 0);
    const cfc::Solution slide{Eigen::Vector4d(0.75, 0, 0, 0.1), cfc::Vector3(1, -0.5, 0), lambda};
    r = cfc::residuals(pb, slide, cfc::Vector::Constant(1, 0.375));
    expect(r.balance < tol && r.bilateral < tol && r.natural_map < tol && r.natural_map_coulomb < tol && r.complementarity < tol, "residuals: sliding contact solved exactly");
    r = cfc::residuals(pb, slide, cfc::Vector::Zero(1));
    expect(r.natural_map > 0.1 && r.natural_map_coulomb < tol, "residuals: the unshifted relaxation does not admit sliding");

    cfc::Solution wrong = slide;
    wrong.v[0] += 1e-3;
    wrong.lambda[0] += 1e-3;
    r = cfc::residuals(pb, wrong, cfc::Vector::Constant(1, 0.375));
    expect(std::abs(r.balance - 1e-3) < 1e-12 && std::abs(r.bilateral - 0.5e-3) < 1e-12 && r.natural_map > 0, "residuals: perturbed solution");
}

// The cone programs evaluated at a candidate (v, r, lambda) reproduce the problem's equations: for form V
// with x = (v, lambda_K), z = (-lambda, S^-1 r), the stationarity Px + q + A'z = (M v - f - H r - G lambda, 0)
// and the slacks b - Ax = (G^T v + C lambda + b, S u_hat); for form F with y = (v, lambda, S^-1 r) the
// zero rows A y - b_cp = M v - G lambda - H r - f.
void testConePrograms() {
    for (const bool diagonal : {true, false}) {
        const cfc::Record rec = randomRecord(3, 7, 6, true, diagonal, 3);
        const cfc::Problem& pb = rec.problem;
        const cfc::Solution& s = *rec.reference;
        const int n = pb.n(), p = pb.p();
        const cfc::ContactOrder order = cfc::contactOrder(pb);
        std::vector<std::pair<int, int>> comp;
        for (int c : order.frictionless) comp.emplace_back(c, 0);
        for (int c : order.frictional)
            for (int j = 0; j < 3; ++j) comp.emplace_back(c, j);
        const int mc = static_cast<int>(comp.size());
        auto scale = [&](int c, int j) { return j == 0 ? 1.0 : pb.mu[c]; };

        cfc::Vector force = pb.M * s.v - pb.f - pb.H * s.r - pb.G * s.lambda;
        cfc::Vector u = pb.H.transpose() * s.v + pb.w;
        const double tol = 1e-12;

        const cfc::ConeProgram V = cfc::buildConeProgram(pb, pb.shift, cfc::Form::V);
        const int na = static_cast<int>(V.P.cols()) - n;
        expect(V.dims.zero == p && V.A.rows() == p + mc && na == p / 2, "form V: sizes");
        cfc::Vector x(n + na), z(p + mc);
        x.head(n) = s.v;
        x.tail(na) = s.lambda.head(na);  // the compliant rows are the first p / 2
        z.head(p) = -s.lambda;
        for (int row = 0; row < mc; ++row) z[p + row] = s.r[3 * comp[(size_t)row].first + comp[(size_t)row].second] / scale(comp[(size_t)row].first, comp[(size_t)row].second);
        const cfc::SparseMatrix Pfull = V.P.selfadjointView<Eigen::Upper>();
        const cfc::Vector stat = Pfull * x + V.q + V.A.transpose() * z;
        expect((stat.head(n) - force).lpNorm<Eigen::Infinity>() < tol && stat.tail(na).lpNorm<Eigen::Infinity>() < tol, "form V: stationarity");
        const cfc::Vector slack = V.b - V.A * x;
        const cfc::Vector bil = pb.G.transpose() * s.v + pb.C * s.lambda + pb.b;
        expect((slack.head(p) + bil).lpNorm<Eigen::Infinity>() < tol, "form V: bilateral rows");
        double err = 0;
        for (int row = 0; row < mc; ++row) {
            const auto [c, j] = comp[(size_t)row];
            err = std::max(err, std::abs(slack[p + row] - scale(c, j) * (u[3 * c + j] + (j == 0 ? pb.shift[c] : 0.0))));
        }
        expect(err < tol, "form V: contact rows are S u_hat");
        const cfc::Solution sv = cfc::recoverSolution(pb, cfc::Form::V, x, z);
        expect((sv.v - s.v).norm() == 0 && (sv.lambda - s.lambda).norm() == 0 && (sv.r - s.r).lpNorm<Eigen::Infinity>() < tol, "form V: recovered solution");

        const cfc::ConeProgram F = cfc::buildConeProgram(pb, pb.shift, cfc::Form::F);
        expect(F.dims.zero == n && F.A.rows() == n + mc && F.P.cols() == n + p + mc, "form F: sizes");
        cfc::Vector y(n + p + mc);
        y.head(n) = s.v;
        y.segment(n, p) = s.lambda;
        for (int col = 0; col < mc; ++col) y[n + p + col] = z[p + col];
        const cfc::Vector rowsF = F.A * y - F.b;
        expect((rowsF.head(n) - force).lpNorm<Eigen::Infinity>() < tol, "form F: balance rows");
        const cfc::Solution sf = cfc::recoverSolution(pb, cfc::Form::F, y, cfc::Vector::Zero(F.A.rows()));
        expect((sf.v - s.v).norm() == 0 && (sf.lambda - s.lambda).norm() == 0 && (sf.r - s.r).lpNorm<Eigen::Infinity>() < tol, "form F: recovered solution");
        // The bilateral block of form F: C in P, b in q.
        const cfc::SparseMatrix PF = F.P.selfadjointView<Eigen::Upper>();
        expect((cfc::SparseMatrix(PF.block(n, n, p, p)) - pb.C).norm() == 0 && (F.q.segment(n, p) - pb.b).norm() == 0, "form F: compliance block");
    }
}

void testCollection(const std::string& path) {
    try {
        const cfc::Record rec = cfc::read(path);
        const cfc::Problem& pb = rec.problem;
        expect(!rec.meta.present && !pb.hasCompliance() && pb.shift.isZero(0), path + ": plain fclib file");
        if (rec.reference) {
            const cfc::Residuals r = cfc::residuals(pb, *rec.reference, pb.shift);
            std::cout << path << ": n " << pb.n() << ", nc " << pb.nc() << ", p " << pb.p() << ", mass asymmetry " << rec.mass_asymmetry << ", reference natural map (Coulomb) "
                      << r.natural_map_coulomb << '\n';
        }
    } catch (const std::exception& e) {
        expect(false, path + ": " + e.what());
    }
}

}  // namespace

int main(int argc, char** argv) {
    const std::string dir = (std::filesystem::temp_directory_path() / ("test_cfc_" + std::to_string(::getpid()))).string();
    std::filesystem::create_directories(dir);
    testRoundTrip(dir);
    testFclib(dir);
    testResiduals();
    testConePrograms();
    for (int i = 1; i < argc; ++i) testCollection(argv[i]);
    std::filesystem::remove_all(dir);

    // Reading a file that is not there fails with an exception, not an exit.
    bool threw = false;
    try {
        cfc::read(dir + "/missing.h5");
    } catch (const std::runtime_error&) {
        threw = true;
    }
    expect(threw, "missing file throws");

    std::cout << (failures ? "test_cfc: FAILED " : "test_cfc: all passed") << (failures ? std::to_string(failures) : "") << '\n';
    return failures ? 1 : 0;
}
