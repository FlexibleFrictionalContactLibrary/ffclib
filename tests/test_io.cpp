// Reading and writing files:
//   1. round trip write -> read, bit-exact, with and without compliance;
//   2. fclib reads our files, and we read files that fclib wrote;
//   3. check();
//   4. invalid files fail with an exception.

#include <cstdlib>
#include <filesystem>
#include <stdexcept>

#include <unistd.h>

#include <hdf5_hl.h>

#include "common.hpp"

extern "C" {
#include <fclib.h>
}

namespace {

using test::expect;

bool same(const ffclib::SparseMatrix& A, const ffclib::SparseMatrix& B) {
    return A.rows() == B.rows() && A.cols() == B.cols() && ffclib::SparseMatrix(A - B).norm() == 0;
}
bool same(const ffclib::Vector& a, const ffclib::Vector& b) {
    return a.size() == b.size() && (a - b).cwiseAbs().sum() == 0;
}
bool same(const ffclib::Problem& a, const ffclib::Problem& b) {
    return same(a.M, b.M) && same(a.H, b.H) && same(a.G, b.G) && same(a.C, b.C) && same(a.f, b.f) && same(a.w, b.w) && same(a.b, b.b) && same(a.mu, b.mu);
}
bool same(const ffclib::Solution& a, const ffclib::Solution& b) {
    return same(a.v, b.v) && same(a.r, b.r) && same(a.lambda, b.lambda);
}

ffclib::SparseMatrix fromFclib(const fclib_matrix& m) {
    std::vector<Eigen::Triplet<double>> t;
    for (int j = 0; j < m.n; ++j)
        for (int k = m.p[j]; k < m.p[j + 1]; ++k) t.emplace_back(m.i[k], j, m.x[k]);
    ffclib::SparseMatrix A(m.m, m.n);
    A.setFromTriplets(t.begin(), t.end());
    return A;
}
ffclib::Vector fromFclib(const double* x, int n) {
    return Eigen::Map<const ffclib::Vector>(x, n);
}

void testRoundTrip(const std::string& dir) {
    for (const bool compliant : {false, true}) {
        const ffclib::Record rec = test::randomRecord(3, 8, 4, compliant, compliant, 11);
        const std::string what = compliant ? "round trip, compliant" : "round trip";
        const std::string path = dir + "/roundtrip.hdf5";
        ffclib::write(path, rec);
        const ffclib::Record back = ffclib::read(path);
        expect(same(rec.problem, back.problem), what + ": problem");
        expect(back.solution && same(*rec.solution, *back.solution), what + ": solution");
        expect(back.guesses.size() == 2 && same(rec.guesses[0], back.guesses[0]) && same(rec.guesses[1], back.guesses[1]), what + ": guesses");
        expect(back.info.title == rec.info.title && back.info.description == rec.info.description && back.info.math_info == rec.info.math_info, what + ": info");
    }
}

void testFclib(const std::string& dir) {
    // fclib reads our file: the same problem without C, the same solution and guesses.
    const ffclib::Record rec = test::randomRecord(4, 9, 6, true, true, 5);
    const ffclib::Problem& pb = rec.problem;
    const std::string ours = dir + "/ours.hdf5";
    ffclib::write(ours, rec);
    if (fclib_global* g = fclib_read_global(ours.c_str())) {
        expect(g->spacedim == 3 && same(fromFclib(*g->M), pb.M) && same(fromFclib(*g->H), pb.H) && g->G && same(fromFclib(*g->G), pb.G), "fclib reads our matrices");
        expect(same(fromFclib(g->f, pb.nv()), pb.f) && same(fromFclib(g->w, 3 * pb.nc()), pb.w) && same(fromFclib(g->mu, pb.nc()), pb.mu) && same(fromFclib(g->b, pb.nb()), pb.b),
               "fclib reads our vectors");
        fclib_delete_global(g);
        std::free(g);
    } else {
        expect(false, "fclib reads our problem");
    }
    if (fclib_solution* s = fclib_read_solution(ours.c_str())) {
        const ffclib::Vector u = pb.H.transpose() * rec.solution->v + pb.w;
        expect(same(fromFclib(s->v, pb.nv()), rec.solution->v) && same(fromFclib(s->r, 3 * pb.nc()), rec.solution->r) && same(fromFclib(s->l, pb.nb()), rec.solution->lambda) &&
                   same(fromFclib(s->u, 3 * pb.nc()), u),
               "fclib reads our solution");
        fclib_delete_solutions(s, 1);
    } else {
        expect(false, "fclib reads our solution");
    }
    int count = 0;
    if (fclib_solution* s = fclib_read_guesses(ours.c_str(), &count)) {
        expect(count == 2 && same(fromFclib(s[1].r, 3 * pb.nc()), rec.guesses[1].r), "fclib reads our guesses");
        fclib_delete_solutions(s, count);
    } else {
        expect(false, "fclib reads our guesses");
    }

    // We read a file that fclib wrote: C = 0.
    const std::string theirs = dir + "/fclib.hdf5";
    std::filesystem::remove(theirs);
    struct Columns {
        std::vector<int> p, i;
        std::vector<double> x;
        fclib_matrix m{};
        explicit Columns(const ffclib::SparseMatrix& A)
            : p(A.outerIndexPtr(), A.outerIndexPtr() + A.cols() + 1), i(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros()), x(A.valuePtr(), A.valuePtr() + A.nonZeros()) {
            m = {static_cast<int>(A.nonZeros()), static_cast<int>(A.rows()), static_cast<int>(A.cols()), p.data(), i.data(), x.data(), -1, nullptr};
        }
    } M(pb.M), H(pb.H), G(pb.G);
    ffclib::Vector f = pb.f, w = pb.w, b = pb.b, mu = pb.mu;
    fclib_global problem{};
    problem.M = &M.m;
    problem.H = &H.m;
    problem.G = &G.m;
    problem.mu = mu.data();
    problem.f = f.data();
    problem.b = b.data();
    problem.w = w.data();
    problem.spacedim = 3;
    expect(fclib_write_global(&problem, theirs.c_str()) == 1, "fclib writes");
    ffclib::Vector v = rec.solution->v, r = rec.solution->r, l = rec.solution->lambda, u = pb.H.transpose() * v + pb.w;
    fclib_solution sol{v.data(), u.data(), r.data(), l.data()};
    expect(fclib_write_solution(&sol, theirs.c_str()) == 1, "fclib writes a solution");
    const ffclib::Record back = ffclib::read(theirs);
    ffclib::Problem rigid = pb;
    rigid.C = ffclib::SparseMatrix(pb.nb(), pb.nb());
    expect(same(rigid, back.problem), "we read fclib's problem");
    expect(back.solution && same(*rec.solution, *back.solution) && back.guesses.empty(), "we read fclib's solution");
}

bool rejected(const ffclib::Problem& pb) {
    try {
        ffclib::check(pb);
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

void testCheck() {
    ffclib::Problem pb = test::randomRecord(2, 3, 4, true, false, 7).problem;
    expect(!rejected(pb), "check accepts a valid problem");
    pb.C.coeffRef(1, 3) = 0.1;
    expect(rejected(pb), "check rejects an asymmetric C");
    pb.C.coeffRef(1, 3) = 0;
    pb.C.prune(0.0);
    pb.mu[1] = -0.1;
    expect(rejected(pb), "check rejects a negative friction coefficient");
    pb.mu[1] = 0.1;
    pb.w.conservativeResize(5);
    expect(rejected(pb), "check rejects a wrong size");
}

template <typename F>
bool throwsRuntimeError(F f) {
    try {
        f();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

void testInvalidFiles(const std::string& dir) {
    expect(throwsRuntimeError([&] { ffclib::read(dir + "/missing.hdf5"); }), "a missing file throws");

    // Decreasing column pointers of M must not be followed.
    const std::string path = dir + "/corrupt.hdf5";
    ffclib::write(path, test::randomRecord(2, 3, 2, false, false, 3));
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    const hid_t M = H5Gopen2(file, "/fclib_global/M", H5P_DEFAULT);
    H5Ldelete(M, "p", H5P_DEFAULT);
    std::vector<int> p(13, 0);
    p[1] = 5;  // then 0 again
    const hsize_t size = p.size();
    H5LTmake_dataset_int(M, "p", 1, &size, p.data());
    H5Gclose(M);
    H5Fclose(file);
    expect(throwsRuntimeError([&] { ffclib::read(path); }), "corrupt column pointers throw");
}

}  // namespace

int main() {
    const std::string dir = (std::filesystem::temp_directory_path() / ("ffclib_test_io_" + std::to_string(::getpid()))).string();
    std::filesystem::create_directories(dir);
    testRoundTrip(dir);
    testFclib(dir);
    testCheck();
    testInvalidFiles(dir);
    std::filesystem::remove_all(dir);
    return test::finish("test_io");
}
