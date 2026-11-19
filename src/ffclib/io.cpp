#include <ffclib/io.hpp>

#include <hdf5.h>

#include <algorithm>
#include <climits>
#include <cstdint>
#include <filesystem>
#include <stdexcept>

// fclib's layout (fclib.h): every value is a dataset, sizes are int datasets of one value, strings are
// fixed-length scalar datasets, and a matrix is a group {nzmax, m, n, nz, p, i, x} in CSparse's convention
// (nz >= 0: nz triplets with row i[k] and column p[k]; nz = -1: compressed columns; nz = -2: compressed
// rows). Arrays are written compressed; fclib's reader decompresses them transparently.
namespace ffclib {

namespace {

// gzip level of the arrays (after HDF5's shuffle filter), which halves the size of typical problem files.
constexpr unsigned kDeflateLevel = 4;
// Values per chunk of a compressed array; HDF5 limits a chunk to 4 GiB.
constexpr hsize_t kChunkValues = hsize_t(1) << 20;

[[noreturn]] void fail(const std::string& what) {
    throw std::runtime_error("ffclib file: " + what);
}

void require(herr_t status, const std::string& what) {
    if (status < 0) fail("HDF5 call failed: " + what);
}

// Owns an HDF5 identifier.
class Hid {
   public:
    Hid(hid_t id, herr_t (*close)(hid_t), const std::string& what) : m_id(id), m_close(close) {
        if (id < 0) fail("HDF5 call failed: " + what);
    }
    ~Hid() { m_close(m_id); }
    Hid(const Hid&) = delete;
    Hid& operator=(const Hid&) = delete;
    operator hid_t() const { return m_id; }

   private:
    hid_t m_id;
    herr_t (*m_close)(hid_t);
};

// HDF5 prints its error stack by default; here failures are exceptions instead.
class QuietErrors {
   public:
    QuietErrors() {
        H5Eget_auto2(H5E_DEFAULT, &m_func, &m_data);
        H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    }
    ~QuietErrors() { H5Eset_auto2(H5E_DEFAULT, m_func, m_data); }
    QuietErrors(const QuietErrors&) = delete;
    QuietErrors& operator=(const QuietErrors&) = delete;

   private:
    H5E_auto2_t m_func{nullptr};
    void* m_data{nullptr};
};

Hid createGroup(hid_t loc, const std::string& name) {
    return Hid(H5Gcreate2(loc, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose, name);
}
Hid openGroup(hid_t loc, const std::string& name) {
    return Hid(H5Gopen2(loc, name.c_str(), H5P_DEFAULT), H5Gclose, name);
}
// Whether the relative path `name` (components separated by '/') exists below `loc`.
bool exists(hid_t loc, const std::string& name) {
    for (std::size_t end = name.find('/');; end = name.find('/', end + 1)) {
        if (H5Lexists(loc, name.substr(0, end).c_str(), H5P_DEFAULT) <= 0) return false;
        if (end == std::string::npos) return true;
    }
}

// ---------- writing ----------

void writeArray(hid_t loc, const char* name, hid_t type, const void* data, std::size_t count) {
    const hsize_t dims[1] = {count};
    Hid space(H5Screate_simple(1, dims, nullptr), H5Sclose, name);
    Hid plist(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, name);
    if (count > 1 && H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0) {
        const hsize_t chunk[1] = {std::min<hsize_t>(count, kChunkValues)};
        require(H5Pset_chunk(plist, 1, chunk), name);
        require(H5Pset_shuffle(plist), name);
        require(H5Pset_deflate(plist, kDeflateLevel), name);
    }
    Hid ds(H5Dcreate2(loc, name, type, space, H5P_DEFAULT, plist, H5P_DEFAULT), H5Dclose, name);
    if (count > 0) require(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), name);
}
void writeVector(hid_t loc, const char* name, const Vector& x) {
    writeArray(loc, name, H5T_NATIVE_DOUBLE, x.data(), static_cast<std::size_t>(x.size()));
}
void writeInt(hid_t loc, const char* name, int value) {
    writeArray(loc, name, H5T_NATIVE_INT, &value, 1);
}
// Fixed-length, null-terminated scalar string, as fclib writes it (H5LTmake_dataset_string).
void writeString(hid_t loc, const char* name, const std::string& value) {
    Hid type(H5Tcopy(H5T_C_S1), H5Tclose, name);
    require(H5Tset_size(type, value.size() + 1), name);
    Hid space(H5Screate(H5S_SCALAR), H5Sclose, name);
    Hid ds(H5Dcreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose, name);
    require(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value.c_str()), name);
}
// Compressed columns (nz = -1).
void writeMatrix(hid_t loc, const char* name, const SparseMatrix& A) {
    SparseMatrix B = A;
    B.makeCompressed();
    Hid g = createGroup(loc, name);
    writeInt(g, "nzmax", static_cast<int>(B.nonZeros()));
    writeInt(g, "m", static_cast<int>(B.rows()));
    writeInt(g, "n", static_cast<int>(B.cols()));
    writeInt(g, "nz", -1);
    writeArray(g, "p", H5T_NATIVE_INT, B.outerIndexPtr(), static_cast<std::size_t>(B.cols() + 1));
    writeArray(g, "i", H5T_NATIVE_INT, B.innerIndexPtr(), static_cast<std::size_t>(B.nonZeros()));
    writeArray(g, "x", H5T_NATIVE_DOUBLE, B.valuePtr(), static_cast<std::size_t>(B.nonZeros()));
}
// fclib's solution: v, l (if nb > 0), u = H^T v + w and r.
void writeSolution(hid_t loc, const Problem& pb, const Solution& s) {
    const Vector u = pb.H.transpose() * s.v + pb.w;
    writeVector(loc, "v", s.v);
    if (pb.nb() > 0) writeVector(loc, "l", s.lambda);
    writeVector(loc, "u", u);
    writeVector(loc, "r", s.r);
}

bool matches(const Problem& pb, const Solution& s) {
    return s.v.size() == pb.nv() && s.r.size() == 3 * pb.nc() && s.lambda.size() == pb.nb();
}

// ---------- reading ----------

template <typename T>
std::vector<T> readArray(hid_t loc, const std::string& name, hid_t memtype) {
    Hid ds(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), H5Dclose, name);
    Hid space(H5Dget_space(ds), H5Sclose, name);
    const hssize_t count = H5Sget_simple_extent_npoints(space);
    if (count < 0) fail(name + ": no simple dataspace");
    std::vector<T> v(static_cast<std::size_t>(count));
    if (count > 0) require(H5Dread(ds, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()), name);
    return v;
}
Vector readVector(hid_t loc, const std::string& name) {
    const auto v = readArray<double>(loc, name, H5T_NATIVE_DOUBLE);
    return Eigen::Map<const Vector>(v.data(), static_cast<Eigen::Index>(v.size()));
}
std::vector<std::int64_t> readInts(hid_t loc, const std::string& name) {
    return readArray<std::int64_t>(loc, name, H5T_NATIVE_INT64);
}
std::int64_t readInt(hid_t loc, const std::string& name) {
    const auto v = readInts(loc, name);
    if (v.size() != 1) fail(name + " must hold one value");
    return v[0];
}
// Fixed-length (fclib) or variable-length (e.g. h5py) string.
std::string readString(hid_t loc, const std::string& name) {
    Hid ds(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), H5Dclose, name);
    Hid ftype(H5Dget_type(ds), H5Tclose, name);
    Hid mtype(H5Tcopy(H5T_C_S1), H5Tclose, name);
    if (H5Tis_variable_str(ftype) > 0) {
        require(H5Tset_size(mtype, H5T_VARIABLE), name);
        char* c = nullptr;
        require(H5Dread(ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &c), name);
        std::string s = c ? c : "";
        H5free_memory(c);
        return s;
    }
    const std::size_t size = H5Tget_size(ftype);
    std::string s(size + 1, '\0');
    require(H5Tset_size(mtype, size + 1), name);
    require(H5Dread(ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, s.data()), name);
    return s.substr(0, s.find('\0'));
}

// All three storage schemes of fclib; duplicate triplets are summed, as in CSparse.
SparseMatrix readMatrix(hid_t loc, const std::string& name) {
    Hid g = openGroup(loc, name);
    const std::int64_t m = readInt(g, "m"), n = readInt(g, "n"), nz = readInt(g, "nz");
    const auto p = readInts(g, "p");
    const auto i = readInts(g, "i");
    const auto x = readArray<double>(g, "x", H5T_NATIVE_DOUBLE);
    if (m < 0 || n < 0 || m > INT_MAX || n >= INT_MAX) fail(name + ": invalid size");
    std::vector<Eigen::Triplet<double>> t;
    auto add = [&](std::int64_t row, std::int64_t col, std::size_t k) {
        if (row < 0 || row >= m || col < 0 || col >= n) fail(name + ": index out of range");
        t.emplace_back(static_cast<int>(row), static_cast<int>(col), x[k]);
    };
    if (nz >= 0) {
        const auto count = static_cast<std::size_t>(nz);
        if (p.size() < count || i.size() < count || x.size() < count) fail(name + ": fewer triplets than nz");
        for (std::size_t k = 0; k < count; ++k) add(i[k], p[k], k);
    } else if (nz == -1 || nz == -2) {
        const auto outer = static_cast<std::size_t>(nz == -1 ? n : m);
        if (p.size() < outer + 1 || p[0] < 0 || p[outer] > static_cast<std::int64_t>(std::min(i.size(), x.size()))) fail(name + ": pointers out of range");
        for (std::size_t j = 0; j < outer; ++j) {
            if (p[j + 1] < p[j]) fail(name + ": decreasing pointers");
            for (auto k = static_cast<std::size_t>(p[j]); k < static_cast<std::size_t>(p[j + 1]); ++k)
                nz == -1 ? add(i[k], static_cast<std::int64_t>(j), k) : add(static_cast<std::int64_t>(j), i[k], k);
        }
    } else {
        fail(name + ": unknown storage nz = " + std::to_string(nz));
    }
    SparseMatrix A(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(n));
    A.setFromTriplets(t.begin(), t.end());
    A.makeCompressed();
    return A;
}

Solution readSolution(hid_t loc, const Problem& pb, const std::string& what) {
    Solution s;
    s.v = readVector(loc, "v");
    s.r = readVector(loc, "r");
    s.lambda = pb.nb() > 0 ? readVector(loc, "l") : Vector();
    if (!matches(pb, s)) fail(what + " does not match the problem");
    return s;
}

}  // namespace

Record read(const std::string& path) {
    const QuietErrors quiet;
    if (!std::filesystem::is_regular_file(path)) fail(path + " does not exist");
    Hid file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose, "open " + path);
    if (!exists(file, "fclib_global")) fail(path + " holds no global problem (/fclib_global)");
    Hid g = openGroup(file, "fclib_global");
    if (exists(g, "spacedim") && readInt(g, "spacedim") != 3) fail(path + ": only 3D problems (spacedim = 3) are supported");

    Record rec;
    Problem& pb = rec.problem;
    pb.M = readMatrix(g, "M");
    pb.H = readMatrix(g, "H");
    pb.G = exists(g, "G") ? readMatrix(g, "G") : SparseMatrix(pb.M.rows(), 0);
    pb.f = readVector(g, "vectors/f");
    pb.w = readVector(g, "vectors/w");
    pb.mu = readVector(g, "vectors/mu");
    pb.b = exists(g, "vectors/b") ? readVector(g, "vectors/b") : Vector();
    pb.C = exists(g, "C") ? readMatrix(g, "C") : SparseMatrix(pb.b.size(), pb.b.size());
    try {
        check(pb);
    } catch (const std::invalid_argument& e) {
        fail(path + ": " + e.what());
    }
    if (exists(g, "info/title")) rec.info.title = readString(g, "info/title");
    if (exists(g, "info/description")) rec.info.description = readString(g, "info/description");
    if (exists(g, "info/math_info")) rec.info.math_info = readString(g, "info/math_info");

    if (exists(file, "solution")) rec.solution = readSolution(openGroup(file, "solution"), pb, path + ": /solution");
    if (exists(file, "guesses")) {
        Hid guesses = openGroup(file, "guesses");
        const std::int64_t count = readInt(guesses, "number_of_guesses");
        for (std::int64_t k = 1; k <= count; ++k) rec.guesses.push_back(readSolution(openGroup(guesses, std::to_string(k)), pb, path + ": /guesses/" + std::to_string(k)));
    }
    return rec;
}

void write(const std::string& path, const Record& rec) {
    const Problem& pb = rec.problem;
    check(pb);
    if (pb.nc() == 0) throw std::invalid_argument("ffclib file: fclib needs at least one contact");
    if (pb.M.nonZeros() > INT_MAX || pb.H.nonZeros() > INT_MAX || pb.G.nonZeros() > INT_MAX || pb.H.cols() >= INT_MAX)
        throw std::invalid_argument("ffclib file: the problem is too large for fclib's 32-bit indices");
    bool sizes = !rec.solution || matches(pb, *rec.solution);
    for (const Solution& s : rec.guesses) sizes = sizes && matches(pb, s);
    if (!sizes) throw std::invalid_argument("ffclib file: a solution or guess does not match the problem");

    const QuietErrors quiet;
    Hid file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose, "create " + path);
    {
        Hid g = createGroup(file, "fclib_global");
        writeInt(g, "spacedim", 3);
        writeMatrix(g, "M", pb.M);
        writeMatrix(g, "H", pb.H);
        if (pb.nb() > 0) writeMatrix(g, "G", pb.G);
        if (pb.C.nonZeros() > 0) writeMatrix(g, "C", pb.C);
        Hid vectors = createGroup(g, "vectors");
        writeVector(vectors, "f", pb.f);
        writeVector(vectors, "w", pb.w);
        writeVector(vectors, "mu", pb.mu);
        if (pb.nb() > 0) writeVector(vectors, "b", pb.b);
        Hid info = createGroup(g, "info");
        writeString(info, "title", rec.info.title);
        writeString(info, "description", rec.info.description);
        writeString(info, "math_info", rec.info.math_info);
    }
    if (rec.solution) writeSolution(createGroup(file, "solution"), pb, *rec.solution);
    if (!rec.guesses.empty()) {
        Hid guesses = createGroup(file, "guesses");
        writeInt(guesses, "number_of_guesses", static_cast<int>(rec.guesses.size()));
        for (std::size_t k = 0; k < rec.guesses.size(); ++k) writeSolution(createGroup(guesses, std::to_string(k + 1)), pb, rec.guesses[k]);
    }
}

}  // namespace ffclib
