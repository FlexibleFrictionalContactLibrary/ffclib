#include <cfc/file.hpp>

#include <hdf5.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <vector>

// The fclib part is written in fclib's layout (fclib.h of FrictionalContactLibrary/fclib): every
// value a dataset, integers as native int, strings as fixed-length scalar datasets, matrices as
// groups {nzmax, m, n, nz, p, i, x} in CSparse's convention (nz >= 0: nz triplets with row i and
// column p; nz = -1: compressed columns; nz = -2: compressed rows). Large datasets are chunked and
// compressed, which fclib's reader handles transparently.
namespace cfc {

namespace {

[[noreturn]] void fail(const std::string& what) {
    throw std::runtime_error("cfc file: " + what);
}

void check(herr_t status, const std::string& what) {
    if (status < 0) fail("HDF5 call failed: " + what);
}

// Owns an HDF5 identifier.
class Hid {
   public:
    Hid(hid_t id, herr_t (*close)(hid_t), const std::string& what) : m_id(id), m_close(close) {
        if (id < 0) fail("HDF5 call failed: " + what);
    }
    ~Hid() {
        if (m_id >= 0) m_close(m_id);
    }
    Hid(const Hid&) = delete;
    Hid& operator=(const Hid&) = delete;
    Hid(Hid&& other) noexcept : m_id(other.m_id), m_close(other.m_close) { other.m_id = -1; }
    Hid& operator=(Hid&&) = delete;
    operator hid_t() const { return m_id; }

   private:
    hid_t m_id;
    herr_t (*m_close)(hid_t);
};

// HDF5 prints its error stack by default; failures are reported as exceptions here instead.
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

Hid createGroup(hid_t loc, const char* name) {
    return Hid(H5Gcreate2(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose, name);
}
Hid openGroup(hid_t loc, const std::string& name) {
    return Hid(H5Gopen2(loc, name.c_str(), H5P_DEFAULT), H5Gclose, name);
}
// Whether the relative path `name` (components separated by '/') exists below `loc`.
bool exists(hid_t loc, const std::string& name) {
    std::size_t end = 0;
    do {
        end = name.find('/', end + 1);
        if (H5Lexists(loc, name.substr(0, end).c_str(), H5P_DEFAULT) <= 0) return false;
    } while (end != std::string::npos);
    return true;
}

// ---------- datasets ----------

std::vector<hsize_t> datasetDims(hid_t ds) {
    Hid space(H5Dget_space(ds), H5Sclose, "dataspace");
    const int rank = H5Sget_simple_extent_ndims(space);
    std::vector<hsize_t> dims(static_cast<std::size_t>(std::max(rank, 0)));
    if (rank > 0) H5Sget_simple_extent_dims(space, dims.data(), nullptr);
    return dims;
}

class Writer {
   public:
    explicit Writer(const WriteOptions& options) : m_options(options) {}

    // Dataset of prod(dims) values; chunked and compressed from 1024 values on.
    void dataset(hid_t loc, const char* name, hid_t type, const void* data, const std::vector<hsize_t>& dims) const {
        hsize_t count = 1;
        for (hsize_t d : dims) count *= d;
        Hid space(H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr), H5Sclose, name);
        Hid plist(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, name);
        if (m_options.compression > 0 && count >= 1024 && H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0) {
            std::vector<hsize_t> chunk = dims;
            const hsize_t inner = count / dims[0];
            chunk[0] = std::max<hsize_t>(1, std::min<hsize_t>(dims[0], 65536 / std::max<hsize_t>(1, inner)));
            check(H5Pset_chunk(plist, static_cast<int>(chunk.size()), chunk.data()), name);
            check(H5Pset_shuffle(plist), name);
            check(H5Pset_deflate(plist, static_cast<unsigned>(std::min(m_options.compression, 9))), name);
        }
        Hid ds(H5Dcreate2(loc, name, type, space, H5P_DEFAULT, plist, H5P_DEFAULT), H5Dclose, name);
        if (count > 0) check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), name);
    }
    void doubles(hid_t loc, const char* name, const double* x, std::size_t n) const { dataset(loc, name, H5T_NATIVE_DOUBLE, x, {static_cast<hsize_t>(n)}); }
    void doubles(hid_t loc, const char* name, const Vector& x) const { doubles(loc, name, x.data(), static_cast<std::size_t>(x.size())); }
    void ints(hid_t loc, const char* name, const int* x, std::size_t n) const { dataset(loc, name, H5T_NATIVE_INT, x, {static_cast<hsize_t>(n)}); }
    void int64s(hid_t loc, const char* name, const std::vector<std::int64_t>& x) const { dataset(loc, name, H5T_NATIVE_INT64, x.data(), {static_cast<hsize_t>(x.size())}); }
    // One int in a dataset of size 1, as fclib stores sizes.
    void scalar(hid_t loc, const char* name, int value) const { ints(loc, name, &value, 1); }

    // Fixed-length, null-terminated scalar string (H5LTmake_dataset_string).
    void string(hid_t loc, const char* name, const std::string& value) const {
        Hid type(H5Tcopy(H5T_C_S1), H5Tclose, name);
        check(H5Tset_size(type, value.size() + 1), name);
        Hid space(H5Screate(H5S_SCALAR), H5Sclose, name);
        Hid ds(H5Dcreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose, name);
        check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value.c_str()), name);
    }

    // A matrix in fclib's layout, compressed columns.
    void matrix(hid_t loc, const char* name, const SparseMatrix& A) const {
        SparseMatrix B = A;
        B.makeCompressed();
        if (B.nonZeros() > INT_MAX || B.rows() > INT_MAX || B.cols() >= INT_MAX) fail(std::string(name) + " is too large for fclib's 32-bit indices");
        Hid g = createGroup(loc, name);
        scalar(g, "nzmax", static_cast<int>(B.nonZeros()));
        scalar(g, "m", static_cast<int>(B.rows()));
        scalar(g, "n", static_cast<int>(B.cols()));
        scalar(g, "nz", -1);
        ints(g, "p", B.outerIndexPtr(), static_cast<std::size_t>(B.cols() + 1));
        ints(g, "i", B.innerIndexPtr(), static_cast<std::size_t>(B.nonZeros()));
        doubles(g, "x", B.valuePtr(), static_cast<std::size_t>(B.nonZeros()));
    }

   private:
    WriteOptions m_options;
};

template <typename T>
std::vector<T> readDataset(hid_t loc, const std::string& name, hid_t memtype) {
    Hid ds(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), H5Dclose, name);
    hsize_t count = 1;
    for (hsize_t d : datasetDims(ds)) count *= d;
    std::vector<T> v(static_cast<std::size_t>(count));
    if (count > 0) check(H5Dread(ds, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()), name);
    return v;
}
Vector readVector(hid_t loc, const std::string& name) {
    const auto v = readDataset<double>(loc, name, H5T_NATIVE_DOUBLE);
    return Eigen::Map<const Vector>(v.data(), static_cast<Eigen::Index>(v.size()));
}
std::vector<std::int64_t> readInts(hid_t loc, const std::string& name) {
    return readDataset<std::int64_t>(loc, name, H5T_NATIVE_INT64);
}
std::int64_t readScalar(hid_t loc, const std::string& name) {
    const auto v = readInts(loc, name);
    if (v.empty()) fail(name + " is empty");
    return v[0];
}
std::string readString(hid_t loc, const std::string& name) {
    Hid ds(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), H5Dclose, name);
    Hid ftype(H5Dget_type(ds), H5Tclose, name);
    if (H5Tis_variable_str(ftype) > 0) {
        Hid mtype(H5Tcopy(H5T_C_S1), H5Tclose, name);
        check(H5Tset_size(mtype, H5T_VARIABLE), name);
        char* c = nullptr;
        check(H5Dread(ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &c), name);
        std::string s = c ? c : "";
        H5free_memory(c);
        return s;
    }
    const std::size_t size = H5Tget_size(ftype);
    std::string s(size + 1, '\0');
    Hid mtype(H5Tcopy(H5T_C_S1), H5Tclose, name);
    check(H5Tset_size(mtype, size + 1), name);
    check(H5Dread(ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, s.data()), name);
    return s.substr(0, s.find('\0'));
}

// A matrix in fclib's layout (all three storage schemes); duplicate triplets are summed as in CSparse.
SparseMatrix readMatrix(hid_t loc, const std::string& name) {
    Hid g = openGroup(loc, name);
    const std::int64_t m = readScalar(g, "m"), n = readScalar(g, "n"), nz = readScalar(g, "nz");
    const auto p = readInts(g, "p");
    const auto i = readInts(g, "i");
    const auto x = readDataset<double>(g, "x", H5T_NATIVE_DOUBLE);
    std::vector<Eigen::Triplet<double>> t;
    auto add = [&](std::int64_t row, std::int64_t col, std::size_t k) {
        if (row < 0 || row >= m || col < 0 || col >= n) fail(name + ": index out of range");
        t.emplace_back(static_cast<int>(row), static_cast<int>(col), x[k]);
    };
    if (nz >= 0) {
        if (p.size() < (std::size_t)nz || i.size() < (std::size_t)nz || x.size() < (std::size_t)nz) fail(name + ": fewer triplets than nz");
        for (std::size_t k = 0; k < (std::size_t)nz; ++k) add(i[k], p[k], k);
    } else if (nz == -1 || nz == -2) {
        const std::int64_t outer = nz == -1 ? n : m;
        if ((std::int64_t)p.size() < outer + 1) fail(name + ": too few pointers");
        if (p[(size_t)outer] > (std::int64_t)std::min(i.size(), x.size())) fail(name + ": pointers beyond the data");
        for (std::int64_t j = 0; j < outer; ++j)
            for (std::int64_t k = p[(size_t)j]; k < p[(size_t)j + 1]; ++k) nz == -1 ? add(i[(size_t)k], j, (size_t)k) : add(j, i[(size_t)k], (size_t)k);
    } else {
        fail(name + ": unknown storage nz = " + std::to_string(nz));
    }
    SparseMatrix A(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(n));
    A.setFromTriplets(t.begin(), t.end());
    A.makeCompressed();
    return A;
}

// ---------- attributes (/cfc) ----------

// Variable-length UTF-8 string, which h5py reads as str.
Hid stringAttrType() {
    Hid t(H5Tcopy(H5T_C_S1), H5Tclose, "string type");
    check(H5Tset_size(t, H5T_VARIABLE), "string size");
    check(H5Tset_cset(t, H5T_CSET_UTF8), "string cset");
    return t;
}
template <typename T>
void writeAttr(hid_t loc, const char* name, hid_t type, const T& value) {
    Hid space(H5Screate(H5S_SCALAR), H5Sclose, "scalar space");
    Hid a(H5Acreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose, name);
    check(H5Awrite(a, type, &value), name);
}
void attr(hid_t loc, const char* name, const std::string& value) {
    const char* c = value.c_str();
    writeAttr(loc, name, stringAttrType(), c);
}
void attr(hid_t loc, const char* name, double value) {
    writeAttr(loc, name, H5T_NATIVE_DOUBLE, value);
}
void attr(hid_t loc, const char* name, std::int64_t value) {
    writeAttr(loc, name, H5T_NATIVE_INT64, value);
}

bool hasAttr(hid_t loc, const char* name) {
    return H5Aexists(loc, name) > 0;
}
std::string readStringAttr(hid_t loc, const char* name) {
    if (!hasAttr(loc, name)) return {};
    Hid a(H5Aopen(loc, name, H5P_DEFAULT), H5Aclose, name);
    Hid ftype(H5Aget_type(a), H5Tclose, name);
    if (H5Tis_variable_str(ftype) > 0) {
        char* c = nullptr;
        check(H5Aread(a, stringAttrType(), &c), name);
        std::string s = c ? c : "";
        H5free_memory(c);
        return s;
    }
    const std::size_t size = H5Tget_size(ftype);
    std::string s(size + 1, '\0');
    Hid mtype(H5Tcopy(H5T_C_S1), H5Tclose, name);
    check(H5Tset_size(mtype, size + 1), name);
    check(H5Aread(a, mtype, s.data()), name);
    return s.substr(0, s.find('\0'));
}
double readDoubleAttr(hid_t loc, const char* name, double fallback) {
    if (!hasAttr(loc, name)) return fallback;
    double v = 0;
    Hid a(H5Aopen(loc, name, H5P_DEFAULT), H5Aclose, name);
    check(H5Aread(a, H5T_NATIVE_DOUBLE, &v), name);
    return v;
}
std::int64_t readIntAttr(hid_t loc, const char* name, std::int64_t fallback) {
    if (!hasAttr(loc, name)) return fallback;
    std::int64_t v = 0;
    Hid a(H5Aopen(loc, name, H5P_DEFAULT), H5Aclose, name);
    check(H5Aread(a, H5T_NATIVE_INT64, &v), name);
    return v;
}

// ---------- fclib part ----------

// fclib's solution (v, l, u, r) with u = H^T v + w.
void writeSolution(const Writer& out, hid_t loc, const Problem& pb, const Solution& s) {
    if (s.v.size() != pb.n() || s.r.size() != 3 * pb.nc() || s.lambda.size() != pb.p()) throw std::invalid_argument("cfc file: a solution does not match the problem");
    Vector u = pb.w;
    u.noalias() += pb.H.transpose() * s.v;
    out.doubles(loc, "v", s.v);
    if (pb.p() > 0) out.doubles(loc, "l", s.lambda);
    out.doubles(loc, "u", u);
    out.doubles(loc, "r", s.r);
}

Solution readSolution(hid_t loc, const Problem& pb, const std::string& path) {
    Solution s;
    s.v = exists(loc, "v") ? readVector(loc, "v") : Vector::Zero(pb.n());
    s.r = readVector(loc, "r");
    s.lambda = exists(loc, "l") ? readVector(loc, "l") : Vector::Zero(pb.p());
    if (s.v.size() != pb.n() || s.r.size() != 3 * pb.nc() || s.lambda.size() != pb.p()) fail(path + ": a solution does not match the problem");
    return s;
}

// (A + A^T) / 2 if A is not symmetric; returns the relative asymmetry max|A - A^T| / max|A|.
double symmetrize(SparseMatrix& A) {
    const SparseMatrix At = A.transpose();
    const SparseMatrix D = A - At;
    double asym = 0, amax = 0;
    for (Eigen::Index j = 0; j < A.outerSize(); ++j)
        for (SparseMatrix::InnerIterator it(A, j); it; ++it) amax = std::max(amax, std::abs(it.value()));
    for (Eigen::Index j = 0; j < D.outerSize(); ++j)
        for (SparseMatrix::InnerIterator it(D, j); it; ++it) asym = std::max(asym, std::abs(it.value()));
    if (asym > 0) {
        A = (A + At) * 0.5;
        A.prune(0.0);
        A.makeCompressed();
    }
    return amax > 0 ? asym / amax : 0;
}

// ---------- metadata (/cfc) ----------

std::vector<double> flatten(const std::vector<Vector3>& x) {
    std::vector<double> out;
    out.reserve(3 * x.size());
    for (const auto& v : x) out.insert(out.end(), v.data(), v.data() + 3);
    return out;
}
// Frames row by row: frame[c](i, j) at index 9 c + 3 i + j (C order of an nc x 3 x 3 array).
std::vector<double> flatten(const std::vector<Matrix3>& x) {
    std::vector<double> out;
    out.reserve(9 * x.size());
    for (const auto& F : x)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out.push_back(F(i, j));
    return out;
}

void writeMetadata(const Writer& out, hid_t file, const Record& rec) {
    const Metadata& md = rec.meta;
    Hid cfc = createGroup(file, "cfc");
    attr(cfc, "format", std::string(kFormatName));
    attr(cfc, "version", std::string(kFormatVersion));
    attr(cfc, "producer", md.producer);
    attr(cfc, "producer_version", md.producer_version);
    attr(cfc, "scene", md.scene);
    attr(cfc, "step", md.step);
    attr(cfc, "time", md.time);
    attr(cfc, "h", md.h);
    attr(cfc, "theta", md.theta);
    attr(cfc, "units", md.units);

    auto ints = [&](hid_t g, const char* name, const std::vector<std::int64_t>& v) {
        if (!v.empty()) out.int64s(g, name, v);
    };
    auto doubles = [&](hid_t g, const char* name, const Vector& v) {
        if (v.size() > 0) out.doubles(g, name, v);
    };
    {
        Hid g = createGroup(cfc, "contacts");
        attr(g, "convention", std::string("per contact (column triple of H): frame columns N, T1, T2 in the world frame, N from body_a to body_b, r acts on body_b; body -1: static"));
        ints(g, "body_a", md.contact_body_a);
        ints(g, "body_b", md.contact_body_b);
        if (!md.contact_point.empty()) out.dataset(g, "point", H5T_NATIVE_DOUBLE, flatten(md.contact_point).data(), {md.contact_point.size(), 3});
        if (!md.contact_frame.empty()) out.dataset(g, "frame", H5T_NATIVE_DOUBLE, flatten(md.contact_frame).data(), {md.contact_frame.size(), 3, 3});
        doubles(g, "gap", md.contact_gap);
        doubles(g, "restitution_normal", md.restitution_normal);
        doubles(g, "restitution_tangential", md.restitution_tangential);
    }
    {
        Hid g = createGroup(cfc, "bilateral");
        attr(g, "convention", std::string("per column of G: kind 0 spring (position level), 1 damper; C_kk = law_compliance * law_scale; body -1: static"));
        ints(g, "kind", md.bilateral_kind);
        ints(g, "constraint", md.bilateral_constraint);
        ints(g, "constraint_row", md.bilateral_constraint_row);
        ints(g, "body_a", md.bilateral_body_a);
        ints(g, "body_b", md.bilateral_body_b);
        doubles(g, "law_compliance", md.bilateral_law_compliance);
        doubles(g, "law_scale", md.bilateral_law_scale);
    }
    {
        Hid g = createGroup(cfc, "bodies");
        attr(g, "convention", std::string("per block of v: type 0 rigid body (v, omega in the body frame), 1 point mass, 2 other"));
        ints(g, "dof_offset", md.body_dof_offset);
        ints(g, "dof_count", md.body_dof_count);
        ints(g, "type", md.body_type);
    }
    if (rec.reference) {
        Hid g = createGroup(cfc, "reference");
        attr(g, "solver", md.solver);
        attr(g, "law", std::string(toString(md.law)));
        attr(g, "iterations", md.iterations);
        attr(g, "solve_seconds", md.solve_seconds);
        attr(g, "natural_map", md.producer_natural_map);
        attr(g, "natural_map_coulomb", md.producer_natural_map_coulomb);
        attr(g, "bilateral", md.producer_bilateral);
    }
    Hid g = createGroup(cfc, "sequence");
    attr(g, "previous", md.previous_file);
}

void readMetadata(hid_t file, const std::string& path, const Problem& pb, Metadata& md) {
    Hid cfc = openGroup(file, "cfc");
    if (readStringAttr(cfc, "format") != kFormatName) fail(path + ": /cfc is not a cfc metadata group");
    const std::string version = readStringAttr(cfc, "version");
    if (version != kFormatVersion) fail(path + ": cfc version " + version + " is not supported (reader " + kFormatVersion + ")");
    md.present = true;
    md.producer = readStringAttr(cfc, "producer");
    md.producer_version = readStringAttr(cfc, "producer_version");
    md.scene = readStringAttr(cfc, "scene");
    md.step = readIntAttr(cfc, "step", -1);
    md.time = readDoubleAttr(cfc, "time", 0);
    md.h = readDoubleAttr(cfc, "h", 0);
    md.theta = readDoubleAttr(cfc, "theta", 0);
    md.units = readStringAttr(cfc, "units");

    // Optional per-item arrays: read if present, checked against the expected count.
    auto ints = [&](const std::string& name, std::size_t count, std::vector<std::int64_t>& v) {
        if (!exists(cfc, name)) return;
        v = readInts(cfc, name);
        if (v.size() != count) fail(path + ": /cfc/" + name + " has a wrong size");
    };
    auto doubles = [&](const std::string& name, Eigen::Index count, Vector& v) {
        if (!exists(cfc, name)) return;
        v = readVector(cfc, name);
        if (v.size() != count) fail(path + ": /cfc/" + name + " has a wrong size");
    };
    const auto nc = static_cast<std::size_t>(pb.nc()), p = static_cast<std::size_t>(pb.p());
    ints("contacts/body_a", nc, md.contact_body_a);
    ints("contacts/body_b", nc, md.contact_body_b);
    if (exists(cfc, "contacts/point")) {
        const Vector point = readVector(cfc, "contacts/point");
        if ((std::size_t)point.size() != 3 * nc) fail(path + ": /cfc/contacts/point has a wrong size");
        for (std::size_t c = 0; c < nc; ++c) md.contact_point.emplace_back(point.segment<3>(3 * (Eigen::Index)c));
    }
    if (exists(cfc, "contacts/frame")) {
        const Vector frame = readVector(cfc, "contacts/frame");
        if ((std::size_t)frame.size() != 9 * nc) fail(path + ": /cfc/contacts/frame has a wrong size");
        for (std::size_t c = 0; c < nc; ++c) {
            Matrix3 F;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) F(i, j) = frame[9 * (Eigen::Index)c + 3 * i + j];
            md.contact_frame.push_back(F);
        }
    }
    doubles("contacts/gap", pb.nc(), md.contact_gap);
    doubles("contacts/restitution_normal", pb.nc(), md.restitution_normal);
    doubles("contacts/restitution_tangential", pb.nc(), md.restitution_tangential);

    ints("bilateral/kind", p, md.bilateral_kind);
    ints("bilateral/constraint", p, md.bilateral_constraint);
    ints("bilateral/constraint_row", p, md.bilateral_constraint_row);
    ints("bilateral/body_a", p, md.bilateral_body_a);
    ints("bilateral/body_b", p, md.bilateral_body_b);
    doubles("bilateral/law_compliance", pb.p(), md.bilateral_law_compliance);
    doubles("bilateral/law_scale", pb.p(), md.bilateral_law_scale);

    if (exists(cfc, "bodies/dof_offset")) md.body_dof_offset = readInts(cfc, "bodies/dof_offset");
    ints("bodies/dof_count", md.body_dof_offset.size(), md.body_dof_count);
    ints("bodies/type", md.body_dof_offset.size(), md.body_type);

    if (exists(cfc, "reference")) {
        Hid g = openGroup(cfc, "reference");
        md.solver = readStringAttr(g, "solver");
        md.law = readStringAttr(g, "law") == "convexified" ? Law::Convexified : Law::Coulomb;
        md.iterations = readIntAttr(g, "iterations", -1);
        md.solve_seconds = readDoubleAttr(g, "solve_seconds", 0);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        md.producer_natural_map = readDoubleAttr(g, "natural_map", nan);
        md.producer_natural_map_coulomb = readDoubleAttr(g, "natural_map_coulomb", nan);
        md.producer_bilateral = readDoubleAttr(g, "bilateral", nan);
    }
    if (exists(cfc, "sequence")) md.previous_file = readStringAttr(openGroup(cfc, "sequence"), "previous");
}

}  // namespace

void write(const std::string& path, const Record& rec, const WriteOptions& options) {
    const Problem& pb = rec.problem;
    validate(pb);
    if (pb.nc() == 0) throw std::invalid_argument("cfc file: problems without contacts cannot be written (fclib needs contact rows)");
    const Writer out(options);
    const QuietErrors quiet;
    std::error_code ec;
    std::filesystem::remove(path, ec);
    Hid file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose, "create " + path);

    {
        Hid g = createGroup(file, "fclib_global");
        out.scalar(g, "spacedim", 3);
        out.matrix(g, "M", pb.M);
        out.matrix(g, "H", pb.H);
        if (pb.p() > 0) out.matrix(g, "G", pb.G);
        if (pb.hasCompliance()) out.matrix(g, "C", pb.C);
        {
            Hid v = createGroup(g, "vectors");
            out.doubles(v, "f", pb.f);
            out.doubles(v, "w", pb.w);
            out.doubles(v, "mu", pb.mu);
            if (pb.p() > 0) out.doubles(v, "b", pb.b);
            out.doubles(v, "desaxce_shift", pb.shift);
        }
        Hid info = createGroup(g, "info");
        const Metadata& md = rec.meta;
        out.string(info, "title", !rec.title.empty() ? rec.title : (md.scene.empty() ? std::string("cfc problem") : md.scene + " step " + std::to_string(md.step)));
        std::string description = "Global frictional contact problem";
        if (!md.producer.empty()) description += " from " + md.producer + (md.producer_version.empty() ? "" : " " + md.producer_version);
        if (md.h > 0) description += ", one time step h = " + std::to_string(md.h) + ", theta = " + std::to_string(md.theta);
        if (rec.reference && !md.solver.empty()) description += "; reference solution by " + md.solver + " (" + toString(md.law) + " law)";
        out.string(info, "description", description + ". Extensions (format cfc " + kFormatVersion + "): /fclib_global/C, /fclib_global/vectors/desaxce_shift, /cfc.");
        out.string(info, "math_info",
                   pb.hasCompliance() ? "M symmetric positive definite; compliant bilateral rows G^T v + b + C lambda = 0 with C in /fclib_global/C" : "M symmetric positive definite");
    }
    if (rec.reference) {
        Hid g = createGroup(file, "solution");
        writeSolution(out, g, pb, *rec.reference);
    }
    if (rec.guess) {
        Hid g = createGroup(file, "guesses");
        out.scalar(g, "number_of_guesses", 1);
        Hid one = createGroup(g, "1");
        writeSolution(out, one, pb, *rec.guess);
    }
    if (rec.meta.present) writeMetadata(out, file, rec);
}

Record read(const std::string& path) {
    const QuietErrors quiet;
    if (!std::filesystem::is_regular_file(path)) fail(path + " does not exist");
    Hid file(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose, "open " + path);
    if (!exists(file, "fclib_global")) fail(path + " holds no fclib global problem");
    Hid g = openGroup(file, "fclib_global");
    if (exists(g, "spacedim") && readScalar(g, "spacedim") != 3) fail(path + ": only 3D problems (spacedim = 3) are supported");

    Record rec;
    Problem& pb = rec.problem;
    pb.M = readMatrix(g, "M");
    pb.H = readMatrix(g, "H");
    pb.G = exists(g, "G") ? readMatrix(g, "G") : SparseMatrix(pb.M.rows(), 0);
    pb.f = readVector(g, "vectors/f");
    pb.w = readVector(g, "vectors/w");
    pb.mu = readVector(g, "vectors/mu");
    pb.b = exists(g, "vectors/b") ? readVector(g, "vectors/b") : Vector::Zero(pb.G.cols());
    pb.shift = exists(g, "vectors/desaxce_shift") ? readVector(g, "vectors/desaxce_shift") : Vector::Zero(pb.mu.size());
    pb.C = exists(g, "C") ? readMatrix(g, "C") : SparseMatrix(pb.b.size(), pb.b.size());
    if (exists(g, "info/title")) rec.title = readString(g, "info/title");
    rec.mass_asymmetry = symmetrize(pb.M);
    symmetrize(pb.C);
    try {
        validate(pb);
    } catch (const std::invalid_argument& e) {
        fail(path + ": " + e.what());
    }

    if (exists(file, "solution")) rec.reference = readSolution(openGroup(file, "solution"), pb, path);
    if (exists(file, "guesses/1")) rec.guess = readSolution(openGroup(file, "guesses/1"), pb, path);
    if (exists(file, "cfc")) readMetadata(file, path, pb, rec.meta);
    return rec;
}

}  // namespace cfc
