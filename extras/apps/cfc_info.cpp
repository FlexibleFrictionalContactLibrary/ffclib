// cfc_info: one JSON line per cfc or fclib file with its sizes, metadata and the residuals of the stored
// reference solution (FORMAT.md, "Residuals").
//
//   cfc_info [--check] FILE...
//
// --check: exit code 1 if the residuals evaluated from a file differ from the producer's own
// (/cfc/reference), the round-trip test of an export.

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <cfc/file.hpp>
#include <cfc/json.hpp>
#include <cfc/residuals.hpp>

namespace {

std::string problemName(const std::string& path) {
    std::string stem = std::filesystem::path(path).filename().string();
    for (const char* ext : {".cfc.h5", ".hdf5", ".h5"})
        if (stem.size() > std::strlen(ext) && stem.compare(stem.size() - std::strlen(ext), std::string::npos, ext) == 0) return stem.substr(0, stem.size() - std::strlen(ext));
    return stem;
}

// Bilateral rows with a nonzero column in C.
long long compliantRows(const cfc::SparseMatrix& C) {
    long long count = 0;
    for (Eigen::Index k = 0; k < C.outerSize(); ++k)
        for (cfc::SparseMatrix::InnerIterator it(C, k); it; ++it)
            if (it.value() != 0) {
                ++count;
                break;
            }
    return count;
}

// Agreement to 1e-9 relative, or to 1e-12 of the magnitude of the terms of the residual: residuals at
// round-off level are summed in a different order by the producer.
bool agree(double a, double b, double terms) {
    return std::isnan(b) || std::abs(a - b) <= 1e-9 * std::max(std::abs(a), std::abs(b)) + 1e-12 * terms;
}

}  // namespace

int main(int argc, char** argv) {
    bool check = false;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--check")
            check = true;
        else if (a == "-h" || a == "--help" || (!a.empty() && a[0] == '-')) {
            std::cout << "usage: cfc_info [--check] FILE...\n";
            return a == "-h" || a == "--help" ? 0 : 2;
        } else
            files.push_back(a);
    }
    if (files.empty()) {
        std::cout << "usage: cfc_info [--check] FILE...\n";
        return 2;
    }

    int failures = 0;
    for (const auto& path : files) {
        cfc::json::Object o;
        o.str("problem", problemName(path)).str("file", path);
        cfc::Record rec;
        try {
            rec = cfc::read(path);
        } catch (const std::exception& e) {
            o.str("error", e.what());
            std::cout << o.dump() << std::endl;
            ++failures;
            continue;
        }
        const cfc::Problem& pb = rec.problem;
        const cfc::Metadata& md = rec.meta;
        o.str("title", rec.title).str("family", md.present ? md.scene : std::filesystem::path(path).parent_path().filename().string()).integer("step", md.step);
        o.str("producer", md.producer).str("producer_version", md.producer_version);
        o.integer("n", pb.n()).integer("nc", pb.nc()).integer("m", 3 * pb.nc()).integer("p", pb.p());
        o.integer("n_compliant", compliantRows(pb.C)).integer("n_frictionless", (pb.mu.array() == 0).count());
        o.integer("nnz_M", pb.M.nonZeros()).integer("nnz_H", pb.H.nonZeros()).integer("nnz_G", pb.G.nonZeros()).integer("nnz_C", pb.C.nonZeros());
        o.num("mass_asymmetry", rec.mass_asymmetry).boolean("has_metadata", md.present).boolean("has_reference", rec.reference.has_value());
        if (rec.reference) {
            const cfc::Residuals r = cfc::residuals(pb, *rec.reference, pb.shift);
            o.str("solver", md.solver).str("law", cfc::toString(md.law)).integer("iterations", md.iterations).num("solve_seconds", md.solve_seconds);
            o.obj("residuals", cfc::json::residuals(r));
            cfc::json::Object producer;
            producer.num("natural_map", md.producer_natural_map).num("natural_map_coulomb", md.producer_natural_map_coulomb).num("bilateral", md.producer_bilateral);
            o.obj("producer_residuals", producer);
            // The producer's natural map is that of the law its solver solved.
            const double own = md.law == cfc::Law::Coulomb ? r.natural_map_coulomb : r.natural_map;
            cfc::Vector bilateralTerms = pb.G.transpose() * rec.reference->v;
            double terms = std::max(bilateralTerms.lpNorm<Eigen::Infinity>(), pb.b.lpNorm<Eigen::Infinity>());
            if (pb.hasCompliance()) terms = std::max(terms, (pb.C * rec.reference->lambda).lpNorm<Eigen::Infinity>());
            const bool ok = agree(own, md.producer_natural_map, r.velocity_scale) && agree(r.natural_map_coulomb, md.producer_natural_map_coulomb, r.velocity_scale) &&
                            agree(r.bilateral, md.producer_bilateral, terms);
            o.boolean("agree", ok);
            failures += check && !ok ? 1 : 0;
        }
        std::cout << o.dump() << std::endl;
    }
    return failures > 0 ? 1 : 0;
}
