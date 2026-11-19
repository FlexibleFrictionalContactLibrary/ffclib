#include <cfc/problem.hpp>

#include <stdexcept>
#include <string>

namespace cfc {

const char* toString(Law law) {
    return law == Law::Coulomb ? "coulomb" : "convexified";
}

void validate(const Problem& pb) {
    const Eigen::Index n = pb.M.rows(), nc = pb.mu.size(), p = pb.b.size();
    auto require = [](bool ok, const char* what) {
        if (!ok) throw std::invalid_argument(std::string("cfc problem: ") + what);
    };
    require(pb.M.cols() == n, "M must be square");
    require(pb.H.rows() == n && pb.H.cols() == 3 * nc, "H must be n x 3 nc");
    require(pb.G.rows() == n && pb.G.cols() == p, "G must be n x p");
    require(pb.C.rows() == p && pb.C.cols() == p, "C must be p x p");
    require(pb.f.size() == n, "f must have n entries");
    require(pb.w.size() == 3 * nc, "w must have 3 nc entries");
    require(pb.shift.size() == nc, "shift must have nc entries");
    if (pb.hasCompliance()) {
        const SparseMatrix Ct = pb.C.transpose();
        require((pb.C - Ct).norm() == 0, "C must be symmetric");
    }
}

}  // namespace cfc
