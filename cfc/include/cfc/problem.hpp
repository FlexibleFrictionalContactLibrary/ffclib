#pragma once

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "types.hpp"

// The global frictional contact problem of fclib, extended by compliant bilateral rows (FORMAT.md):
//
//   M v = H r + G lambda + f
//   G^T v + b + C lambda = 0                          C = 0: perfect constraints (plain fclib)
//   u = H^T v + w,  u_hat = u + (s, 0, 0) per contact
//   K_mu* ∋ u_hat ⊥ r ∈ K_mu                          per contact, local frame (N, T1, T2)
//
// with K_mu = {r : ||r_T|| <= mu r_N}. Coulomb's law has the De Saxce shift s = mu ||u_T||; the
// convexified law freezes it at Problem::shift.
namespace cfc {

inline constexpr const char* kFormatName = "cfc";
inline constexpr const char* kFormatVersion = "0.2";

// Contact law that a solution satisfies (and that the residuals measure).
enum class Law { Coulomb, Convexified };
const char* toString(Law law);

struct Problem {
    SparseMatrix M;  // n x n, symmetric positive definite, both triangles stored
    SparseMatrix H;  // n x 3 nc, contact force directions, columns (N, T1, T2) per contact
    SparseMatrix G;  // n x p, bilateral force directions
    SparseMatrix C;  // p x p, symmetric positive semidefinite compliance of the bilateral rows, both
                     // triangles stored; no entries: perfect constraints
    Vector f;        // n
    Vector w;        // 3 nc
    Vector b;        // p
    Vector mu;       // nc; 0: frictionless contact
    // nc, the frozen De Saxce shift of the convexified law; zero if the file has none.
    Vector shift;

    int n() const { return static_cast<int>(M.rows()); }
    int nc() const { return static_cast<int>(mu.size()); }
    int p() const { return static_cast<int>(b.size()); }
    bool hasCompliance() const { return C.nonZeros() > 0; }
};

// Throws std::invalid_argument if the sizes of the problem do not match.
void validate(const Problem& pb);

struct Solution {
    Vector v;       // n
    Vector r;       // 3 nc
    Vector lambda;  // p
};

// What the producer knows beyond the problem (group /cfc, optional): provenance, contacts, bodies,
// bilateral rows and how the reference solution was computed. Default values when absent.
struct Metadata {
    bool present{false};  // the file has a /cfc group

    // Provenance
    std::string producer;          // program that wrote the file, e.g. "CardilloCxx"
    std::string producer_version;  // its version (git describe)
    std::string scene;
    std::int64_t step{-1};  // index of the solved time step (0: from t = 0 to h)
    double time{0};         // time at the start of the step
    double h{0};
    double theta{0};
    std::string units;

    // Per contact, in the column order of H. Bodies are blocks of v (-1: static partner).
    std::vector<std::int64_t> contact_body_a, contact_body_b;
    std::vector<Vector3> contact_point;  // world frame
    std::vector<Matrix3> contact_frame;  // columns N, T1, T2 (world frame), N from A to B
    Vector contact_gap;                  // signed distance (negative: penetration)
    Vector restitution_normal, restitution_tangential;

    // Per body (block of v): first index and size; type 0 rigid body (v, omega in the body frame),
    // 1 point mass, 2 other.
    std::vector<std::int64_t> body_dof_offset, body_dof_count, body_type;

    // Per bilateral row (column of G): kind 0 spring (position-level force law), 1 damper; the
    // constraint it belongs to and its row there; the constraint's bodies; the compliance of the force
    // law and the factor of the time discretization with C_kk = law_compliance_k * law_scale_k.
    std::vector<std::int64_t> bilateral_kind, bilateral_constraint, bilateral_constraint_row, bilateral_body_a, bilateral_body_b;
    Vector bilateral_law_compliance, bilateral_law_scale;

    // Reference solution: solver, law and the producer's own residuals (NaN if not evaluated).
    std::string solver;
    Law law{Law::Coulomb};
    std::int64_t iterations{-1};
    double solve_seconds{0};
    double producer_natural_map{std::numeric_limits<double>::quiet_NaN()};
    double producer_natural_map_coulomb{std::numeric_limits<double>::quiet_NaN()};
    double producer_bilateral{std::numeric_limits<double>::quiet_NaN()};

    // File of the previous time step of the same run, if it was written, else empty.
    std::string previous_file;
};

// Everything one file holds.
struct Record {
    Problem problem;
    Metadata meta;
    std::optional<Solution> reference;  // /solution
    std::optional<Solution> guess;      // /guesses/1, e.g. the solver's warm start
    // Set by read(): title of /fclib_global/info, and the relative asymmetry ||M - M^T||_max / ||M||_max
    // of the stored mass matrix, which read() replaces by (M + M^T) / 2 (FCLIB collection files store M
    // in single precision).
    std::string title;
    double mass_asymmetry{0};
};

}  // namespace cfc
