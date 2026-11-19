#pragma once

#include <memory>

#include "problem.hpp"

namespace cfc {

// Applies M^-1 (diagonal M directly, otherwise through a sparse LDL^T factorization).
class MassInverse {
   public:
    explicit MassInverse(const SparseMatrix& M);
    ~MassInverse();
    Vector apply(const Vector& x) const;
    SparseMatrix apply(const SparseMatrix& X) const;

   private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

// Residuals of a candidate solution (v, r, lambda), the common accuracy measure of all solvers
// (FORMAT.md, "Residuals"; the Python twin is python/cfc/residuals.py). All values are in velocity
// units and take the maximum over rows or contacts, so that they do not depend on the masses:
//
//   balance      ||v - M^-1 (H r + G lambda + f)||_inf
//   bilateral    ||G^T v + b + C lambda||_inf
//   natural_map  max_c ||r_c - proj_K(r_c - rho_c u_hat_c)|| / rho_c, u_hat with the shift of the law
//   cone_primal  max_c dist(r_c, K_mu) / rho_c
//   cone_dual    max_c dist(u_hat_c, K_mu*)
//   complementarity  max_c sqrt(|r_c^T u_hat_c| / rho_c)
//
// with rho_c = 1 / lambda_max(H_c^T M^-1 H_c) per contact (contacts with H_c = 0 skipped), the norm
// max(|x_N|, ||x_T||) of the natural map and the Euclidean norm for the distances. `natural_map` uses
// the given frozen shift (convexified law), `natural_map_coulomb` the shift mu ||u_T|| of Coulomb's law;
// the cone and complementarity terms are those of the convexified law (the KKT conditions of the cone
// program). velocity_scale = max_c ||H_c^T M^-1 f + w_c||, the largest contact velocity of the motion
// without contact and bilateral forces (1 if it vanishes), turns them into relative residuals.
struct Residuals {
    double balance{0};
    double bilateral{0};
    double natural_map{0};
    double natural_map_coulomb{0};
    double cone_primal{0};
    double cone_dual{0};
    double complementarity{0};
    double velocity_scale{1};

    // Largest residual of the convexified law: balance, bilateral and its natural map; NaN if any is NaN.
    double convexifiedMax() const;
    // The same for Coulomb's law.
    double coulombMax() const;
};

// `shift`: the frozen De Saxce shift per contact of the convexified law (Problem::shift, or zero for the
// relaxation without shift).
Residuals residuals(const Problem& pb, const Solution& sol, const Vector& shift);

}  // namespace cfc
