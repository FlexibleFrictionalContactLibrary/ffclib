#pragma once

#include <cmath>

#include "types.hpp"

// The Coulomb cone of one contact in the component order (N, T1, T2):
//
//   K_mu  = { x : ||x_T|| <= mu x_N }          (percussions)
//   K_mu* = { v : mu ||v_T|| <= v_N }          (shifted contact velocities)
//
// The polar cone of K_mu is -K_mu*, so x = proj_K(x) + proj_{-K*}(x) with orthogonal parts (Moreau).
namespace cfc {

// proj_{K_mu}(x). The polar test comes first: at mu = 0 it catches every point with x_N <= 0.
inline Vector3 projectCoulombCone(const Vector3& x, double mu) {
    const double tNorm = std::sqrt(x[1] * x[1] + x[2] * x[2]);
    if (mu * tNorm <= -x[0]) return Vector3::Zero();
    if (tNorm <= mu * x[0]) return x;
    // Boundary point (r, mu r e_T) closest to x, r = (x_N + mu ||x_T||) / (1 + mu^2); tNorm > 0 here.
    const double r = (x[0] + mu * tNorm) / (1.0 + mu * mu);
    const double scaleT = mu * r / tNorm;
    return Vector3(r, scaleT * x[1], scaleT * x[2]);
}

// Normal component mu ||u_T|| of the De Saxce shift of Coulomb's law.
inline double deSaxce(const Vector3& u, double mu) {
    return mu * std::sqrt(u[1] * u[1] + u[2] * u[2]);
}

}  // namespace cfc
