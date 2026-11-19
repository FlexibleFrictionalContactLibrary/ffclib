#pragma once

#include <cmath>

#include <Eigen/Core>

// The friction cone of one contact, components (N, T1, T2):
//
//   K_mu  = { x : ||x_T|| <= mu x_N }     percussions
//   K_mu* = { x : mu ||x_T|| <= x_N }     velocities (the dual cone)
//
// Its polar cone is -K_mu*, so x = proj_K(x) + proj_{-K*}(x) with orthogonal parts (Moreau).
namespace ffclib {

// proj_{K_mu}(x). The polar test comes first: at mu = 0 it catches every point with x_N <= 0.
inline Eigen::Vector3d projectOntoCone(const Eigen::Vector3d& x, double mu) {
    const double tNorm = std::hypot(x[1], x[2]);
    if (mu * tNorm <= -x[0]) return Eigen::Vector3d::Zero();
    if (tNorm <= mu * x[0]) return x;
    // The closest boundary point (r, mu r x_T / ||x_T||) with r = (x_N + mu ||x_T||) / (1 + mu^2); tNorm > 0 here.
    const double r = (x[0] + mu * tNorm) / (1 + mu * mu);
    const double scaleT = mu * r / tNorm;
    return {r, scaleT * x[1], scaleT * x[2]};
}

}  // namespace ffclib
