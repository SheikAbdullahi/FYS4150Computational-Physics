// PenningTrap.cpp
#include "PenningTrap.h"
#include <cmath> // for std::pow and std::sqrt

PenningTrap::PenningTrap(double B, double V, double dim, double k, const std::vector<Particle>& p)
    : B0(B), V0(V), d(dim), ke(k), particles(p) {}

Eigen::Vector3d PenningTrap::magneticForce(const Particle& particle) {
    // Lorentz force only due to magnetic field: F = q(v x B)
    Eigen::Vector3d B_field(0, 0, B0); // Assuming uniform magnetic field in z-direction
    return particle.charge * particle.velocity.cross(B_field);
}

Eigen::Vector3d PenningTrap::electricForce(const Particle& particle) {
    // Electric force: F = qE
    // Assuming quadrupole electric field
    double E_x = V0 / (d*d) * particle.position[0];
    double E_y = V0 / (d*d) * particle.position[1];
    double E_z = -2 * V0 / (d*d) * particle.position[2];
    Eigen::Vector3d E_field(E_x, E_y, E_z);

    return particle.charge * E_field;
}

Eigen::Vector3d PenningTrap::coulombForce(const Particle& particle) {
    // Initializing total Coulomb force
    Eigen::Vector3d total_force(0, 0, 0);

    for (const auto& other : particles) {
        if (&other == &particle) continue; // no self-interaction

        Eigen::Vector3d r = particle.position - other.position;
        double r_magnitude = r.norm();

        if (r_magnitude > 0) { // To avoid division by zero
            // Coulomb's Law: F = ke * q1 * q2 * r / |r|^3
            total_force += ke * particle.charge * other.charge * r / std::pow(r_magnitude, 3);
        }
    }

    return total_force;
}

Eigen::Vector3d PenningTrap::totalForce(const Particle& particle) {
    // Summing up all forces
    return magneticForce(particle) + electricForce(particle) + coulombForce(particle);
}

void PenningTrap::updateParticles(double dt) {
    // Simple Euler's method is used for updating particles' positions and velocities
    // A more accurate method like Runge-Kutta might be needed for a precise simulation

    for (auto& particle : particles) {
        Eigen::Vector3d acceleration = totalForce(particle) / particle.mass;
        particle.updateVelocity(acceleration, dt);
        particle.updatePosition(dt);
    }
}
