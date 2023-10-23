// PenningTrap.cpp
#include "PenningTrap.h"

// Constants
const double ke = 8.9875517873681764e9;  // Coulomb's constant in vacuum

// Constructor implementation
PenningTrap::PenningTrap(double B0, double V0, double d) : B0(B0), V0(V0), d(d) {}

// Method to add a particle to the trap
void PenningTrap::add_particle(Particle& p) {
    particles.push_back(p);
}

// Method to calculate and return the external electric field at a point in space
arma::vec PenningTrap::external_electric_field(arma::vec r) {
    // Based on the formula for the quadrupole field in a Penning trap
    return {V0 / (d * d) * r[0], V0 / (d * d) * r[1], -2 * V0 / (d * d) * r[2]};
}

// Method to calculate and return the external magnetic field at a point in space
arma::vec PenningTrap::external_magnetic_field(arma::vec r) {
    // The magnetic field in a Penning trap is uniform
    return {0, 0, B0};
}

// Method to calculate the force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j) {
    arma::vec r = particles[i].r - particles[j].r;  // relative position vector
    double r_norm = arma::norm(r, 2);  // the Euclidean norm of r
    return ke * particles[i].q * particles[j].q / (r_norm * r_norm * r_norm) * r;  // Coulomb's law
}

// Method to calculate the total force on a particle from the external fields
arma::vec PenningTrap::total_force_external(int i) {
    arma::vec force_electric = particles[i].q * external_electric_field(particles[i].r);
    arma::vec force_magnetic = particles[i].q * arma::cross(particles[i].v, external_magnetic_field(particles[i].r));
    return force_electric + force_magnetic;
}

// Method to calculate the total force on a particle from other particles
arma::vec PenningTrap::total_force_particles(int i) {
    arma::vec total_force = {0, 0, 0};
    for (int j = 0; j < particles.size(); j++) {
        if (i != j) {
            total_force += force_particle(i, j);
        }
    }
    return total_force;
}

// Method to calculate the total force on a particle from both external fields and other particles
arma::vec PenningTrap::total_force(int i) {
    return total_force_external(i) + total_force_particles(i);
}

void PenningTrap::evolve_forward_euler(double dt) {
    std::vector<arma::vec> forces;
    for (int i = 0; i < particles.size(); i++) {
        forces.push_back(total_force(i));
    }

    for (int i = 0; i < particles.size(); i++) {
        particles[i].r += dt * particles[i].v;
        particles[i].v += dt * forces[i] / particles[i].m;
    }
}

void PenningTrap::evolve_RK4(double dt) {
    std::vector<Particle> original_particles = particles;  // Make a copy of all particles

    // Store k values
    std::vector<std::vector<arma::vec>> k_r(4, std::vector<arma::vec>(particles.size()));
    std::vector<std::vector<arma::vec>> k_v(4, std::vector<arma::vec>(particles.size()));

    for (int step = 0; step < 4; step++) {
        std::vector<arma::vec> forces;
        for (int i = 0; i < particles.size(); i++) {
            forces.push_back(total_force(i));
        }

        for (int i = 0; i < particles.size(); i++) {
            k_r[step][i] = dt * (step == 0 ? particles[i].v : particles[i].v + 0.5 * k_v[step - 1][i]);
            k_v[step][i] = dt * forces[i] / particles[i].m;

            // Update positions and velocities for the next k calculation
            particles[i].r = original_particles[i].r + (step == 2 ? k_r[step][i] : 0.5 * k_r[step][i]);
            particles[i].v = original_particles[i].v + (step == 2 ? k_v[step][i] : 0.5 * k_v[step][i]);
        }
    }

    // Final update
    for (int i = 0; i < particles.size(); i++) {
        particles[i].r = original_particles[i].r + (1.0/6.0) * (k_r[0][i] + 2.0 * k_r[1][i] + 2.0 * k_r[2][i] + k_r[3][i]);
        particles[i].v = original_particles[i].v + (1.0/6.0) * (k_v[0][i] + 2.0 * k_v[1][i] + 2.0 * k_v[2][i] + k_v[3][i]);
    }
}
