// PenningTrap.h
#ifndef PENNINGTRAP_H
#define PENNINGTRAP_H

#include "Particle.h"
#include <vector>

class PenningTrap {
public:
    // Member variables
    double B0;  // Magnetic field strength
    double V0;  // Applied potential
    double d;   // Characteristic dimension
    std::vector<Particle> particles;  // Container for particles

    // Constructor
    PenningTrap(double B0, double V0, double d);

    // Method to add a particle to the trap
    void add_particle(Particle& p);

    // Method to calculate and return the external electric field at a point in space
    arma::vec external_electric_field(arma::vec r);

    // Method to calculate and return the external magnetic field at a point in space
    arma::vec external_magnetic_field(arma::vec r);

    // Method to calculate the force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // Method to calculate the total force on a particle from the external fields
    arma::vec total_force_external(int i);

    // Method to calculate the total force on a particle from other particles
    arma::vec total_force_particles(int i);

    // Method to calculate the total force on a particle from both external fields and other particles
    arma::vec total_force(int i);
    // Method to evolve the system one time step forward using the Forward Euler method
    void evolve_forward_euler(double dt);

    // Method to evolve the system one time step forward using the RK4 method
    void evolve_RK4(double dt);
};

#endif
