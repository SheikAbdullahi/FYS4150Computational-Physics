#include "PenningTrap.h"
#include "Particle.h"
#include <iostream>

int main() {
    // Create a Penning trap instance
    PenningTrap myTrap(1e4, 1e8, 1e2);  // Example values for B0, V0, and d

    // Create some particles
    Particle p1(1.602e-19, 9.109e-31, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    Particle p2(-1.602e-19, 9.109e-31, {0.1, 0.0, 0.0}, {0.0, 0.0, 0.0});

    // Add particles to the trap
    myTrap.add_particle(p1);
    myTrap.add_particle(p2);

    // Calculate and display forces, for example, the total force on the first particle
    arma::vec force = myTrap.total_force(0);
    force.print("Total force on particle 1:");

    // Evolve the system using Forward Euler
    myTrap.evolve_forward_euler(1e-5);  // for a 0.00001 second time step

    // Evolve the system using RK4
    myTrap.evolve_RK4(1e-5);  // for a 0.00001 second time step

    return 0;
}
