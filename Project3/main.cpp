// main.cpp
#include <iostream>
#include <vector>
#include "Particle.h"
#include "PenningTrap.h"

int main() {
    // Constants for the PenningTrap
    double B0 = 9.65e1;  // magnetic field strength
    double V0 = 9.65e6;  // electric potential
    double d = 500;      // characteristic dimension
    double ke = 1.38935333e5;  // Coulomb constant

    // Create a few particles
    Particle p1(1, 1, {0, 0, 0}, {0, 0, 0});  // A particle with charge 1e, mass 1u, position (0, 0, 0), and velocity (0, 0, 0)
    Particle p2(1, 1, {10, 0, 0}, {0, 1, 0});  // Another particle, different position and velocity

    // Store them in a vector
    std::vector<Particle> particles = {p1, p2};

    // Create a PenningTrap with those particles
    PenningTrap trap(B0, V0, d, ke, particles);

    // Time evolution parameters
    double dt = 0.1;  // time step
    int steps = 1000;  // number of steps to simulate

    // Simulate some time evolution
    for (int i = 0; i < steps; ++i) {
        trap.updateParticles(dt);

        // Output the state at each step, or every few steps, as needed
        std::cout << "After step " << i+1 << ":" << std::endl;
        for (size_t j = 0; j < particles.size(); ++j) {
            std::cout << "  Particle " << j+1 << " position: " 
                      << particles[j].position[0] << ", "
                      << particles[j].position[1] << ", "
                      << particles[j].position[2] << std::endl;
        }
    }

    return 0;
}
