// Particle.cpp
#include "Particle.h"
#include <iostream>

// Constructor implementation
Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity) : q(charge), m(mass), r(position), v(velocity) {}

// Method to update the particle's position with time step dt
void Particle::move(double dt) {
    r += v * dt;  // Simple Euler integration for motion
}

// Method to print particle's properties
void Particle::print() {
    r.print("Position:");
    v.print("Velocity:");
    std::cout << "Charge: " << q << "\nMass: " << m << std::endl;
}
