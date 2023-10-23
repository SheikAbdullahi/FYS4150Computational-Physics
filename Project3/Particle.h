// Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H

#include <armadillo>

class Particle {
public:
    // Member variables
    double q;  // charge
    double m;  // mass
    arma::vec r;  // position
    arma::vec v;  // velocity

    // Constructor
    Particle(double charge, double mass, arma::vec position, arma::vec velocity);

    // Method to update the particle's position with time step dt
    void move(double dt);

    // Method to print particle's properties
    void print();
};

#endif
