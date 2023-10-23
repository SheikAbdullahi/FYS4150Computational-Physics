// PenningTrap.h
#ifndef PENNINGTRAP_H
#define PENNINGTRAP_H

#include <vector>
#include "Particle.h"

class PenningTrap {
public:
    double B0;
    double V0;
    double d;
    double ke;
    std::vector<Particle> particles;

    PenningTrap(double B, double V, double dim, double k, const std::vector<Particle>& p);
    Eigen::Vector3d magneticForce(const Particle& particle);
    Eigen::Vector3d electricForce(const Particle& particle);
    Eigen::Vector3d coulombForce(const Particle& particle);
    Eigen::Vector3d totalForce(const Particle& particle);
    void updateParticles(double dt);
};

#endif
