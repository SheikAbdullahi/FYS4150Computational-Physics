// Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Dense>

class Particle {
public:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double mass;
    double charge;

    Particle(const Eigen::Vector3d& pos, const Eigen::Vector3d& vel, double m, double q);
    void updatePosition(double dt);
    void updateVelocity(const Eigen::Vector3d& acceleration, double dt);
};

#endif
