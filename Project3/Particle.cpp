// Particle.cpp
#include "Particle.h"

Particle::Particle(const Eigen::Vector3d& pos, const Eigen::Vector3d& vel, double m, double q)
    : position(pos), velocity(vel), mass(m), charge(q) {}

void Particle::updatePosition(double dt) {
    position += velocity * dt;
}

void Particle::updateVelocity(const Eigen::Vector3d& acceleration, double dt) {
    velocity += acceleration * dt;
}
