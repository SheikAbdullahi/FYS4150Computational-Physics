import numpy as np

class Particle:
    def __init__(self, position, velocity, mass, charge):
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.mass = mass
        self.charge = charge

    def update_position(self, dt):
        self.position += self.velocity * dt

    def update_velocity(self, acceleration, dt):
        self.velocity += acceleration * dt


class PenningTrap:
    def __init__(self, B0, V0, d, particles):
        self.B0 = B0
        self.V0 = V0
        self.d = d
        self.particles = particles

    def magnetic_field(self):
        return np.array([0, 0, self.B0])

    def electric_field(self, particle):
        x, y, z = particle.position
        return np.array([x, y, -(2 * z)]) * self.V0 / (self.d ** 2)

    def force_due_to_fields(self, particle):
        q = particle.charge
        v = particle.velocity
        E = self.electric_field(particle)
        B = self.magnetic_field()
        return q * (E + np.cross(v, B))

    def total_force(self, particle):
        # Assuming no interactions between particles for simplicity
        return self.force_due_to_fields(particle)

    def acceleration(self, particle):
        return self.total_force(particle) / particle.mass

    def update(self, dt):
        for particle in self.particles:
            acceleration = self.acceleration(particle)
            particle.update_velocity(acceleration, dt)
            particle.update_position(dt)
