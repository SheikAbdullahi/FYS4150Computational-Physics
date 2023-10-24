import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
B0 = 1  
V0 = 1  
d = 1  
dt = 1e-6  
total_time = 50e-6  
num_steps = int(total_time / dt)

# Initialize particles
particles = [
    Particle(position=[20e-6, 0, 20e-6], velocity=[0, 25, 0], mass=1, charge=1),  # Example values
    
]

# Initialize PenningTrap
trap = PenningTrap(B0, V0, d, particles)

# Storage for time evolution data
positions = [[] for _ in particles]
velocities = [[] for _ in particles]

# Simulation loop
for _ in range(num_steps):
    trap.update(dt)
    for i, particle in enumerate(particles):
        positions[i].append(particle.position.copy())
        velocities[i].append(particle.velocity.copy())

# Convert positions and velocities to NumPy arrays for easier slicing
positions = [np.array(pos) for pos in positions]
velocities = [np.array(vel) for vel in velocities]

# Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for pos in positions:
    ax.plot3D(pos[:, 0], pos[:, 1], pos[:, 2])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
