# visualize_data.py

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data from file
positions = np.load('positions.npy')
velocities = np.load('velocities.npy')

# Constants for plotting
total_time = 50e-6
num_steps = positions.shape[1]  # assuming the second dimension is the time step
time = np.linspace(0, total_time, num_steps)

# 1. Plot of the motion in the z direction as a function of time
plt.figure()
for i, z in enumerate(positions[:, :, 2]):  # assuming positions are structured as [particle][time_step][coordinate]
    plt.plot(time, z, label=f'Particle {i+1}')
plt.xlabel('Time (s)')
plt.ylabel('Z position (m)')
plt.title('Z position over time')   
plt.legend()
plt.show()


# For 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for pos in positions:
    ax.plot(pos[:, 0], pos[:, 1], pos[:, 2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('3D trajectory of particles')
plt.show()
