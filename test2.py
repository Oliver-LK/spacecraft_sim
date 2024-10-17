import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Orbit:
    def __init__(self, eccentricity, semimajor_axis_km, inclination_deg, arg_periapsis_deg, long_asc_node_deg):
        self.eccentricity = eccentricity
        self.semimajor_axis_km = semimajor_axis_km
        self.inclination_deg = inclination_deg
        self.arg_periapsis_deg = arg_periapsis_deg
        self.long_asc_node_deg = long_asc_node_deg

    def calculate_position(self, true_anomaly_deg):
        # Convert angles to radians
        true_anomaly = np.radians(true_anomaly_deg)
        inclination = np.radians(self.inclination_deg)
        arg_periapsis = np.radians(self.arg_periapsis_deg)
        long_asc_node = np.radians(self.long_asc_node_deg)
        
        # Calculate distance from central body
        r = self.semimajor_axis_km * (1 - self.eccentricity**2) / (1 + self.eccentricity * np.cos(true_anomaly))
        
        # Position in the orbital plane
        x_orb = r * np.cos(true_anomaly)
        y_orb = r * np.sin(true_anomaly)
        z_orb = 0

        # Rotation matrices to account for inclination, argument of periapsis, and ascending node
        # Rotation around z-axis by longitude of ascending node
        x_rot_1 = x_orb * np.cos(long_asc_node) - y_orb * np.sin(long_asc_node)
        y_rot_1 = x_orb * np.sin(long_asc_node) + y_orb * np.cos(long_asc_node)
        z_rot_1 = z_orb

        # Rotation around x-axis by inclination
        x_rot_2 = x_rot_1
        y_rot_2 = y_rot_1 * np.cos(inclination) - z_rot_1 * np.sin(inclination)
        z_rot_2 = y_rot_1 * np.sin(inclination) + z_rot_1 * np.cos(inclination)

        # Rotation around z-axis by argument of periapsis
        x_final = x_rot_2 * np.cos(arg_periapsis) - y_rot_2 * np.sin(arg_periapsis)
        y_final = x_rot_2 * np.sin(arg_periapsis) + y_rot_2 * np.cos(arg_periapsis)
        z_final = z_rot_2

        return x_final, y_final, z_final

    def plot_sun(self, ax):
        # Plot the Sun as a point at the origin
        ax.scatter(0, 0, 0, color='yellow', s=100, label='Central Body (Sun)')  # Size can be adjusted with 's'

    def plot_orbit(self, ax):
        # Generate points for true anomaly from 0 to 360 degrees
        true_anomalies = np.linspace(0, 360, 500)
        x_vals = []
        y_vals = []
        z_vals = []

        for ta in true_anomalies:
            x, y, z = self.calculate_position(ta)
            x_vals.append(x)
            y_vals.append(y)
            z_vals.append(z)

        # Plotting the orbit path
        ax.plot(x_vals, y_vals, z_vals, label=f'Orbit (a={self.semimajor_axis_km} km, e={self.eccentricity})')

# Create a plot for multiple orbits
def plot_multiple_orbits(orbit_list):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the Sun
    if orbit_list:
        orbit_list[0].plot_sun(ax)  # Plot the Sun only once

    # Plot each orbit
    for orbit in orbit_list:
        orbit.plot_orbit(ax)

    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title('3D Orbit Visualization')
    ax.legend()
    plt.show()

# Example usage
orbit1 = Orbit(eccentricity=0.5, semimajor_axis_km=10000, inclination_deg=30, arg_periapsis_deg=45, long_asc_node_deg=60)
orbit2 = Orbit(eccentricity=0.3, semimajor_axis_km=12000, inclination_deg=45, arg_periapsis_deg=30, long_asc_node_deg=90)
orbit3 = Orbit(eccentricity=0.1, semimajor_axis_km=8000, inclination_deg=15, arg_periapsis_deg=60, long_asc_node_deg=150)

# Plot multiple orbits
plot_multiple_orbits([orbit1, orbit2, orbit3])
