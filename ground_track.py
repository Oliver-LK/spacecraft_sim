import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Constants
MU_EARTH = 398600.4418  # Earth's gravitational parameter, km^3/s^2
EARTH_RADIUS = 6378.137  # Earth's radius, km

def orbital_elements_to_position(a, e, i, omega, raan, M, num_points=500):
    """Calculate the ground track positions using given orbital elements."""
    # Convert angles from degrees to radians
    i = np.deg2rad(i)
    omega = np.deg2rad(omega)
    raan = np.deg2rad(raan)

    # Solve Kepler's Equation to get the eccentric anomaly (E)
    E = M  # Initial guess for eccentric anomaly
    for _ in range(10):  # Iterate to refine solution
        E = M + e * np.sin(E)

    # True anomaly (theta)
    theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

    # Radius in the orbital plane
    r = a * (1 - e * np.cos(E))

    # Position in the orbital plane (Perifocal coordinate system)
    x_orb = r * np.cos(theta)
    y_orb = r * np.sin(theta)
    z_orb = 0  # In the orbital plane, z-coordinate is zero

    # Rotation matrices to convert from orbital plane to ECI coordinates
    R3_W = np.array([[np.cos(-raan), np.sin(-raan), 0],
                     [-np.sin(-raan), np.cos(-raan), 0],
                     [0, 0, 1]])

    R1_i = np.array([[1, 0, 0],
                     [0, np.cos(-i), np.sin(-i)],
                     [0, -np.sin(-i), np.cos(-i)]])

    R3_w = np.array([[np.cos(-omega), np.sin(-omega), 0],
                     [-np.sin(-omega), np.cos(-omega), 0],
                     [0, 0, 1]])

    rotation_matrix = R3_W @ R1_i @ R3_w
    r_eci = rotation_matrix @ np.array([x_orb, y_orb, z_orb])

    return r_eci

# Generate the map for plotting
def plot_ground_track(orbital_parameters):
    # Setup the Basemap projection
    m = Basemap(projection='mill', lon_0=0)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 60), labels=[0, 0, 0, 1])

    lons, lats = [], []

    for params in orbital_parameters:
        r_eci = orbital_elements_to_position(*params)
        lon, lat = np.rad2deg(np.arctan2(r_eci[1], r_eci[0])), r_eci[2]
        x, y = m(lon, lat)
        lons.append(x)
        lats.append(y)

        #Plot the trajectory, simulated at finer grid.
        m.plot(lons,lats)

    plt.title('Satellite Ground Track Projection')
    plt.show()