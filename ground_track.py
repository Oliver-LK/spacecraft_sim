import numpy as np
import matplotlib.pyplot as plt

# Constants for the Moon
MU_MOON = 4902.800066  # Moon's gravitational parameter, km^3/s^2
MOON_RADIUS = 1737.4  # Moon's radius, km
MOON_ROTATION_PERIOD = 27.32 * 24 * 3600  # Moon's sidereal rotation period in seconds

def orbital_elements_to_position(a, e, i, omega, raan, M):
    """Calculate the satellite position in Moon-centered inertial (MCI) coordinates."""
    # Convert angles from degrees to radians
    i = np.deg2rad(i)
    omega = np.deg2rad(omega)
    raan = np.deg2rad(raan)

    # Solve Kepler's Equation to get the eccentric anomaly (E)
    E = M  # Initial guess for eccentric anomaly
    for _ in range(10):  # Iterate to refine the solution
        E = M + e * np.sin(E)

    # True anomaly (theta)
    theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

    # Radius in the orbital plane
    r = a * (1 - e * np.cos(E))

    # Position in the orbital plane (Perifocal coordinate system)
    x_orb = r * np.cos(theta)
    y_orb = r * np.sin(theta)
    z_orb = 0  # In the orbital plane, z-coordinate is zero

    # Rotation matrices to convert from the orbital plane to MCI coordinates
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
    r_mci = rotation_matrix @ np.array([x_orb, y_orb, z_orb])

    return r_mci

def plot_ground_track(a, e, i, omega, raan, M0, num_orbits=1, num_points_per_orbit=500):
    # Setup a blank plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_aspect('auto')

    # Set up a blank grid (no surface, just a white background)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")

    lons, lats = [], []

    # Calculate the orbital period
    T = 2 * np.pi * np.sqrt(a**3 / MU_MOON)  # Orbital period in seconds

    # Sweep through multiple orbits by adjusting the range of the mean anomaly
    for orbit in range(num_orbits):
        for idx, M in enumerate(np.linspace(0, 2 * np.pi, num_points_per_orbit)):
            # Calculate elapsed time since start of the first orbit
            elapsed_time = (M + orbit * 2 * np.pi) / (2 * np.pi) * T

            # Calculate the satellite position in MCI coordinates
            r_mci = orbital_elements_to_position(a, e, i, omega, raan, M + orbit * 2 * np.pi)
            lon = np.rad2deg(np.arctan2(r_mci[1], r_mci[0]))
            lat = np.rad2deg(np.arcsin(r_mci[2] / np.linalg.norm(r_mci)))

            # Adjust longitude for the Moon's rotation over elapsed time
            moon_rotation = (360 * elapsed_time / MOON_ROTATION_PERIOD) % 360
            lon = (lon - moon_rotation + 360) % 360  # Adjust for Moon's rotation

            # Convert longitude to range [-180, 180]
            if lon > 180:
                lon -= 360

            lons.append(lon)
            lats.append(lat)

    # Plot the ground track on the plain background
    ax.plot(lons, lats, marker='.', color='red', markersize=2, linestyle='-')
    ax.set_title(f'Satellite Ground Track over the Moon ({num_orbits} Orbits)')
    plt.grid(True)
    plt.show()

# Example orbital parameters to test the function (for a lunar orbit)
a = 30579  # Semi-major axis in km
e = 0.545  # Eccentricity
i = 90  # Inclination in degrees
omega = 104.6  # Argument of perigee in degrees
raan = 205.9  # Right ascension of the ascending node in degrees
M0 = 0  # Mean anomaly at epoch in degrees

# Plotting the ground track for a specified number of orbits
plot_ground_track(a, e, i, omega, raan, M0, num_orbits=5)
