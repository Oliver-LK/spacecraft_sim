# ========================================================================
#   Sending a spacecraft to the moon
#   Date: 3/10/2024
#   About:  Handles the orbital parameters and plotting
# 
#   Author: Oliver Clements
#           
# ========================================================================

# Library Imports 
from __future__ import annotations  # Needed for forward references
import math
from typing import Union
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Module imports
from initial_conditions import *
from constatnts import *
from fundemental import Velocity, Position, calculate_r_scalar_km, newton_method

class ParentBody:
    """ Class that contains information about the parent body"""

    def __init__(self) -> None:
        self.radius_km  = None
        self.mue_kmps2  = None
        self.soi_km     = None
        self.name       = None

class Orbit:
    """ Class containing orbital elements"""
    # ============================== Orbital elements and related calculations ===================================

    def __init__(self, parent_body: ParentBody) -> None:
        # Orbit parameters
        self.eccentricity       = None
        self.semimajor_axis_km  = None
        self.inclination_deg    = None
        self.arg_periapsis_deg  = None
        self.long_asc_node_deg  = None
        self.true_anomaly_deg   = None
        self.parent_body        = parent_body
        self.angular_momentum   = None

        # Stuff for Plotting
        self.label              = None
        self.max_range          = None


    def __str__(self) -> str:
        return f"e: {self.eccentricity}, a: {self.semimajor_axis_km}km, i: {self.inclination_deg}deg, h: {self.angular_momentum}"


    def calculate_position(self, true_anomaly_deg) -> Position:
        """ Calculates the position vector for a given anomaly"""
        pos_km  = Position()
        r_km    = calculate_r_scalar_km(self.semimajor_axis_km, self.eccentricity, true_anomaly_deg)
        pos_km.calculate_position(true_anomaly_deg, self.inclination_deg, self.arg_periapsis_deg, self.long_asc_node_deg, r_km)

        return pos_km
    

    def calculate_velocity(self, true_anomaly_deg) -> Velocity:
        """ Calculates the velocity vector for a given anomaly"""
        vel_km  = Velocity()
        r_km    = calculate_r_scalar_km(self.semimajor_axis_km, self.eccentricity, true_anomaly_deg)
        vel_km.calculate_velocity_vector(true_anomaly_deg, self.parent_body.mue_kmps2, self.eccentricity, self.angular_momentum, r_km)

        return vel_km
    

    def calculate_orbital_period(self) -> float:
        """ Finds the period of the orbit"""
        return 2 * math.pi * math.sqrt(self.semimajor_axis_km ** 3 / self.parent_body.mue_kmps2)
    

    def calculate_mean_motion(self) -> float:
        """ Calculates the mean motion of an object"""
        return math.sqrt(self.parent_body.mue_kmps2 / abs(self.semimajor_axis_km ** 3))
    

    def calculate_eccentric_anomaly(self) -> float:
        """ Calculates the eccentric anomaly"""
        true_anomaly_radians = np.deg2rad(self.true_anomaly_deg)    # Need as radians for calculations
        eccentric_anomaly = 2 * np.arctan(np.sqrt(abs(1 - self.eccentricity) / (1 + self.eccentricity)) * np.tan(true_anomaly_radians / 2))
        return eccentric_anomaly
    

    def calculate_mean_anomaly(self, eccentric_anomaly) -> float:
        """ Calculates the mean anomaly"""
        mean_anomaly = eccentric_anomaly - self.eccentricity * np.sin(eccentric_anomaly)
        return mean_anomaly
    

    def calculate_flight_path_angle_deg(self) -> float:
        """ Calculates the flight path angle"""
        true_anomaly_rad = np.deg2rad(self.true_anomaly_deg)
        flightpath_angle_rad = self.eccentricity * math.sin(true_anomaly_rad) / (1 - self.eccentricity * math.cos(true_anomaly_rad))

        return np.rad2deg(flightpath_angle_rad)
    

    def calculate_apoapsis(self) -> float:
        """ Finds the apoaspsis on a orbit"""
        return self.semimajor_axis_km * (1 + self.eccentricity)
    

    def calculate_periapsis(self) -> float:
        """ Finds the periapsis on a orbit"""
        return self.semimajor_axis_km * (1 - self.eccentricity)
    
    def calculate_new_semi_major_axis(self, new_r_km: float) -> float:
        """ Finds the new semimajor axis based on new orbit"""
        r_km = calculate_r_scalar_km(self.semimajor_axis_km, self.eccentricity, self.true_anomaly_deg)
        return (r_km + new_r_km) / 2
    

    def calculate_new_eccentricity_with_perapsis(self, new_r_km: float) -> float:
        """ Updates the eccentricity based on a new apoapsis and same periapsis"""
        return 1 - self.calculate_periapsis() / self.calculate_new_semi_major_axis(new_r_km)
    

    def calculate_new_eccentricity_with_apoapsis(self, new_r_km: float) -> float:
        """ Updates the eccentricity based on a new apoapsis and same periapsis"""
        return  self.calculate_apoapsis() / self.calculate_new_semi_major_axis(new_r_km) -1

    
    


    
    # ============================== Forward in time calculations ===================================

    def update_true_anomaly(self, delta_t_s: float) -> None:
        """ Updates the true anomaly for a given time (s)"""

        n = self.calculate_mean_motion()    # Mean motion
        true_anomaly_radians = np.deg2rad(self.true_anomaly_deg)    # Need as radians for calculations

        # Anomaly parameters for orbit
        E0 = self.calculate_eccentric_anomaly()
        M0 = self.calculate_mean_anomaly(E0)

        # Calculate the new mean anomaly (M_new)
        M_new = M0 + n * delta_t_s

        # Solve Kepler's equation for eccentric anomaly (E) using Newton's method
        def kepler_equation(E):
            return E - self.eccentricity * np.sin(E) - M_new

        E_new = newton_method(kepler_equation, E0)  # Initial guess is E0

        # Step 5: Calculate the new true anomaly (theta) from the new eccentric anomaly (E_new)
        true_anomaly_rad = 2 * np.arctan(np.sqrt((1 + self.eccentricity) / (1 - self.eccentricity)) * np.tan(E_new / 2))
        true_anomaly_deg = np.degrees(true_anomaly_rad)

        self.true_anomaly_deg = true_anomaly_deg

    def calculate_new_true_anomaly_for_t(self, delta_t_s: float) -> None:
        """ Updates the true anomaly for a given time (s)"""

        n = self.calculate_mean_motion()    # Mean motion
        true_anomaly_radians = np.deg2rad(self.true_anomaly_deg)    # Need as radians for calculations

        # Anomaly parameters for orbit
        E0 = self.calculate_eccentric_anomaly()
        M0 = self.calculate_mean_anomaly(E0)

        # Calculate the new mean anomaly (M_new)
        M_new = M0 + n * delta_t_s

        # Solve Kepler's equation for eccentric anomaly (E) using Newton's method
        def kepler_equation(E):
            return E - self.eccentricity * np.sin(E) - M_new

        E_new = newton_method(kepler_equation, E0)  # Initial guess is E0

        # Step 5: Calculate the new true anomaly (theta) from the new eccentric anomaly (E_new)
        true_anomaly_rad = 2 * np.arctan(np.sqrt((1 + self.eccentricity) / (1 - self.eccentricity)) * np.tan(E_new / 2))
        true_anomaly_deg = np.degrees(true_anomaly_rad)

        return true_anomaly_deg


    def calc_time_to_reach_true_anomaly(self, target_true_anomaly_deg: float) -> float:
        """ Calculates the time to reach a specific anomaly"""
        # Step 1: Calculate the mean motion (n)
        mean_motion = self.calculate_mean_motion()

        # Step 2: Convert target true anomaly to eccentric anomaly
        target_true_anomaly_rad = np.radians(target_true_anomaly_deg)
        eccentric_anomaly_target = 2 * np.arctan(
            np.sqrt(abs(1 - self.eccentricity) / (1 + self.eccentricity)) * np.tan(target_true_anomaly_rad / 2)
        )

        # Calculate the current eccentric anomaly (E0) from the current true anomaly
        E0 = self.calculate_eccentric_anomaly()

        # Calculate the mean anomalies (current and target)
        M0 = self.calculate_mean_anomaly(E0)
        M1 = eccentric_anomaly_target - self.eccentricity * np.sin(eccentric_anomaly_target)

        # Calculate the time difference between the target and current mean anomaly
        dM = M1 - M0

        # Ensure delta_mean_anomaly is within the range [0, 2π] to get the positive time difference
        dM = (dM + 2 * np.pi) % (2 * np.pi)

        # Calculate the time to reach the target anomaly
        time_to_target = dM / mean_motion  # Time in seconds

        return time_to_target
    

    # ============================== Manuever calculations ===================================

    def elliptical_v_at_anomaly(self, true_anomaly_deg: float) -> float:
        """ Finds the velocity at some anomaly"""
        r_km = calculate_r_scalar_km(self.semimajor_axis_km, self.eccentricity, true_anomaly_deg)
        return math.sqrt((2 * self.parent_body.mue_kmps2 / r_km) - (self.parent_body.mue_kmps2 / self.semimajor_axis_km))
    

    def elliptical_v_at_anomaly_for_new_r(self, new_r_km: float, new_semimajor_km: float) -> float:
        """ Finds the velocity at some anomaly"""
        
        return math.sqrt((2 * self.parent_body.mue_kmps2 / new_r_km) - (self.parent_body.mue_kmps2 / new_semimajor_km))
    

    def elliptical_v_at_anomaly_new_semi(self, true_anomaly_deg: float, a_t: float) -> float:
        """ Finds the velocity at some anomaly"""
        r_km = calculate_r_scalar_km(self.semimajor_axis_km, self.eccentricity, true_anomaly_deg)
        return math.sqrt((2 * self.parent_body.mue_kmps2 / r_km) - (self.parent_body.mue_kmps2 / a_t))
    

    def required_dv_for_tangential_manoeuver(self, new_r_km) -> float:
        """ Finds the required dv for a single tangential manuever. 
            Assumes manoeuver takes place at current anomaly 
        """
        new_semimajor_km = self.calculate_new_semi_major_axis(new_r_km)
        v1_kmps = self.elliptical_v_at_anomaly(self.true_anomaly_deg)
        v2_kmps = self.elliptical_v_at_anomaly_for_new_r(new_r_km, new_semimajor_km)

        dv = v2_kmps - v1_kmps

        return dv


    def required_dv_for_coplanar_transfer(self, a2_km: float, e2: float) -> float:
        """ Calculates the dv required for a coplanar transfer"""
        rp1_km = self.semimajor_axis_km * (1 - self.eccentricity)
        rp2_km = a2_km * (1 - e2)
        at = self.calculate_new_semi_major_axis(rp2_km)

        v1 = self.elliptical_v_at_anomaly(self.true_anomaly_deg)
        v2 = self.elliptical_v_at_anomaly_new_semi(self.true_anomaly_deg, at)

        dv = v2 - v1

        return dv


    def required_dv_for_inc_change(self, desired_inclination_deg) -> float:
        """ Calculates the velocity for a change in inclination"""
        di_rad = np.deg2rad(self.inclination_deg - desired_inclination_deg)
        
        velocity_kmps = self.elliptical_v_at_anomaly(self.true_anomaly_deg)

        flight_path_angle_deg = self.calculate_flight_path_angle_deg()
        dv_kmps = 2 * velocity_kmps * math.cos(np.deg2rad(flight_path_angle_deg)) * math.sin(di_rad/2)

        return dv_kmps
    

# ============================== Coordinate Systems ===================================
    def convert_from_perifocal_to_cartesian(self, vec: np.ndarray) -> np.ndarray:
        """ Converts some vector in the perifocal reference frame to the cartesian reference frame"""

        # Step 3: Rotation matrix from the perifocal to the Cartesian coordinate system
        cos_Omega = np.cos(np.deg2rad(self.long_asc_node_deg))
        sin_Omega = np.sin(np.deg2rad(self.long_asc_node_deg))
        cos_omega = np.cos(np.deg2rad(self.arg_periapsis_deg))
        sin_omega = np.sin(np.deg2rad(self.arg_periapsis_deg))
        cos_i = np.cos(np.deg2rad(self.inclination_deg))
        sin_i = np.sin(np.deg2rad(self.inclination_deg))


        rotation_matrix = np.array([
            [cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i, -cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i, sin_Omega * sin_i],
            [sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i, -sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i, -cos_Omega * sin_i],
            [sin_omega * sin_i, cos_omega * sin_i, cos_i]
        ])

        return rotation_matrix @ vec
    

    def convert_from_perifocal_to_cartesian_v(self) -> np.ndarray:
        """ Converts velocity vector in the perifocal reference frame to the cartesian reference frame"""
        self.angular_momentum = math.sqrt(self.parent_body.mue_kmps2 * self.semimajor_axis_km * (1 - self.eccentricity ** 2))
        v_p_kmps = - self.parent_body.mue_kmps2 / self.angular_momentum * math.sin(np.deg2rad(self.true_anomaly_deg))
        v_q_kmps = self.parent_body.mue_kmps2 / self.angular_momentum * (self.eccentricity + math.cos(np.deg2rad(self.true_anomaly_deg)))

        vel_vec = np.array([v_p_kmps, v_q_kmps, 0]).T

        return self.convert_from_perifocal_to_cartesian(vel_vec)
    

    def reassign_orbital_elements_from_cartesian(self, pos_vec: np.ndarray, vel_vec: np.ndarray) -> None:
        """ Reassign all the orbital elements from cartesian reference frame"""
        # normalising
        pos_norm = np.linalg.norm(pos_vec)
        vel_norm = np.linalg.norm(vel_vec)

        # Angular Momentum
        h_vec = np.cross(pos_vec, vel_vec)
        h_norm = np.linalg.norm(h_vec)
        self.angular_momentum = h_norm

        # Eccentricity
        e_vec = np.cross(vel_vec, h_vec) / self.parent_body.mue_kmps2 - (pos_vec / pos_norm)
        e_norm = np.linalg.norm(e_vec)
        self.eccentricity = e_norm

        # Semimajor Axis
        semi_major_axis = h_norm ** 2 / (self.parent_body.mue_kmps2 * (1 - e_norm ** 2))
        self.semimajor_axis_km = semi_major_axis

        # True Anomnaly
        true_anomaly_rad = np.arccos(np.dot(e_vec, pos_vec) / (e_norm * pos_norm))
        self.true_anomaly_deg = np.rad2deg(true_anomaly_rad)

        # Inclination
        inclination_rad = np.arccos(h_vec[2] / h_norm)
        self.inclination_deg = np.rad2deg(inclination_rad)

        # Right accession of right ascending node
        k = np.array([0, 0, 1]).T
        N = np.cross(k, h_vec)
        raan_rad = np.arccos(N[0]/np.linalg.norm(N))
        self.long_asc_node_deg = np.rad2deg(raan_rad)

        # Argument of periapsis
        arg_peri_rad = np.arccos(np.dot(N, e_vec) / (np.linalg.norm(N) * e_norm))
        self.arg_periapsis_deg = np.rad2deg(arg_peri_rad)






    
    

    # ============================== Plotting and related ===================================
    
    def plot_earth(self, ax) -> None:
        """ Create a sphere to represent the earth"""
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = self.sun_radius_km * np.outer(np.cos(u), np.sin(v))
        y = self.sun_radius_km * np.outer(np.sin(u), np.sin(v))
        z = self.sun_radius_km * np.outer(np.ones(np.size(u)), np.cos(v))

        # Plot the Sun at the origin
        ax.plot_surface(x, y, z, color='cyan', alpha=0.6)


    def plot_earth(self, ax) -> None:
        """ Plots the earth at the origin"""
        ax.scatter(0, 0, 0, color='Cyan', s=100, label='Earth (not to scale)')  # Size can be adjusted with 's'

    def plot_moon(self, ax) -> None:
        """ Plots the earth at the origin"""
        ax.scatter(0, 0, 0, color='Grey', s=100, label='Moon (not to scale)')  # Size can be adjusted with 's'


    def plot_orbit(self, ax):
        """ Defines a single instance for a plot (but does not plot) """
        # Generate points for true anomaly from 0 to 360 degrees
        true_anomalies = np.linspace(0, 360, 500)
        x_vals = []
        y_vals = []
        z_vals = []

        for ta in true_anomalies:
            pos_km = self.calculate_position(ta)
            x_vals.append(pos_km.x_km)
            y_vals.append(pos_km.y_km)
            z_vals.append(pos_km.z_km)

        # Plotting the orbit path
        ax.plot(x_vals, y_vals, z_vals, label=self.label)
        current_position = self.calculate_position(self.true_anomaly_deg)
        ax.scatter(current_position.x_km, current_position.y_km, current_position.z_km, color='black', s=20, label=f'{self.label} Current Position')

        return x_vals, y_vals, z_vals  # Return the data for axis scaling calculations


    def plot_orbits(self, orbit_list):
        """ Plots multiple orbits (does plot) for a given list """
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        max_range = 0  # Initialize the max range variable
        x_limits = [0, 0]
        y_limits = [0, 0]
        z_limits = [0, 0]

        # Plot each orbit and update axis limits
        for orbit in orbit_list:
            x_vals, y_vals, z_vals = orbit.plot_orbit(ax)
            
            # Update axis limits to ensure equal scaling
            x_limits = [min(x_limits[0], min(x_vals)), max(x_limits[1], max(x_vals))]
            y_limits = [min(y_limits[0], min(y_vals)), max(y_limits[1], max(y_vals))]
            z_limits = [min(z_limits[0], min(z_vals)), max(z_limits[1], max(z_vals))]

            max_range = max(max_range, np.ptp(x_vals), np.ptp(y_vals), np.ptp(z_vals))

        # Center the plot based on the largest range
        if max(x_limits) < 1e7 and max(y_limits) < 1e7 and max(z_limits) < 1e7:
            mid_x = (x_limits[0] + x_limits[1]) / 2
            mid_y = (y_limits[0] + y_limits[1]) / 2
            mid_z = (z_limits[0] + z_limits[1]) / 2
            ax.set_xlim(mid_x - max_range / 2, mid_x + max_range / 2)
            ax.set_ylim(mid_y - max_range / 2, mid_y + max_range / 2)
            ax.set_zlim(mid_z - max_range / 2, mid_z + max_range / 2)
        
            ax.set_xlim(mid_x - max_range / 2, mid_x + max_range / 2)
            ax.set_ylim(mid_y - max_range / 2, mid_y + max_range / 2)
            ax.set_zlim(mid_z - max_range / 2, mid_z + max_range / 2)

        else:
            ax.set_xlim(-7e4, 7e4)
            ax.set_ylim(-7e4, 7e4)
            ax.set_zlim(-7e4, 7e4)


        if self.parent_body.name == "Earth":
            orbit_list[0].plot_earth(ax)  # Plot the Earth only once

        elif self.parent_body.name == "Moon":
            orbit_list[0].plot_moon(ax)


        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
        ax.set_title('3D Orbit Visualization')
        ax.legend(loc='upper right')
        plt.show()



def give_Earth() -> ParentBody:
    """ Returns earth"""
    earth = ParentBody()
    earth.radius_km = EARTH_RADIUS_KM
    earth.mue_kmps2 = EARTH_G_C_KMS
    earth.soi_km    = 1e12
    earth.name      = "Earth"
    return earth

def give_Moon() -> ParentBody:
    """ Returns the moon"""
    moon = ParentBody()
    moon.radius_km  = LUNA_RADIUS
    moon.mue_kmps2  = LUNA_G_C_KMS
    moon.soi_km     = LUNA_SOI_KM 
    moon.name       = "Moon"
    return moon


def assign_initial_conditions(craft: Orbit) -> Orbit:
    """ Assigns teh initial conditions of the craft"""
    craft.eccentricity      = ECCENTRICITY_0
    craft.semimajor_axis_km = SEMIMAJOR_AXIS_0_KM
    craft.inclination_deg   = INCLINATION_0_DEG
    craft.arg_periapsis_deg = ARG_PERIAPSIS_0_DEG
    craft.long_asc_node_deg = LONG_ASC_NODE_0_DEG
    craft.true_anomaly_deg  = TRUE_ANOMALY_0_DEG
    craft.angular_momentum  = math.sqrt(craft.parent_body.mue_kmps2 * craft.semimajor_axis_km * (1 - craft.eccentricity ** 2))

    return craft


def assign_moon_orbit(moon: Orbit) -> Orbit:
    """ Sets the moon orbit"""
    moon.eccentricity       = 0     # as circular orbit assumption
    moon.semimajor_axis_km  = LUNA_ORBIT_R_KM
    moon.inclination_deg    = LUNA_INC_DEG
    moon.arg_periapsis_deg  = LUNA_ARG_PERI
    moon.long_asc_node_deg  = LUNA_RAAN
    moon.true_anomaly_deg   = 0
    moon.angular_momentum   =  moon.angular_momentum  = math.sqrt(moon.parent_body.mue_kmps2 * moon.semimajor_axis_km * (1 - moon.eccentricity ** 2))

    moon.label              = "Moons Orbit"

    return moon

    

def update_all_anomalies(orbit_list: Orbit, delta_t_s: float) -> None:
    """ Updates the anomalies for a list of orbits"""

    for orbit in orbit_list:
        orbit.update_true_anomaly(delta_t_s)
