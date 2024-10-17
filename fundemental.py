# ========================================================================
#   Sending a spacecraft to the moon
#   Date: 3/10/2024
#   About:  Velocity classes
# 
#   Author: Oliver Clements
#           
# ========================================================================

from __future__ import annotations  # Needed for forward references
import copy
import math
import numpy as np

def calculate_r_scalar_km(semimajor_axis_km: float, eccentricity: float, true_anomaly_deg: float) -> float:
    """ Calculates the scalar for the position vector from orbital body"""
    return semimajor_axis_km * (1 - eccentricity**2) / (1 + eccentricity * math.cos(np.deg2rad(true_anomaly_deg)))


def calculate_r_vector_km(semimajor_axis_km: float, eccentricity: float, true_anomaly_deg: float) -> np.ndarray:
    """ Calculates the r vector"""
    r_vec_km    = np.zeros(3)
    r_km = calculate_r_scalar_km(semimajor_axis_km, eccentricity, true_anomaly_deg)
    r_vec_km[0] = r_km * math.cos(math.radians(true_anomaly_deg))
    r_vec_km[1] = r_km * math.sin(math.radians(true_anomaly_deg))
    r_vec_km[2] = 0

    return r_vec_km


def calculate_r_hat_km(emimajor_axis_km: float, eccentricity: float, true_anomaly_deg: float) -> np.ndarray:
    """ Finds the rhat of the orbit"""

def calculate_v_hat_km(emimajor_axis_km: float, eccentricity: float, true_anomaly_deg: float) -> np.ndarray:
    """ Finds the vhat of the orbit"""

def calculate_h_hat_km(emimajor_axis_km: float, eccentricity: float, true_anomaly_deg: float) -> np.ndarray:
    """ Finds the hhat of the orbit"""


class Position:
    """ Class that handles the position vector"""

    def __init__(self) -> None:
        self.x_km   = None
        self.y_km   = None
        self.z_km   = None

    def __str__(self) -> str:
        return f"x: {self.x_km:.3f} km\ty: {self.y_km:.3f} km\tz: {self.z_km:.3f} km"


    def calculate_position(self, true_anomaly_deg: float, inclination_deg: float, arg_periapsis_deg, long_asc_node_deg, r_km: float)  -> None:
        # Convert angles to radians
        true_anomaly = np.radians(true_anomaly_deg)
        inclination = np.radians(inclination_deg)
        arg_periapsis = np.radians(arg_periapsis_deg)
        long_asc_node = np.radians(long_asc_node_deg)
        
        # Position in the orbital plane
        x_orb = r_km * np.cos(true_anomaly)
        y_orb = r_km * np.sin(true_anomaly)
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
        self.x_km = x_rot_2 * np.cos(arg_periapsis) - y_rot_2 * np.sin(arg_periapsis)
        self.y_km = x_rot_2 * np.sin(arg_periapsis) + y_rot_2 * np.cos(arg_periapsis)
        self.z_km = z_rot_2

        return self
    

    def calculate_position_magnitude(self):
        """ Calculates the magnitude of the position vector"""
        return math.sqrt(self.x_km ** 2 + self.y_km ** 2 + self.z_km ** 2)
    
    def give_position_vector(self) -> np.ndarray:
        new_position_vector_kmps = np.array([
            self.x_km,
            self.y_km,
            self.z_km
            ])
        
        return new_position_vector_kmps


class Velocity:
    """ Class handles all the velocity related tasks"""

    def __init__(self) -> None:
        self.v_r_kmps       = None
        self.v_theta_kmps   = None


    def __str__(self) -> str:
        return f"v_r: {self.v_r_kmps:.3f} km/s,\tv_theta: {self.v_theta_kmps:.3f} km/s"

    def calculate_velocity_vector(self, anomaly_deg: float, mue_kmps: float, eccentricity: float, angular_momentum: float, radial_distance_km) -> None:
        """ Calculates each component for the velocity vector"""
        self.v_r_kmps           = mue_kmps / angular_momentum * eccentricity * math.sin(math.radians(anomaly_deg))
        self.v_theta_kmps       = angular_momentum / radial_distance_km

        return self
    
    def calculate_velocity_magnitude(self) -> float:
        """ Calculates the scalar magnitude of the velocity vector"""
        return math.sqrt(self.v_r_kmps ** 2 + self.v_theta_kmps ** 2)
    

    def give_velocity_vector(self) -> np.ndarray:
        new_velocity_vector_kmps = np.array([
            self.v_r_kmps,
            self.v_theta_kmps,
            ])
        return new_velocity_vector_kmps
    


def newton_method(f, x0, tol=1e-8, max_iter=100):
    """ Newton-Raphson method to solve for the root of an equation f(x) = 0 """
    x = x0
    for _ in range(max_iter):
        fx = f(x)
        dfx = (f(x + tol) - fx) / tol  # Numerical derivative
        if abs(fx) < tol:
            return x
        x -= fx / dfx
    raise ValueError("Newton's method did not converge.")


def subtract_vector_pos(pos1: Position, pos2: Position) -> Position:
    """ Subtracts 2 position vectors"""
    pos = Position()
    pos.x_km = pos1.x_km - pos2.x_km
    pos.y_km = pos1.y_km - pos2.y_km
    pos.z_km = pos1.z_km - pos2.z_km

    return pos

    
def subtract_vector_vel(vel1: Velocity, vel2: Velocity) -> Velocity:
    """ Subtracts 2 position vectors"""
    vel = Velocity()
    vel.v_r_kmps = vel1.v_r_kmps - vel2.v_r_kmps
    vel.v_theta_kmps = vel1.v_theta_kmps - vel2.v_theta_kmps

    return vel




def get_h_vector(self, pos: Position, vel: Velocity) -> np.ndarray:
    """ Calculates the h vector"""
    vel_vec = vel.give_velocity_vector()