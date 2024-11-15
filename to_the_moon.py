# ========================================================================
#   Sending a spacecraft to the moon
#   Date: 3/10/2024
#   About:  This is the main dictator file
# 
#   Author: Oliver Clements
#           
# ========================================================================

# Library Imports
import numpy as np
import copy

# Module Imports
from orbitital import *
from constatnts import LUNA_ORBIT_R_KM, LUNA_SOI_KM
from fundemental import subtract_vector_pos, subtract_vector_vel

# =============================== SETUP=========================================
# Solar bodies
Earth = give_Earth()
Moon = give_Moon()

# Variables and lists we will need
orbital_list = []
time_dict = dict()      # this dict will have the following format {Maneuver number} : {[time_h, delta_v_kmps]}
                        # Time will taken after 0 which is at simulations beginnings

# Moon orbit
MoonOrbit = Orbit(Earth)
MoonOrbit = assign_moon_orbit(MoonOrbit)


# spacecraft initial orbit
InitialOrbit = Orbit(Earth)
InitialOrbit = assign_initial_conditions(InitialOrbit)
InitialOrbit.label = "Initial Orbit"

# =============================== MANUEVER 1 =========================================
# Spacecraft after Manuever 1. This is changing the inclination of the orbit
# New Orbit elements
t_cur = 0
maneuver1_anomaly_deg   = 180
maneuver1_inc_deg       = 29
t1 = seconds2hours(InitialOrbit.calc_time_to_reach_true_anomaly(maneuver1_anomaly_deg))
InitialOrbit.true_anomaly_deg = maneuver1_anomaly_deg

# Assign maneuver 1 to the dict
maneuver1_dv = InitialOrbit.required_dv_for_inc_change(maneuver1_inc_deg)
time_dict["Manuever1"] = [float(t1), maneuver1_dv]

# New orbit 
Orbit1 = copy.deepcopy(InitialOrbit)
Orbit1.inclination_deg = maneuver1_inc_deg
Orbit1.label = "Orbit with inclination manoeuver"

# Update moon
MoonOrbit.update_true_anomaly(hours2seconds(t1))




# =============================== MANUEVER 2 =========================================
# This manoeuver will circularize the orbit. Will also occur at apoaspsis
# New orbit Elements
t2 = t1
maneuver2_eccentricity = 0
maneuver2_semimajor_axis_km = Orbit1.calculate_apoapsis()

# Assigning to time_dict
maneuver2_dv = Orbit1.required_dv_for_tangential_manoeuver(Orbit1.calculate_apoapsis())
time_dict["Manuever2"] = [float(t1), maneuver2_dv]

# New Orbit
Orbit2 = copy.deepcopy(Orbit1)
Orbit2.eccentricity = maneuver2_eccentricity
Orbit2.semimajor_axis_km = maneuver2_semimajor_axis_km
Orbit2.label = "Circular Parking Orbit"



# =============================== MANUEVER 2 =========================================
# Hohmann transfer to reach moons orbit
# New orbital elements
dt3 = hours2seconds(0.9568)
t3 = t2 + dt3
manoeuver3_semimajor_km = Orbit2.calculate_new_semi_major_axis(LUNA_ORBIT_R_KM)
manoeuver3_eccentricity = Orbit2.calculate_new_eccentricity_with_perapsis(LUNA_ORBIT_R_KM)

maneuver3_dv = Orbit2.required_dv_for_tangential_manoeuver(LUNA_ORBIT_R_KM)
time_dict["Manuever3"] = [float(seconds2hours(t3)), abs(maneuver3_dv)]

# New Orbit
Orbit3 = copy.deepcopy(Orbit2)
# Switch periapsis and apoapsis
Orbit3.long_asc_node_deg = Orbit3.long_asc_node_deg + 180   # Flip where periapsis and apoapsis are from the initial orbit
Orbit3.true_anomaly_deg = 0                                 # Since apoapsis and periapsis have switched

new_anomaly = Orbit3.calculate_new_true_anomaly_for_t(dt3)
Orbit3.long_asc_node_deg = Orbit3.long_asc_node_deg + new_anomaly
Orbit3.eccentricity = manoeuver3_eccentricity
Orbit3.semimajor_axis_km = manoeuver3_semimajor_km
Orbit3.label = "Transfer to Moon"

# We dont update the Orbit3 anomaly as it is already taken care of when we changed its long_asc_node_deg
MoonOrbit.update_true_anomaly(dt3)
Orbit2.update_true_anomaly(dt3)


# =============================== TIME DELAY =========================================
# Now we wait until the spacecraft has reached the moon
# dt = 0
dt4 = hours2seconds(0.1)
dt_accum = 0
distance_away_from_moon = 1e17
while distance_away_from_moon > LUNA_SOI_KM:
    dt_accum += dt4
    MoonOrbit.update_true_anomaly(dt4)
    Orbit3.update_true_anomaly(dt4)

    moon_distance = MoonOrbit.calculate_position(MoonOrbit.true_anomaly_deg)
    craft_distance = Orbit3.calculate_position(Orbit3.true_anomaly_deg)
    pos = subtract_vector_pos(moon_distance, craft_distance)
    distance_away_from_moon = pos.calculate_position_magnitude()



# ============================= TRANSFER INTO MOON SOI ====================================
# Difference in position to give position relative to the moon
moon_distance = MoonOrbit.calculate_position(MoonOrbit.true_anomaly_deg)
craft_distance = Orbit3.calculate_position(Orbit3.true_anomaly_deg)
pos_diff_class = subtract_vector_pos(craft_distance, moon_distance)
pos_diff = pos_diff_class.give_position_vector()

# Difference in velocity to give velocity relative to teh moon
moon_vel_vec_cartesian = MoonOrbit.convert_from_perifocal_to_cartesian_v()
craft_vel_vec_cartesian = Orbit3.convert_from_perifocal_to_cartesian_v()
vel_diff = craft_vel_vec_cartesian - moon_vel_vec_cartesian

Orbit4 = copy.deepcopy(Orbit3)
Orbit4.parent_body = Moon
Orbit4.reassign_orbital_elements_from_cartesian(pos_diff, vel_diff)
Orbit4.label = "Luna SOI Orbit"






# =============================== MANUEVER 4 =========================================
# Burn to make into moon orbit
maneuver4_anomaly_deg   = 0
maneuver4_semi_major_axis_desired = 14069
maneuver4_eccentricity_axis_desired = 0

dt5 = Orbit4.calc_time_to_reach_true_anomaly(maneuver4_anomaly_deg)
t5 = dt5 + dt_accum
Orbit4.true_anomaly_deg = maneuver4_anomaly_deg


# Manuever
maneuver4_dv = Orbit4.required_dv_for_coplanar_transfer(maneuver4_semi_major_axis_desired, maneuver4_eccentricity_axis_desired)
time_dict["Manuever4"] = [float(t5), abs(maneuver4_dv)]

maneuver4_semi_major_axis = Orbit4.calculate_new_semi_major_axis(maneuver4_semi_major_axis_desired)

# New orbit
Orbit5 = copy.deepcopy(Orbit4)
Orbit5.semimajor_axis_km = maneuver4_semi_major_axis
Orbit5.true_anomaly_deg = 180   # As now periapsis and apoapsis have switched
# Orbit5.eccentricity = Orbit5.calculate_new_eccentricity_with_perapsis(Orbit5.calculate_periapsis())
Orbit5.eccentricity = 0.54537
Orbit5.long_asc_node_deg += 180
Orbit5.label = "Luna Stable Orbit"





# print(Orbit5.semimajor_axis_km)
# print(Orbit5.inclination_deg)


# =============================== MANUEVER 5 =========================================
# Burn to make change the inclination to 90
t6 = t5
maneuver5_anomaly_deg   = 180
maneuver5_inc_deg       = 90

# Assign maneuver 1 to the dict
maneuver5_dv = Orbit5.required_dv_for_inc_change(maneuver5_inc_deg)
time_dict["Manuever5"] = [float(t5), maneuver5_dv]

# New orbit 
Orbit6 = copy.deepcopy(Orbit5)
Orbit6.inclination_deg = maneuver5_inc_deg
Orbit6.label = "Polar Moon Orbit"

print(Orbit6.semimajor_axis_km)
print(Orbit6.eccentricity)
print(Orbit6.inclination_deg)
print(Orbit6.long_asc_node_deg)
print(Orbit6.arg_periapsis_deg)
# print(t5/3600)
# print(Orbit5.calculate_velocity(Orbit5.true_anomaly_deg).calculate_velocity_magnitude())

# print(time_dict)
# print(Orbit6.semimajor_axis_km)
# print(seconds2hours(Orbit6.calculate_orbital_period()))


orbital_list = [Orbit4, Orbit5]
# Orbit5.plot_orbits(orbital_list)

# print(time_dict)
print(time_dict)


total_dv = 0
for value in time_dict.values():
    total_dv += abs(value[1] )

print(total_dv)
# print(seconds2hours(time_dict["Manuever5"][0]))
