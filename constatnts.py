# ========================================================================
#   Sending a spacecraft to the moon
#   Date: 3/10/2024
#   About:  Global Constants
# 
#   Author: Oliver Clements
#           
# ========================================================================


EARTH_G_C_KMS   = 398600 
EARTH_RADIUS_KM = 6371

LUNA_G_C_KMS    = 4905
LUNA_RADIUS     = 1737
LUNA_SOI_KM     = 66200
LUNA_ORBIT_R_KM = 384400
LUNA_INC_DEG    = 28.58
LUNA_RAAN       = 90
LUNA_ARG_PERI   = 270


def seconds2hours(seconds: float) -> float:
    """ Converts seconds to hours"""
    return seconds / 3600


def hours2seconds(hours: float) -> float:
    """ Converts hours to seconds"""
    return hours * 3600
