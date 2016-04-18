#groups the tle values
import re
from math import radians, pi

def group(arr):
    tle1 = arr[1]
    tle2 = arr[2]
    tle1 = re.sub(r"\s{1,}", ",", tle1)
    tle1_split = tle1.split(",")

    #epoch
    epoch = tle1_split[3]
    
    #drag value
    bstar = tle1_split[6] 

    if bstar[0]=='+':
        bstar = float("0."+bstar[1:6])/(10**float(bstar[7]))
    else:
        bstar = float("0."+bstar[1:6])/(10**float(bstar[7]))*-1
        
    tle2 = re.sub(r"\s{1,}", ",", tle2)
    tle2_split = tle2.split(",")

    #mean motion
    n = tle2_split[7] 

    n = float(n[:11])
    n = n*((2*pi)/(1440))

    #inclination
    i = radians(float(tle2_split[2]))
    
    #raan
    node = float(tle2_split[3])
    
    #eccentricity
    e = float("0."+tle2_split[4])

    #argument of perigee
    omega = radians(float(tle2_split[5]))

    #mean anomaly
    m = radians(float(tle2_split[6]))
    
    return epoch, n, i, node, e, omega, m, bstar
