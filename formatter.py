#groups the tle values
import re
from math import radians, pi

def group(arr):
    tle1 = arr[1]
    tle2 = arr[2]
    tle1 = re.sub(r"\s{1,}", ",", tle1)
    tle1_split = tle1.split(",")
    bstar = tle1_split[6]

    if bstar[0]=='+':
        bstar = float("0."+bstar[1:6])/(10**float(bstar[7]))
    
    tle2 = re.sub(r"\s{1,}", ",", tle2)
    tle2_split = tle2.split(",")
    n = tle2_split[7]

    n = float(n[:11])
    n = n*((2*pi)/(1440))
    
    i = float(tle2_split[2])
    e = float("0."+tle2_split[4])

    i = radians(i)
    omega = radians(float(tle2_split[5]))
    
    return n,i,e,bstar,omega
