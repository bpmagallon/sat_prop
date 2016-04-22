#julian date
from math import copysign

def julian(y, m, d, ut):

    jd = (367*y)-((7*(y+((m+9)/12)))/4)+((275*m)/9)+d+1721013.5+ut/24-copysign(0.5,100*y+m-190002.5)+0.5

    return jd
