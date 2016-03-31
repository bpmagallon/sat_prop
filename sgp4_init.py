import urllib
import formatter as fr
import constants as con
from math import cos, sin, pi

#get TLE
#celestrak = urllib.urlopen("http://www.celestrak.com/NORAD/elements/stations.txt")
celestrak = open("sample_tle.txt", "r") #for sample data
tle_data = celestrak.readlines()

#get Orbital Elements
tle_ele = fr.group(tle_data)

#compute semi major and mean motion
#initialization

###tle elements###
n = tle_ele[0]
i = tle_ele[1]
e = tle_ele[2]
bstar = tle_ele[3]
omega = tle_ele[4]
##################

a1 = (con.ke/n)**(con.exp23)

angle = cos(i)
e2 = e**2

factor = ((3*angle*angle)-1)/((1-e2)**con.exp32)

gamma1 = con.exp32*(con.k2/(a1*a1))*factor

a2 = a1*(1-(con.exp13*gamma1)-(gamma1*gamma1)-(con.exp13481*gamma1*gamma1*gamma1))

gamma0 = con.exp32*(con.k2/(a2*a2))*factor

n01 = tle_ele[0]/(1+gamma0)

a01 = (con.ke/n01)**con.exp23

perigee = (a01*(1-e)-1)*con.ae
apogee = (a01*(1+e)-1)*con.ae
period = (2*pi*1440/1440)/n01
semimajoraxis = (con.mu/((n/60)**2))**(con.exp13)

#initialization for secular effects of atmospheric drag

if perigee < 220:
    #simplified eqn
    mode = 'simplified'

else:
    mode = 'normal'

if perigee < 156:
    s4 = perigee - 78

    if s4 < 20:
        s4 = 20

    q = ((120-s4)*(1/con.ae))**4
    s4 = s4/(con.ae+1)

else:
    s4 = con.ks
    q = con.q

tsi = 1/(a01-s4)
beta0 = (1-e2)**0.5
eta = a01*e*tsi

coeff1 = a01*(1+(1.5*eta*eta)+(4*eta*e)+(e*eta*eta*eta))
coeff2 = 1.5*con.k2*(tsi/(1-(eta*eta)))*(-0.5+(1.5*(angle*angle)))*(8+(24*eta*eta)+(3*eta*eta*eta*eta))

c2 = q*(tsi**4)*n01*((1-(eta*eta))**-3.5)*(coeff1+coeff2)

c1 = bstar*c2

coeff1 = 2*eta*(1+(e*eta))+(0.5*e)+(0.5*eta*eta*eta)
coeff2 =2*con.k2*tsi/(a01*(1-(eta*eta)))
coeff3 = 3*(1-(3*angle*angle))*(1+(1.5*eta*eta)-(2*e*eta)-(0.5*e*eta*eta*eta))
coeff4 = con.exp34*(1-(angle*angle))*((2*eta*eta)-(e*eta)-(e*eta*eta*eta))*(cos(2*omega))

c4 = 2*n01*q*(tsi**4)*a01*(beta0*beta0)*((1-(eta*eta))**-3.5)*(coeff1-coeff2*(coeff3+coeff4))

if mode =='simplified':
    c3 = 0
    c5 = 0
    deltaomega = 0
    d2 = 0
    d3 = 0
    d4 = 0
else:
    if e > con.ecc_all:
        c3 = q*(tsi**5)*con.a30*n01*con.ae*sin(i)/(con.k2*e)
    else:
        c3 = 0
    c5 = 2*q*(tsi**4)*a01*beta0*beta0*((1-(eta*eta))**-3.5)*(1+(con.exp114*eta*(eta+e))+(e*eta*eta*eta))
    deltaomega = bstar*c3*cos(omega)
    d2 = 4*a01*tsi*c1*c1
    d3 = con.exp43*a01*tsi*tsi*(17*a01+s4)*(c1**3)
    d4 = con.exp23*a01*a01*(tsi**3)*(221*a01+31*s4)*(c1**4)
    
#initialization for secular effects of earth zonal harmonics
coeff1 = 1.5*con.k2*(-1+(3*angle*angle))/(a01*a01*(beta0**3))
coeff2 = con.exp316*con.k2*con.k2*(13-78*(angle*angle)+137*(angle**4))/((a01**4)*(beta0**7))

mdot = (coeff1 + coeff2)*n01

coeff1 = -1.5*con.k2*(1-5*(angle*angle))/(a01*a01*(beta0**4))
coeff2 = con.exp316*con.k2*con.k2*(7-114*(angle*angle)+395*(angle**4))/((a01**4)*(beta0**8))
coeff3 = 1.25*con.k4*(3-36*(angle*angle)+49*(angle**4))/((a01**4)*(beta0**8))

omegadot = (coeff1 + coeff2 + coeff3)*n01

coeff1 = -3*con.k2*angle/(a01*a01*(beta0**4))
coeff2 = 1.5*con.k2*con.k2*(4*angle-19*(angle**3))/((a01**4)*(beta0**8))
coeff3 = 2.5*con.k4*angle*(3-7*(angle*angle))/((a01**4)*(beta0**8))

comega = (coeff1 + coeff2 + coeff3)*n01

if e > con.ecc_all:
    deltam = -1*con.exp23*q*bstar*(tsi**4)*con.ae/(e*eta)
else:
    deltam = 0

comegacoeff = -10.5*n01*con.k2*angle*c1/(a01*a01*beta0*beta0)
