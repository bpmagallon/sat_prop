from sgp4_init import epoch, n, i, bstar, omega, m , node, e, a01, n01, mdot, nodedot, omegadot, deltaomega, deltam, nodecoeff, mode, c1, c3, c4, c5, d2, d3, d4
import constants as con
from math import sin, cos, fmod, pi, copysign, atan
import julian_date

#time since tle epoch
year = 2000+int(epoch[0:2])
day_sub = ((year - 1901)*1461)/4+365
day = epoch[3:]
jd1900 = julian_date.julian(1900, 1, 0, 0)
sgdp = jd1900 + day_sub + (day-1)
ep = (year - 1900)*1000+day

jd = 0.0 #date for propagation (year, month, day, time ut)
tse = (jd - sgdp)*con.mday


#propagation
maxi = 10
mold = m

#update for earth zonal gravity and partial atmospheric drag effects
mdf = m + tse*(n01+mdot)
wdf = omega + (omegadot*tse)
nodedf = node + (nodedot*tse)

node = nodedf - (nodecoeff*tse*tse)

if mode == 'simplified':
    a = a01*(1-(tse*tse))
    e = e - (bstar*tse*c4)
    xl = mdf + wdf + node + (n01*tse*tse*con.exp32)
else:
    m = mdf + deltaomega +deltam
    w = wdf = deltaomega - deltam
    
    coeff1 = 1-(c1*tse)-(d2*tse*tse)-(d3*tse*tse*tse)-(d4*tse*tse*tse*tse)
    a = a01*coeff1*coeff1

    coeff2 = bstar*(c4*tse+c5*(sin(m)-sin(mold)))
    e = e - coeff2

    coeff3 = (con.exp32*tse*tse)+((d2+(2*c1*c1))*tse*tse*tse)+(con.exp14*(3*d3+12*c1*d2+10*c1*c1*c1)*tse*tse*tse*tse)
    coeff4 = con.exp15*(3*d4+12*c1*d3+6*d2*d2+30*c1*c1*d2+15*c1*c1*c1*c1)*(tse**5)
    xl = m + w + node + (n01*coeff3*coeff4)
                        
beta2 = 1 - (e*e)

#update for long-period periodic effects of lunar and solar gravity
#skip, for deep-space only

#updage for long-period periodic effects of earth gravity

coeff5 = 1/(a*beta2)
a30vk2 = -1*con.j3/(con.k2*(con.ae**3))
xlcoeff = 0.125*a30vk2*sin(i)*((3)+(5*cos(i)))/(1+cos(i))
aycoeff = (0.25)*a30vk2*sin(i)

axn = e * cos(w)
ayn = e*sin(w)+aycoeff*coeff5
xlt = xl + axn*coeff5*xlcoeff


elsq =axn*axn + ayn*ayn

#newton-rapshson
epw = capu = fmod(xlt-node, 2*pi)
maxnr = elsq**0.5

for i in range(maxi):
    nr = 0.0
    f = 0.0
    df = 0.0

    ecosE = axn*cos(epw)+ayn*sin(epw)
    esinE = axn*sin(epw)-ayn*cos(epw)

    f = capu - epw + esinE
    if abs(f)<0.000000000001:
        break

    df = 1 - ecosE
    nr = f/df

    if (i==0) and (abs(nr)>(1.25*maxnr)):
        nr = copysign(maxnr, nr)

    else:
        nr = f/(df+0.5*esinE*nr)

    epw+=nr

#short period preliminary quantities
coeff = 1 - elsq
betal = coeff**0.5
pl = a*coeff
r = a*(1-ecosE)
invR = float(1)/r
coeff2 = a*invR
coeff3 = float(1)/(1+betal)
cosu = coeff2*(cos(epw)-axn+ayn*esinE*coeff3)
sinu = coeff2*(sin(epw)-ayn-axn*esinE*coeff3)
u = atan(sinu/cosu)
sin2u = 2*sinu*cosu
cos2u = 2*cosu*cosu-1

coeff = float(1)/pl
coeff2 = con.k2*coeff
coeff3 = coeff2*coeff1
coeff4 = 3*cos(i)*cos(i)-1
coeff5 = 1 -(cos(i)*cos(i))
coeff6 = 7*cos(i)*cos(i)-1

#update for short term periodics
rk = r*(1-1.5*coeff3*betal*coeff4)+0.5*coeff2*cos2u*coeff5
uk = u - 0.25*coeff3*sin2u*coeff6
xnodek = node + 1.5*coeff3*cos(i)*sin2u
xinck = i + (1.5*coeff3*cos(i)*sin(i)*cos2u)


radius = rk*con.xkmper/con.ae
theta = uk
eqinc = xinck
ascn = xnodek
argp = omega
smjaxs = a*con.xkmper/con.ae

#skip velocity computation
#transform kepler to xyz form

xmx = -sin(ascn)*cos(i)
xmy = cos(ascn)*cos(i)

ux = xmx * sin(theta) + cos(ascn)*cos(theta)
uy = xmy * sin(theta) + sin(ascn)*cos(theta)
uz = sin(ascn)*sin(theta)

#position
x = radius * ux
y = radius * uy
z = radius * uz

