#constants
mday = 1440
mu = 398600 #gravitational parameter
xkmper = 6378.135 #earth radius
ae = 1 #earth unit
ke = 0.0743669161331734132 # sqrt of GM
q0 = 120
q = ((q0-78)*ae/xkmper)**4
ks = ae*(1+78/xkmper)
j2 = 0.001082616
j3 = -0.00000253881
j4 = -0.00000165597
k2 = 0.5*(ae*ae)*j2
k4 = -(float(3)/float(8))*j4*(ae**4)
a30 = -j3*(ae**3)
ecc_all = 0.0001

exp15 = float(1)/float(5)
exp32 = float(3)/float(2)
exp14 = float(1)/float(4)
exp316 = float(3)/float(16)
exp23 = float(2)/float(3)
exp43 = float(4)/float(3)
exp13 = float(1)/float(3)
exp34 = float(3)/float(4)
exp114 = float(11)/float(4)
exp13481 = float(134)/float(81)
