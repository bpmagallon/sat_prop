#constants

mu = 398600 #gravitational parameter
ae = 6378.135 #earth radius
ke = 0.0743669161 # sqrt of GM
q0 = 120
q = ((q0-78)*1/ae)**4
ks = 1*(1+(78/ae))
j2 = 0.001082616
j3 = -0.00000253881
j4 = -0.00000165597
k2 = 0.5*(ae*ae)*j2
k4 = -(float(3)/float(8))*j4*(ae**4)
a30 = -j3*(ae**3)
ecc_all = 0.0001
exp32 = float(3)/float(2)
exp316 = float(3)/float(16)
exp23 = float(2)/float(3)
exp43 = float(4)/float(3)
exp13 = float(1)/float(3)
exp34 = float(3)/float(4)
exp114 = float(11)/float(4)
exp13481 = float(134)/float(81)
