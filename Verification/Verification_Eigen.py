# Citation 550 - Linear simulation

import scipy as np
from math import pi,sin,cos
from math import radians as rad
# xcg = 0.25 * c

# Stationary flight condition

hp0    =  1530.1    	      # pressure altitude in the stationary flight condition [m]
V0     =  134.56        # true airspeed in the stationary flight condition [m/sec]
alpha0 =  rad(1.8)          # angle of attack in the stationary flight condition [rad]
th0    =  rad(2.2)       # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      =  6612           # mass [kg]

# aerodynamic properties
e      =   0.838         # Oswald factor [ ]
CD0    =   0.023          # Zero lift drag coefficient [ ]
CLa    =   4.910          # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    =   0.07620         # longitudinal stabilty [ ]
Cmde   =   0.02540          # elevator effectiveness [ ]

# Aircraft geometry

S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3] 
lamda = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# air density [kg/m^3]  
rho    = rho0 * pow( ((1+(lamda * hp0 / Temp0))), (-((g / (lamda*R)) + 1)))   

W      = m * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.02792
CXa    = +0.47966		# Positive! (has been erroneously negative since 1993) 
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu    = -0.37616
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415

CYb    = -0.7500
CYbdot =  0     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0     
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939


symcor = V0/c               #assymetric correction factor
asymcor= V0/b               #assymetric correction factor


#ABC formula with imaginary part
def abcimag(A,B,C,cor):
    return (-B*cor+cor*np.absolute(np.sqrt(4*A*C-B**2))*1j)/(2*A),(-B*cor-cor*np.absolute(np.sqrt(4*A*C-B**2))*1j)/(2*A)


#SHORT PERIOD
ac1    = 4*muc**2*KY2
bc1    = -2*muc*(KY2*CZa+Cmadot+Cmq)
cc1    = CZa*Cmq -2*muc*Cma
lamdac1= abcimag(ac1,bc1,cc1,symcor)

#PHUGOID
ac2    = 2*muc*(CZa*Cmq - 2*muc*Cma)
bc2    = 2*muc*(CXu*Cma-Cmu*CXa) + Cmq*(CZu*CXa-CXu*CZa)
cc2    = CZ0 * (Cmu*CZa - CZu*Cma)
lamdac2= abcimag(ac2,bc2,cc2,symcor)

#APERIODIC ROLL
lamdab1 = (Clp/(4*mub*KX2))*asymcor

#DUTCH ROLL
ab2    = 8 * mub**2*KZ2
bb2    = -2*mub*(Cnr+2*KZ2*CYb)
cb2    = 4*mub*Cnb + CYb*Cnr
lamdab2= abcimag(ab2,bb2,cb2,asymcor)

#SPIRAL
lamdab3 = (2*CL*(Clb*Cnr - Cnb*Clr)/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr + 4* mub*Clb)))*asymcor

