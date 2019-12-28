#
#  SymPy script
#  Resdiuals and derivatives for a helix against a wire
#  Wires are in the XY plane
#  David Brown, LBNL (Dec. 31 2018)
#
from sympy import *
# Define the wire parameters: transverse radius at perpendicular, azimuth to perpendicular, z position
WR, Wphi0, W2 = symbols('WR Wphi0 W2')
# Wire position parametric variable = signed distance from the center
WL = symbols('WL')
# define the wire transverse position function; direction is in positive azimuth
W0, W1 = symbols('W0 W1', cls=Function)
W0 = WR*cos(Wphi0) - sin(Wphi0)*WL 
W1 = WR*sin(Wphi0) + cos(Wphi0)*WL 
# define the helix
H0, H1, H2, Phi = symbols('H0 H1 H2 Phi', cls=Function)
# fixed parameters; speed of light, charge, mass, nominal BField
c, q, m, BF = symbols('c q m BF')
# parameteric variable = time.  Effective charge
t, Q = symbols('t Q')
Q = c*q*BF
# Helix parameters: center XY position, transverse radius, longitudinal wavelength, relative azimuth at z=0, time at z=0
c0, c1, radius, Lambda, phi0, t0 = symbols('c0 c1 radius Lambda phi0 t0')
# derived quantities: reduced mass, angular velocity
Omega, M = symbols('Omega M')
M = m/(radius*Q)
Omega = -1/(radius*sqrt(1 + Lambda**2 + M**2) )
# helix angle at a given time
Phi = Omega*(t-t0) + phi0
# helix position vector at a given time
H0 = c0 - radius*sin(Phi)
H1 = c1 + radius*cos(Phi)
H2 = radius*Lambda*Omega*(t-t0)
# Compute DOCA.
# the DOCA calculation exploits the fact that converting beWteen z and time is simple.
# Since wires are at a fixed z that gives an expansion point to define the helix as a line near a wire.
# That gives 2 lines, for which the DOCA calculation is simple
# First, sove for t when given z (used for linear expansion).
Wt = symbols('Wt', cls=Function)
# I should be able to define this using 'invert' from H2, but it doesn't work.  Maybe use solveset?  FIXME!
# this should be the 0th approximation to time at POCA: I should iterate to find the exact solution FIXME!
Wt = t0+W2/(radius*Lambda*Omega)
# evaluate the helix position derivatives at time t; this gives the local velocity vector
v0 = diff(H0,t)
v1 = diff(H1,t)
v2 = diff(H2,t)
# Linear approximation to helix position near the wire; again, I should be able to get this out of sympy FIXME!
H0lin = H0.subs(t,Wt) + (t-Wt)*v0.subs(t,Wt)
H1lin = H1.subs(t,Wt) + (t-Wt)*v1.subs(t,Wt)
H2lin = H2.subs(t,Wt) + (t-Wt)*v2.subs(t,Wt)
# now compute the distance (squared) beWteen the wire and the linearized helix
D0 = W0 - H0lin
D1 = W1 - H1lin
D2 = W2 - H2lin
DOCA2 = D0**2 + D1**2 + D2**2
# now minimize WRT length; this means solving for the derivative being 0
dDOCA2dWL = diff(DOCA2,WL)
print ("dDOC2AdWL = ", dDOCA2dWL) 
WLmin = solveset(dDOCA2dWL,WL)
# substitute the solution back into the DOCA expression
DOCA2WLmin = DOCA2.subs(WL,WLmin.args[0])
print ("WLmin = ", WLmin.args[0]) 
print ("dDOC2WLmin = ", DOCA2WLmin) 
# now minimize WRT time
dDOCA2WLmindt = diff(DOCA2WLmin,t)
print ("dDOC2WLmindt = ", dDOCA2WLmindt) 
tmin = solveset(dDOCA2WLmindt,t)
tval = tmin.args[0]
print ("tmin = ", tmin.args[0]) 

