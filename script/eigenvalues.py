## \cond

import numpy as np

def flux(state):
    pres = .4*(state[2] - .5*state[0]**2/state[1])
    f = np.zeros(3)
    f[0] = state[0]**2/state[1] + pres
    f[1] = state[0]
    f[2] = state[0]/state[1]*(state[2] + pres)
    return f

mass = 1.225
pres = 101325
sound_speed = (1.4*pres/mass)**.5
#veloc = sound_speed
veloc = 12.
state = np.array([veloc*mass, mass, pres/.4 + .5*mass*veloc**2])
d_mass = 1
d_veloc = sound_speed/mass*d_mass
d_pres = sound_speed**2*d_mass
d_state = np.array([
    d_mass*veloc + mass*d_veloc,
    d_mass,
    d_pres/.4 + .5*d_mass*veloc**2 + mass*veloc*d_veloc,
])
d = 1e-6
d_flux = (flux(state + d*d_state) - flux(state))/d
print(d_state)
print(d_flux)
print(d_flux/d_state)
print(veloc + sound_speed)

## \endcond
