import sympy as sp

# setup flux and jacobian
heat_rat = sp.Symbol("gamma")
state = sp.symbols("q_0, q_1, q_2, q_3") # conserved variables
mass = state[2]
veloc_norm = state[0]/state[2]
veloc_tang = state[1]/state[2]
pres = (heat_rat - 1)*(state[3] - mass/2*(veloc_norm**2 + veloc_tang**2))
sound_speed = sp.sqrt(heat_rat*pres/mass)
flux = [
    mass*veloc_norm**2 + pres,
    mass*veloc_norm*veloc_tang,
    mass*veloc_norm,
    veloc_norm*(heat_rat*pres/(heat_rat - 1) + mass/2*(veloc_norm**2 + veloc_tang**2))
]
flux = [sp.simplify(f) for f in flux]
jacobian = sp.Matrix([[sp.simplify(flux[i_var].diff(state[j_var])) for j_var in range(4)] for i_var in range(4)])

# compute eigensystem
eigvals = list(jacobian.eigenvals())
print("eigenvalues:")
sp.pprint(eigvals)
print()
convective = (jacobian - eigvals[0]*sp.eye(4)).nullspace()
convective[0] += convective[1]*state[0]**2/(2*state[2]**2)
convective[0] = sp.simplify(convective[0]*state[2])
convective[1] *= state[1]
print("convective eigenvectors:")
sp.pprint(convective)
print()
print("accoustic eigenvectors:")
for i in [0, 1]:
    eigval = eigvals[i + 1]
    sign = 2*i - 1
    eigvec = (jacobian - eigval*sp.eye(4)).nullspace()[0]
    eigvec *= state[2]/eigvec[2]
    eigvec = sp.simplify(eigvec)
    sp.pprint(eigvec)
    diff = eigvec - state[2]*sp.Matrix([veloc_norm + sign*sound_speed, veloc_tang, 1, (state[3] + pres)/mass + sign*veloc_norm*sound_speed])
    diff = diff.subs(heat_rat, 1.4)
    vals = [4.6, 3.78, 1.02, 1.401e5]
    for i_var in range(4):
        diff = diff.subs(state[i_var], vals[i_var])
    print("error vs. expected expression: ", diff.evalf().norm())
