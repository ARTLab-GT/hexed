import sympy as sp

heat_rat = sp.Symbol("gamma")
state = sp.symbols("q_0, q_1, q_2, q_3")
mass = state[2]
veloc_norm = state[0]/state[2]
veloc_tang = state[1]/state[2]
pres = (heat_rat - 1)*(state[3] - mass/2*(veloc_norm**2 + veloc_tang**2))
flux = [
    mass*veloc_norm**2 + pres,
    mass*veloc_norm*veloc_tang,
    mass*veloc_norm,
    veloc_norm*(heat_rat*pres/(heat_rat - 1) + mass/2*(veloc_norm**2 + veloc_tang**2))
]
flux = [sp.simplify(f) for f in flux]
jacobian = sp.Matrix([[sp.simplify(flux[i_var].diff(state[j_var])) for j_var in range(4)] for i_var in range(4)])
eigvals = list(jacobian.eigenvals())
convective = (jacobian - eigvals[0]*sp.eye(4)).nullspace()
convective[0] += convective[1]*state[0]**2/(2*state[2]**2)
convective[0] = sp.simplify(convective[0])
convective[1] *= state[1]
print("convective eigenvectors:")
sp.pprint(convective)
"""
accoustic0 = (jacobian - eigvals[1]*sp.eye(4)).nullspace()
accoustic1 = (jacobian - eigvals[2]*sp.eye(4)).nullspace()
sp.pprint(sp.simplify(accoustic0[0]))
sp.pprint(sp.simplify(accoustic1[0]))
"""
