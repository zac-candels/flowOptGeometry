import sympy as sp
import numpy as np
from scipy.optimize import fsolve

theta = sp.symbols('theta ')
a, b, c, phi = sp.symbols('a b c phi', positive=True)

range_cond = sp.And(phi >= 0, phi <= sp.pi)


x = a*sp.cos(theta)*sp.sin(phi)
y = b*sp.sin(theta)*sp.sin(phi)
z = c*sp.cos(phi)

dx_dphi = sp.diff(x, phi)
dy_dphi = sp.diff(y, phi)
dz_dphi = sp.diff(z, phi)

dx_dtheta = sp.diff(x,theta)
dy_dtheta = sp.diff(y,theta)
dz_dtheta = sp.diff(z,theta)


T_theta = sp.Matrix([ dx_dtheta, dy_dtheta, dz_dtheta ])

T_phi = sp.Matrix([ dx_dphi, dy_dphi, dz_dphi ])

cross_prod = T_phi.cross(T_theta)

cross_prod_mag = cross_prod.norm()


chatgpt_expr = a*sp.sin(phi)*sp.sqrt(a**2 * sp.cos(phi)*sp.cos(phi) + c**2 * sp.sin(phi)*sp.sin(phi)   )

equal = sp.simplify( chatgpt_expr - cross_prod_mag )