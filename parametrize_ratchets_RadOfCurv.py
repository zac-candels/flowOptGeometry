import sympy as sym
import numpy as np
from scipy.optimize import fsolve

x, y, z, a, b, c = sym.symbols('x y z a b c')

sym.init_printing(use_unicode=True)

y_fn = b*sym.sqrt(1 - (x/a)**2)

x_fn = a*sym.sqrt(1 - (z/c)**2)

dy = sym.diff(y_fn, x)

ddy = sym.diff(y_fn, x, 2)

ddx = sym.diff(x_fn, z, 2)

x_tilt = -sym.sqrt(3)*a/2

x_curv = -sym.sqrt(7)*a/4

slope_tilt_angle = dy.subs(x, x_tilt)

slope_R_curv = dy.subs(x, x_curv)

concavity_R_curv = ddy.subs(x, x_curv)

R2 = ddx.subs(z, 0)


R_curv = (1 + slope_R_curv**2)**(3/2) / concavity_R_curv


exp_tilt_angle = 60
exp_R_curv = 1
exp_R2 = 0.4

def func_yUp(x):
    angle = np.arctan(np.sqrt(3)*x[1]/x[0]) * 180 / np.pi - exp_tilt_angle
    R_curv = 27*x[0]**2*( 1 + (7/9)*x[1]**2/x[0]**2 )**(3/2) / (64*x[1]) - exp_R_curv
    R2 = x[2]**2 / x[0] - exp_R2
    return [angle, R_curv, R2]

[a, b, c] = fsolve(func_yUp, [1, 1, 1])


def func_zUp(x):
    angle = np.arctan(np.sqrt(3)*x[2]/x[1]) * 180 / np.pi - 60
    R_curv = 27*x[1]**2*( 1 + (7/9)*x[2]**2/x[1]**2 )**(3/2) / (64*x[2]) - 0.5
    R2 = x[0]**2 / x[1] - 0.4
    return [angle, R_curv, R2]

root2 = fsolve(func_zUp, [1, 1, 1])






