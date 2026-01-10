import sympy as sym
import numpy as np
from scipy.optimize import fsolve

from scipy.integrate import dblquad


exp_tilt_angle = 60
exp_R_curv = 1
exp_R2 = 0.4
exp_Sa = 0.674


def integrand(theta, phi, a=1.0, b=1, c=1.0):
    """
    The function from the image. 
    Note: scipy.integrate.dblquad expects the first argument to be the inner 
    variable (theta) and the second to be the outer variable (phi).
    """
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    
    # term1 = a^2 * sin^4(phi) * |c * sin(theta)|^2
    term1 = (a**2) * (sin_phi**4) * np.abs(c * sin_theta)**2
    
    # term2 = a^2 * sin^4(phi) * |c * cos(theta)|^2
    term2 = (a**2) * (sin_phi**4) * np.abs(c * cos_theta)**2
    
    # term3 = |a^2 * sin(phi) * sin^2(theta) * cos(phi) + a^2 * sin(phi) * cos(phi) * cos^2(theta)|^2
    inner_term3 = (a**2 * sin_phi * (sin_theta**2) * cos_phi) + \
                  (a**2 * sin_phi * cos_phi * (cos_theta**2))
    term3 = np.abs(inner_term3)**2
    
    return np.sqrt(term1 + term2 + term3)

def surfaceArea(a, b, c, theta_L):
    
    area = dblquad(integrand, 0, np.pi/3, 0, theta_L, args=(a, b, c) )[0]
    
    return area
    

# NEED TO FIX INTEGRAND. CHANGE THE AREA DIFFERENTIAL 
# SO THAT IT INCLUDES B AS WELL.



def func_zUp(x):
    angle = x[2]*np.sqrt(3) - exp_tilt_angle
    R_curv = 27*x[0]**2*( 1 + (7/9)*x[2]**2/x[0]**2 )**(3/2) / (64*x[2]) - exp_R_curv
    R2 = 64 / (27*x[0]) - exp_R2
    Sa = surfaceArea(x[0], x[1], x[2], x[3]) - exp_Sa
    return [angle, R_curv, R2, Sa]

[a, b, c, theta_L] = fsolve(func_zUp, [1, 1, 1, 1])




