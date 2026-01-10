import sympy as sym
import numpy as np
from scipy.optimize import fsolve

from scipy.integrate import dblquad


exp_tilt_angle = 60
exp_R_curv = 0.75
exp_R2 = 0.4
exp_Sa = 0.2


def integrand(theta, phi, a, b, c):
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
    term1 = (a*c)**2 * (sin_phi**4) * np.abs(sin_theta)**2
    
    # term2 = b^2 * sin^4(phi) * |c * cos(theta)|^2
    term2 = (b*c)**2 * (sin_phi**4) * np.abs( cos_theta)**2
    
    # term3 = |a^2 * sin(phi) * sin^2(theta) * cos(phi) + a^2 * sin(phi) * cos(phi) * cos^2(theta)|^2
    inner_term3 = (a*b * sin_phi * (sin_theta**2) * cos_phi) + \
                  (a*b * sin_phi * cos_phi * (cos_theta**2))
    term3 = np.abs(inner_term3)**2
    
    return np.sqrt(term1 + term2 + term3)

def surfaceArea(a, b, c, theta_L):
    
    phi_lims = [0, np.pi/3]
    
    theta_lims = [np.pi/2-theta_L/2, np.pi/2+theta_L/2]
    
    
    area = dblquad(integrand, phi_lims[0], phi_lims[1],
                   lambda x: theta_lims[0], lambda x: theta_lims[1], args=(a, b, c) )[0]
    
    return area




def func_zUp(x):
    angle = np.arctan(x[2]*np.sqrt(3)/x[1])*180/np.pi - exp_tilt_angle
    
    R_curv = 27*x[1]**2*( 1 + (7/9)*x[2]**2/x[1]**2 )**(3/2) / (64*x[2]) - exp_R_curv
    
    R2 = np.sqrt(3)*x[0]**2 / (2*x[1]) - exp_R2
    
    Sa = surfaceArea(x[0], x[1], x[2], x[3]) - exp_Sa
    return [angle, R_curv, R2, Sa]

[a, b, c, theta_L] = fsolve(func_zUp, [1, 1, 1, 1])


print("theta_L = ", theta_L*180/np.pi)



