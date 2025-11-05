import numpy as np
from scipy.integrate import quad
from scipy.optimize import root

# In this program, S1 and S2 correspond to the arclengths of the ratchet 
# in the xz and xy planes. 

#%% Compute a, b, c, from alpha, S1, S2, theta
def F(vars, alpha, S1, S2):
    a, c, theta = vars

    # Eqn 1
    f1 = np.arctan(c/a) - alpha

    # Eqn 2
    integrand1 = lambda z: np.sqrt(1 + (a**2/c**2 - 1) * np.sin(z)**2)
    I1, _ = quad(integrand1, np.pi/2, 0)
    f2 = -c * I1 - S1

    # Eqn 3
    f3 = a*theta - S2

    return np.array([f1, f2, f3])

# Example parameters
alpha = 48.37*np.pi/180
S1 = 66.82
s2 = 27.92

# Initial guess

x0 = [1, 1, 1]

sol = root(F, x0, args=(alpha, R1, R2))

print("Success or Failure:", sol.success)
print("Solution (a, c, theta):", sol.x[0], sol.x[1], sol.x[2]*180/np.pi)

a = sol.x[0]
c = sol.x[1]
theta = sol.x[2]*180/np.pi


#%% Compute alpha, R1, R2 for given values of a, c, theta

a = 40
c = 45
theta = 40*np.pi/180

alpha = np.arctan(c/a)*180/np.pi

integrand1 = lambda z: np.sqrt(1 + (a**2/c**2 - 1) * np.sin(z)**2)
I1, _ = quad(integrand1, np.pi/2, 0)
S1 = -c * I1

S2 = a*theta

