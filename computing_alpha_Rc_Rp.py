import sympy as sp
import numpy as np
from scipy.optimize import fsolve

from scipy.integrate import dblquad

x, y, z, a, b, c = sp.symbols('x y z a b c')

sp.init_printing(use_unicode=True)


# For \alpha and R_c, need projection onto the yz-plane
z_fn = c*sp.sqrt(1 - (y/b)**2 )

# We'll need the first and second derivatives of z wrt y
dz_dy = sp.diff(z_fn, y)
d2z_dy2 = sp.diff(z_fn, y, 2)

# Let's first do \alplha. Need to evaluate z at the value of y
# for which z = c/2.

r_alpha = sp.Matrix([0, b*sp.sqrt(3)/2, c/2])

deriv_for_alpha = dz_dy.subs(y, r_alpha[1])


# This gives us an anlytical expression for z'(y). 
# Now let's determine an analytical expression for R_c.

# We first need to know where to evaluate 
# the first and second derivatives.

r_Rc = sp.Matrix([0, b*sp.sqrt(7)/4, 3*c/4])

# Recall that radius of curvature defined as (1 + (z')^2 )^(3/2) / z''

deriv_for_Rc = dz_dy.subs(y, r_Rc[1])
second_deriv_for_Rc = d2z_dy2.subs(y, r_Rc[1])

Rcurv = (1 + (deriv_for_Rc)**2 )**(3/2) / second_deriv_for_Rc


# Finally we will do the same procedure for R_p - the radius of
# curvature in the xy-plane. First we need y as a function of x

y_fn = b*sp.sqrt( sp.Rational(3/4) - (x/a)**2 )

# We will evaluate this radius of curvature at r_alpha

r_Rp = r_alpha

# Since y' is 0 at r_alpha, we will simply have R_p = y''.

d2y_dx2 = sp.diff(y_fn, x, 2)

Rp = d2y_dx2.subs(x, r_alpha[0])





