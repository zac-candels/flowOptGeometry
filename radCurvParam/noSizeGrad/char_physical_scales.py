
# The first step is to choose the (physical) length scale. 
# This corresponds to the (physical) distance between lattice points. 
# This is chosen by the user. Here we will take a value motivated 
# by the paper by Yehuda and Bat-El

l_c = 2e-5

rho_c = 1000

# Now we must choose the non-dimensional relaxation time, tau_star. 
# A simple procedure is to first try some value \tau_star > 0.5 such that 
# \tau_star = O(1). We will first try 0.6.

tau = 0.6

# We must know the relevant length, velocity, and time scales. 
# Here (based on experimental data) we have u_{max} = 0.5mm/s. 

u_max = 0.0005

# Further, since nu = 10e-6 m^2/s, we have 
nu = 10e-6

t_c = (1./3.)*(tau - 0.5)*l_c**2 / nu 

print("t_c = ", t_c)

# Now, we have to determine u_{max} in lattice units. Since 
# we now have reference length and time scales, 
# we have 

U_c = l_c / t_c 

# Then, 

u_lu = u_max / U_c

print("u_max in lattice units is", u_lu)

# If u_lu < 1/3, we're ok. Now, we also need
# the conversion factor for the surface tension.
# The conversion factor is \sigma_c = \rho_c * l_c^3 * t_c^{-2}. 

sigma_c = rho_c * l_c**3 * t_c**(-2)

sigma_bar = 0.0728 / sigma_c

print("sigma_bar = ", sigma_bar)
