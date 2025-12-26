
# The first step is to choose the (physical) length scale. 
# This corresponds to the (physical) distance between lattice points. 
# This is chosen by the user. Here we will take a value motivated 
# by the paper by Yehuda and Bat-El

l_c = 7.5e-5

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

# Now, we have to determine u_{max} in lattice units. Since 
# we now have reference length and time scales, 
# we have 

U_c = l_c / t_c 

# Then, 

u_lu = u_max / U_c

