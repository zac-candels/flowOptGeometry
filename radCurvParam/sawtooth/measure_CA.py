import os, math, re, sys
import numpy as np
import struct
from matplotlib import pyplot as plt

import warnings
warnings.filterwarnings("ignore")

def read_params(path):
    params = {}
    # matches:   key   =   integer   (ignoring anything after a '#')
    pat = re.compile(r'^\s*([A-Za-z_]\w*)\s*=\s*([0-9]+)')
    with open(path) as f:
        for line in f:
            m = pat.match(line)
            if m:
                key, val = m.group(1), int(m.group(2))
                params[key] = val
    return params

# Usage
path = "."
path_input = path + "/input.txt"
params = read_params(path_input)
saveInterval = params["saveInterval"]
Num_steps   = params["timesteps"]
postHeight = params["postheight"]


plt.close('all')
def coord_k(k, ly, lz):
    """From a k value, determines its xk, yk, and zk."""
    xk = math.floor(k/(ly*lz))
    yk = math.floor((k - xk*lz*ly)/lz)
    zk = k - xk*lz*ly - yk*lz
    return xk, yk, zk

def interpolate_x_direction(pts_x, pts_phi):
    pts_x = np.asarray(pts_x)
    pts_phi = np.asarray(pts_phi)
    gas_phase_x, gas_phase_phi = [], []
    liquid_phase_x, liquid_phase_phi = [], []
    if len(pts_x) > 2:
        #sorted_indices = np.argsort(pts_phi)
        #pts_phi = pts_phi[sorted_indices]
        #pts_x = pts_x[sorted_indices]
        mask_gas = pts_phi < 0.5
        mask_liquid = pts_phi >= 0.5
        
        gas_phase_phi = pts_phi[mask_gas]
        gas_phase_x = pts_x[mask_gas]
        
        liquid_phase_phi = pts_phi[mask_liquid]
        liquid_phase_x = pts_x[mask_liquid]
        
        if len(liquid_phase_phi) == 0:
            closest_gas_index = np.argmax(gas_phase_phi)
            closest_gas_pt = gas_phase_x[closest_gas_index]
            # Nothing to interpolate
            return closest_gas_pt
        elif len(gas_phase_phi) == 0:
            closest_liquid_index = np.argmin(liquid_phase_phi)
            closest_liquid_pt = liquid_phase_x[closest_liquid_index]
            return closest_liquid_pt
        else:
            closest_gas_index = np.argmax(gas_phase_phi)
            closest_gas_pt = gas_phase_x[closest_gas_index]
            
            closest_liquid_index = np.argmin(liquid_phase_phi)
            closest_liquid_pt = liquid_phase_x[closest_liquid_index]
        
        x0, x1 = [closest_liquid_pt, closest_gas_pt]
        phi0, phi1\
            = [gas_phase_phi[closest_gas_index],liquid_phase_phi[closest_liquid_index]]
            
    if len(pts_x) == 2:
        x0, x1 = pts_x[0], pts_x[1]
        phi0, phi1 = pts_phi[0], pts_phi[1]
    
    m = ( phi1 - phi0 )\
        /( x1 - x0 )
    x_interpolated = x0 + (1/m)*(0.5 - phi0)
    
    return x_interpolated


def interpolate_y_direction(pts_y, pts_phi):
    pts_y = np.asarray(pts_y)
    pts_phi = np.asarray(pts_phi)
    gas_phase_y, gas_phase_phi = [], []
    liquid_phase_y, liquid_phase_phi = [], []
    if len(pts_y) > 2:
        #sorted_indices = np.argsort(pts_phi)
        #pts_phi = pts_phi[sorted_indices]
        #pts_x = pts_x[sorted_indices]
        mask_gas = pts_phi < 0.5
        mask_liquid = pts_phi >= 0.5
        
        gas_phase_phi = pts_phi[mask_gas]
        gas_phase_y = pts_y[mask_gas]
        
        liquid_phase_phi = pts_phi[mask_liquid]
        liquid_phase_y = pts_y[mask_liquid]
        
        if len(liquid_phase_phi) == 0:
            closest_gas_index = np.argmax(gas_phase_phi)
            closest_gas_pt = gas_phase_y[closest_gas_index]
            # Nothing to interpolate
            return closest_gas_pt
        elif len(gas_phase_phi) == 0:
            closest_liquid_index = np.argmin(liquid_phase_phi)
            closest_liquid_pt = liquid_phase_y[closest_liquid_index]
            return closest_liquid_pt
        else:
            closest_gas_index = np.argmax(gas_phase_phi)
            closest_gas_pt = gas_phase_y[closest_gas_index]
            
            closest_liquid_index = np.argmin(liquid_phase_phi)
            closest_liquid_pt = liquid_phase_y[closest_liquid_index]
        
        y0, y1 = [closest_liquid_pt, closest_gas_pt]
        phi0, phi1\
            = [gas_phase_phi[closest_gas_index],liquid_phase_phi[closest_liquid_index]]
            
    if len(pts_y) == 2:
        y0, y1 = pts_y[0], pts_y[1]
        phi0, phi1 = pts_phi[0], pts_phi[1]
    
    m = ( phi1 - phi0 )\
        /( y1 - y0 )
    y_interpolated = y0 + (1/m)*(0.5 - phi0)
    
    return y_interpolated



plt.close('all')

# Directory where data is stored
datadir = "./data/"

# Open the Header file
HeaderFile = open(datadir+"Header.mat", 'rb')

# Simulation size, struct.unpack reads 4 bytes from Header.mat and interprets them as an integer
LX = struct.unpack('=i', HeaderFile.read(4))[0]
# Read the next 4 bytes...
LY = struct.unpack('=i', HeaderFile.read(4))[0]
LZ = struct.unpack('=i', HeaderFile.read(4))[0]

# 2D or 3D
ndim = struct.unpack('=i', HeaderFile.read(4))[0]

# What time to start at, what time to end at, and the time increment
tstart = 0
tend = struct.unpack('=i', HeaderFile.read(4))[0]
tinc = struct.unpack('=i', HeaderFile.read(4))[0]
tinc = 10000

def coord_k(k, ly, lz):
    """From a k value, determines its xk, yk, and zk."""
    xk = math.floor(k/(ly*lz))
    yk = math.floor((k - xk*lz*ly)/lz)
    zk = k - xk*lz*ly - yk*lz
    return xk, yk, zk
#pressure = np.zeros((tend, LX, LY, LZ))

# Slice in the z direction to plot
zslice=0

# Where to save the plots
outDirName = "figures"
os.system("mkdir -p %s"%outDirName)

print("tend=",tend)

t = 400000
# Boundary file
file_name = datadir + f"BoundaryLabels_t{t}.mat"
with open(file_name, 'rb') as f:
    dat = f.read()
solid = np.ndarray((LX, LY, LZ), '=i', dat, 0, (4 * LY * LZ, 4 * LZ, 4))

# Order parameter
file_name = datadir + f"OrderParameter_t{t}.mat"
with open(file_name, 'rb') as f:
    dat = f.read()
C = np.ndarray((LX, LY, LZ), '=d', dat, 0, (8 * LY * LZ, 8 * LZ, 8))
liquid = np.array(C[:,:])
# liquid[np.logical_or(solid == 1, solid == -1)] = 0.5
# liquid[np.logical_or(solid == 3, solid == 2)] = 0.5

phi = liquid[:,:,0]


full_interface_x = []
full_interface_y = []
phi_interface = []
    
for i in range(len(phi[:,0])):
    for j in range(len(phi[0,:])):
        if phi[i,j] > 0.3 and phi[i,j] < 0.7:    
            if j <= 12 and j > 2:
                full_interface_x.append(i)
                full_interface_y.append(j)    
                phi_interface.append(phi[i,j])
            
full_interface_x = np.asarray(full_interface_x)
full_interface_y = np.asarray(full_interface_y)
phi_interface = np.asarray(phi_interface)
mid_pt = np.median(full_interface_x)


plt.figure()
plt.plot(full_interface_x, full_interface_y, 'o', color='r')
plt.title('Full interface')
plt.savefig("./interface.png")


# plt.figure()
# plt.plot(left_interface_x, left_interface_y, 'o')
# plt.title('left interface')


tol = 1e-5

interpolated_interface_x = []
interpolated_interface_y = []

# Go through each point (x,y,phi) in the interface. If for any value of y
# there are multiple x values, interpolate the value of phi between them
# to get a single value. This new value of x will replace the two previous ones.
for i in range(len(full_interface_x)):
    x_i, y_i, phi_i = full_interface_x[i], full_interface_y[i], phi_interface[i]
    pts_x, pts_y, pts_phi = [x_i], y_i, [phi_i]
    for j in range(len(full_interface_x)):
        # for each fixed point P_i, loop through every other 
        # point P_j to see if there are other points that have the 
        # same y-value but different x-values. That is to say,
        # there may be multiple points at the same y-value. We collect all of these into an array,
        # select the two closest to 0.5 - but on opposite sides of 0.5 - and do the same routine as
        # mentioned above - that is to say, interpolate phi between them and return a single point.
        if j == i:
            continue
        x_j, y_j, phi_j\
            = full_interface_x[j], full_interface_y[j], phi_interface[j]
        if abs(y_j - y_i) < tol: # if point P_j is at the same y-level as point P_i
            if abs(x_i - x_j) > 2*tol: # if the points are at different x-levels, create arrays pts_x and pts_phi
                pts_x.append(x_j), pts_phi.append(phi_j)

    if len(pts_x) == 1: # If there is only one point at the given y-level, there's nothing further to do.
                        # Include this point in the new interface.
        interpolated_interface_x.append(x_i)
        interpolated_interface_y.append(y_i)
    else: # If there are multiple points at the same y-level, interpolate the values of phi between them and return a single point, x_new.
        # add x_new to the new interface.
        x_new = interpolate_x_direction(pts_x, pts_phi)
        interpolated_interface_x.append( x_new )
        interpolated_interface_y.append( y_i )
        
        
interpolated_interface_x = np.asarray(interpolated_interface_x)
interpolated_interface_y = np.asarray(interpolated_interface_y)
        


plt.figure()
plt.plot(interpolated_interface_x, interpolated_interface_y, 'o')
plt.savefig("./figures/interface.png", dpi=400, format='png')





