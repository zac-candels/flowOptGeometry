
import os, math, re, sys
import numpy as np
import struct
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("ignore")


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

# Slice in the z direction to plot
zslice=int(LZ/2)

for t in range(tstart + tinc, tend + 1, tinc):

    plt.close('all')

    file_name = datadir + "OrderParameter_t%li.mat"%t
    FileC = open(file_name, 'rb')
    dat = FileC.read()
    phi = np.ndarray((LX, LY, LZ), '=d', dat, 0, (8 * LY * LZ, 8 * LZ, 8))
    liquid = np.array(phi[:,:])
    # Set order parameter in the solid to 0.5 for visualisation
    # liquid[np.where(np.logical_or(solid == 1, solid == -1))[0], np.where(np.logical_or(solid == 1, solid == -1))[1], np.where(np.logical_or(solid == 1, solid == -1))[2]] = 0.
    # liquid[np.where(np.logical_or(solid == 3, solid == 2))[0], np.where(np.logical_or(solid == 3, solid == 2))[1], np.where(np.logical_or(solid == 3, solid == 2))[2]] = 0.
    FileC.close()
    phi = liquid[:,:,0]
    #phi = np.flipud(phi)
    
    
    full_interface_x = []
    full_interface_y = []
    phi_interface = []
    
    droplet_x = []
    droplet_y = []
    droplet_phi = []
        
    for i in range(len(phi[:,0])):
        for j in range(len(phi[0,:])):
            if phi[i,j] > 0.3 and phi[i,j] < 0.7:    
                #if j >= 25:
                full_interface_x.append(i)
                full_interface_y.append(j)    
                phi_interface.append(phi[i,j])
                  
    full_interface_x = np.asarray(full_interface_x)
    full_interface_y = np.asarray(full_interface_y)
    phi_interface = np.asarray(phi_interface)
    mid_pt_x = np.median(full_interface_x)
    
    # plt.figure()
    # plt.plot(droplet_x, droplet_y, 'o', color='r')
    # plt.xlim(0,LX)
    # plt.ylim(0, LY)
    # title_str = "t =" + str(t)
    # plt.title(title_str)
    
    mask_right = full_interface_x > mid_pt_x
    right_interface_x = full_interface_x[mask_right]
    right_interface_y = full_interface_y[mask_right]
    right_interface_phi = phi_interface[mask_right]
    
    mid_pt_y = np.median(full_interface_y)
    mask_lower = right_interface_y < mid_pt_y 
    right_interface_x = right_interface_x[mask_lower]
    right_interface_y = right_interface_y[mask_lower]
    right_interface_phi = right_interface_phi[mask_lower]
    
    CL_pos = 1e7
    CL_candidates = []
    for idx1 in range(len(right_interface_x)):
        right_neighbors = []
        left_neighbors = []
        
        x1 = right_interface_x[idx1]
        for idx2 in range(len(right_interface_x)):
            #print("x1 =", x1, ", x2 =", x2)
            x2 = right_interface_x[idx2]
            if x1 == x2: 
                continue
            if x2 > x1:
                right_neighbors.append(x2)
            elif x2 < x1:
                left_neighbors.append(x2)
        
        if left_neighbors == []:
            continue
        left_neighbors = np.sort(np.asarray(left_neighbors))
        right_neighbors = np.sort(np.asarray(right_neighbors))
        
        if x1 >= np.max(left_neighbors) + 2: 
            CL_candidates.append(x1)
            #CL_pos = x1
    #CL_pos = np.max(np.asarray(CL_candidates))
    
    print("CL_pos = ", CL_pos)
                
    
    
            
    
    plt.figure()
    plt.plot(full_interface_x, full_interface_y, 'o', color='r')
    plt.title('partial interface')
    plt.savefig("./interface.png")