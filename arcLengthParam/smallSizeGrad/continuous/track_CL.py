
import os, math, re, sys
import numpy as np
import struct
from matplotlib import pyplot as plt
import warnings
import re
from pathlib import Path 
warnings.filterwarnings("ignore")


plt.close('all')
def coord_k(k, ly, lz):
    """From a k value, determines its xk, yk, and zk."""
    xk = math.floor(k/(ly*lz))
    yk = math.floor((k - xk*lz*ly)/lz)
    zk = k - xk*lz*ly - yk*lz
    return xk, yk, zk

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

data_path = Path(datadir)
pattern = re.compile(r"^OrderParameter_t(\d+)\.mat$")

tend = max(
    int(pattern.match(p.name).group(1))
    for p in Path(datadir).iterdir()
    if p.is_file() and pattern.match(p.name)
)

tinc = 2000

# Slice in the z direction to plot
zslice=int(LZ/2)

left_CL_init = 0
right_CL_init = 0
left_CL_disp_vec = []
right_CL_disp_vec = []
for t in range(tstart, tend, tinc):

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
    #phi = liquid[:,:,0]
    #phi = np.flipud(phi)
    
    rgbv = (liquid[:, :, zslice])
    phi = rgbv#[:-19, :]
    
    
    full_interface_x = []
    full_interface_y = []
    phi_interface = []
    
    droplet_x = []
    droplet_y = []
    droplet_phi = []
        
    for i in range(len(phi[:,0])):
        for j in range(len(phi[0,:])):
            if phi[i,j] > 0.4 and phi[i,j] < 0.6:
                if j > len(phi[0,:])/10 and j < 0.9*len(phi[0,:]):
                    #if j >= 25:
                    full_interface_x.append(i)
                    full_interface_y.append(j)    
                    phi_interface.append(phi[i,j])
                    
                  
    full_interface_x = np.asarray(full_interface_x)
    full_interface_y = np.asarray(full_interface_y)
    phi_interface = np.asarray(phi_interface)
    mid_pt_x = np.median(full_interface_x)
    
    left_CL = np.min(full_interface_x)
    right_CL = np.max(full_interface_x)
    
    if t == 0:
        left_CL_init = left_CL 
        right_CL_init = right_CL
    
    left_CL_disp = np.abs(left_CL - left_CL_init)
    right_CL_disp = np.abs(right_CL - right_CL_init)
    
    right_CL_disp_vec.append(right_CL_disp)
    left_CL_disp_vec.append(left_CL_disp)
    
    #print(f"Right CL displacement {right_CL_disp}")
    #print(f"Left CL displacement {left_CL_disp}\n\n")
    
            
    
    # plt.figure()
    # plt.plot(full_interface_x, full_interface_y, 'o', color='r')
    # title_str = f"interface, t={t}"
    # plt.title(title_str)
    # plt.savefig("./interface.png")
    # plt.show()
    
    a=1
    
right_CL_disp_vec = np.asarray(right_CL_disp_vec)
left_CL_disp_vec = np.asarray(left_CL_disp_vec)

t = range(tstart, tend, tinc)
plt.figure()
plt.plot(t, left_CL_disp_vec)
title_str = "Contact line spreading, continuous geometry"
plt.title(title_str)
plt.savefig("./spread.png")
