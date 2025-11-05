import numpy as np
from tvtk.api import tvtk, write_data
import struct

file_name = "./data/BoundaryLabels_t0.mat"

# Open the Header file
HeaderFile = open("./data/Header.mat", 'rb')

# Simulation size, struct.unpack reads 4 bytes from Header.mat and interprets them as an integer
LX = struct.unpack('=i', HeaderFile.read(4))[0]
# Read the next 4 bytes...
LY = struct.unpack('=i', HeaderFile.read(4))[0]
LZ = struct.unpack('=i', HeaderFile.read(4))[0]

File = open(file_name, 'rb')
dat=File.read()
solid = np.array(np.ndarray((LX, LY, LZ), '=i', dat, 0, (4* LY * LZ, 4 * LZ, 4)))
solid[solid==-1]=1


#rho=np.concatenate((rho,rho),axis=0)
#rho=np.concatenate((rho,rho),axis=1)
#rho=np.concatenate((rho,rho),axis=2)

grid = tvtk.ImageData(origin=(0,0,0), #spacing=(10, 5, -10)
                      dimensions=solid.shape)
grid.point_data.scalars = solid.ravel(order='F')
grid.point_data.scalars.name = file_name

# Writes legacy ".vtk" format if filename ends with "vtk", otherwise
# this will write data using the newer xml-based format.
write_data(grid, 'solid.vtk')

#######################################################################


file_name = "./data/OrderParameter_t0.mat"


File = open(file_name, 'rb')
dat=File.read()
phi = np.ndarray((LX, LY, LZ), '=d', dat, 0, (8* LY * LZ, 8 * LZ, 8))

#rho=np.concatenate((rho,rho),axis=0)
#rho=np.concatenate((rho,rho),axis=1)
#rho=np.concatenate((rho,rho),axis=2)

grid = tvtk.ImageData(origin=(0,0,0), #spacing=(10, 5, -10)
                      dimensions=phi.shape)
grid.point_data.scalars = phi.ravel(order='F')
grid.point_data.scalars.name = file_name

# Writes legacy ".vtk" format if filename ends with "vtk", otherwise
# this will write data using the newer xml-based format.
write_data(grid, 'phi.vtk')