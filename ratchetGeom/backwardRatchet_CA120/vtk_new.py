import numpy as np
import vtk
from vtk.util import numpy_support
import struct

# -------------------------------------------------------------------------
# Helper: write a numpy array as a legacy .vtk ImageData file
# -------------------------------------------------------------------------
def write_legacy_vtk(filename, array, array_name="scalars"):
    # Convert numpy array (Fortran order) to VTK array
    vtk_array = numpy_support.numpy_to_vtk(
        num_array=array.ravel(order="F"),
        deep=True,
        array_type=vtk.VTK_DOUBLE if array.dtype == np.float64 else vtk.VTK_INT,
    )

    # Create vtkImageData and assign data
    img = vtk.vtkImageData()
    img.SetDimensions(array.shape)
    img.SetOrigin(0.0, 0.0, 0.0)
    img.GetPointData().SetScalars(vtk_array)
    vtk_array.SetName(array_name)

    # Write legacy .vtk format
    writer = vtk.vtkDataSetWriter()
    writer.SetFileName(filename)
    writer.SetInputData(img)
    writer.SetFileTypeToBinary()
    writer.Write()


# -------------------------------------------------------------------------
# Read geometry dimensions
# -------------------------------------------------------------------------
with open("./data/Header.mat", "rb") as f:
    LX = struct.unpack("=i", f.read(4))[0]
    LY = struct.unpack("=i", f.read(4))[0]
    LZ = struct.unpack("=i", f.read(4))[0]

# -------------------------------------------------------------------------
# BoundaryLabels_t0.mat → solid.vtk
# -------------------------------------------------------------------------
with open("./data/BoundaryLabels_t0.mat", "rb") as f:
    dat = f.read()
solid = np.ndarray((LX, LY, LZ), dtype="=i", buffer=dat, offset=0,
                   strides=(4*LY*LZ, 4*LZ, 4)).copy()
solid[solid == -1] = 1

write_legacy_vtk("solid.vtk", solid, array_name="BoundaryLabels_t0")

# -------------------------------------------------------------------------
# OrderParameter_t0.mat → phi.vtk
# -------------------------------------------------------------------------
with open("./data/OrderParameter_t0.mat", "rb") as f:
    dat = f.read()
phi = np.ndarray((LX, LY, LZ), dtype="=d", buffer=dat, offset=0,
                 strides=(8 * LY * LZ, 8 * LZ, 8))

write_legacy_vtk("phi.vtk", phi, array_name="OrderParameter_t0")
