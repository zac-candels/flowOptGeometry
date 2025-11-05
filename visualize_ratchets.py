import numpy as np

# Try both skimage API names for compatibility
try:
    from skimage.measure import marching_cubes
except Exception:
    from skimage.measure import marching_cubes_lewiner as marching_cubes

# -----------------------------
# Parameters
# -----------------------------
a1, b1, c1 = 1.44, 4.23, 2.58
w = 0.02
a2, b2, c2 = a1 + w, b1 + w, c1 + w       # a2 = 6 (z half-length), c2 = 7 (xy radius)
angular_width_deg = 10
theta_half = np.deg2rad(angular_width_deg) / 2.0

x0, y0, z0 = 50.0, 50.0, 0.0  # center

# -----------------------------
# Grid bounds
# -----------------------------
margin = 1.5
x_min, x_max = x0 - (c2 + margin), x0 + (c2 + margin)
y_min, y_max = y0 - (c2 + margin), y0 + (c2 + margin)
z_min, z_max = z0 - (a2 + margin), z0 + (a2 + margin)

nx = ny = 400
nz = 400

x_vals = np.linspace(x_min, x_max, nx)
y_vals = np.linspace(y_min, y_max, ny)
z_vals = np.linspace(z_min, z_max, nz)

dx = x_vals[1] - x_vals[0]
dy = y_vals[1] - y_vals[0]
dz = z_vals[1] - z_vals[0]

X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing="ij")

# -----------------------------
# Implicit fields
# -----------------------------
F_outer = (b2 * c2)**2 * (X - x0)**2 + (a2 * c2)**2 * (Y - y0)**2\
    + (a2 * b2)**2 * (Z - z0)**2 - (a2 * b2 * c2)**2
F_inner = (b1 * c1)**2 * (X - x0)**2 + (a1 * c1)**2 * (Y - y0)**2\
    + (a1 * b1)**2 * (Z - z0)**2 - (a1 * b1 * c1)**2

H_half  = (y0 - Y)                  
H_ang   = np.abs(np.arctan2((X-x0), (Y-y0))) - theta_half  
H_z     = (z0 - Z)                  

# Combined field
G = np.maximum.reduce([F_outer, -F_inner, H_half, H_ang, H_z])

# -----------------------------
# Marching cubes on G = 0
# -----------------------------
verts, faces, normals, values = marching_cubes(G, level=0.0, spacing=(dx, dy, dz))

# Shift vertices into world coordinates
verts[:, 0] += x_min
verts[:, 1] += y_min
verts[:, 2] += z_min

# -----------------------------
# Write VTK (legacy PolyData)
# -----------------------------
vtk_file = "wedge.vtk"
with open(vtk_file, "w") as f:
    f.write("# vtk DataFile Version 3.0\n")
    f.write("Wedge region between ellipsoids\n")
    f.write("ASCII\n")
    f.write("DATASET POLYDATA\n")
    f.write(f"POINTS {len(verts)} float\n")
    for v in verts:
        f.write(f"{v[0]} {v[1]} {v[2]}\n")

    f.write(f"POLYGONS {len(faces)} {len(faces) * 4}\n")
    for face in faces:
        f.write(f"3 {face[0]} {face[1]} {face[2]}\n")

print(f"✅ Wrote VTK file: {vtk_file}")

