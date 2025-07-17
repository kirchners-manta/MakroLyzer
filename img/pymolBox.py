from pymol import cmd

cmd.load("~/phd/plastic/PolyLyzer/polymer_analysis/example_polymers/PE/07.xyz", "molecule")

# Calculate the middle of geometry of the molecule
model = cmd.get_model("molecule")
cx = sum([atom.coord[0] for atom in model.atom])/len(model.atom)
cy = sum([atom.coord[1] for atom in model.atom])/len(model.atom)
cz = sum([atom.coord[2] for atom in model.atom])/len(model.atom)
middle = (cx, cy, cz)

# Get the min and max coordinates for x, y, and z
minX, minY, minZ = cmd.get_extent("molecule")[0]
maxX, maxY, maxZ = cmd.get_extent("molecule")[1]

# find overall absolute max
value = max(abs(maxX - cx), abs(maxY - cy), abs(maxZ - cz), abs(minX - cx), abs(minY -cy), abs(minZ -cz))


# Adjust to middle the cube around the molecule's middle of geometry
corners = [
    (cx - value, cy - value, cz - value),
    (cx - value, cy - value, cz + value),
    (cx - value, cy + value, cz - value),
    (cx - value, cy + value, cz + value),
    (cx + value, cy - value, cz - value),
    (cx + value, cy - value, cz + value),
    (cx + value, cy + value, cz - value),
    (cx + value, cy + value, cz + value)
]

for i, corner in enumerate(corners):
    cmd.pseudoatom(f"corner{i+1}", pos=corner)

# hide pseudoatoms
cmd.hide("everything", "corner*")

# List of tuples representing pairs of corners to connect
edges = [
    (1, 2), (1, 3), (1, 5), (2, 4), (2, 6),
    (3, 4), (3, 7), (4, 8), (5, 6), (5, 7),
    (6, 8), (7, 8)
]

for edge in edges:
    cmd.distance(f"edge{edge[0]}_{edge[1]}", f"corner{edge[0]}", f"corner{edge[1]}")

# Customize the appearance
cmd.hide("labels")
cmd.set("dash_width", 20.0)
cmd.set("dash_color", "yellow", "edge*")


# The cube should be devided into 10x10x10 sections
# Create the grid lines
grid_size = 10
for i in range(grid_size + 1):
    # Create lines parallel to the x-axis
    cmd.pseudoatom(f"x_line_{i}", pos=(cx - value + i * (2 * value / grid_size), cy - value, cz - value))
    cmd.pseudoatom(f"x_line_{i + grid_size + 1}", pos=(cx - value + i * (2 * value / grid_size), cy + value, cz - value))
    cmd.pseudoatom(f"x_line_{i + 2 * (grid_size + 1)}", pos=(cx - value + i * (2 * value / grid_size), cy - value, cz + value))
    cmd.pseudoatom(f"x_line_{i + 3 * (grid_size + 1)}", pos=(cx - value + i * (2 * value / grid_size), cy + value, cz + value))

    # Create lines parallel to the y-axis
    cmd.pseudoatom(f"y_line_{i}", pos=(cx - value, cy - value + i * (2 * value / grid_size), cz - value))
    cmd.pseudoatom(f"y_line_{i + grid_size + 1}", pos=(cx + value, cy - value + i * (2 * value / grid_size), cz - value))
    cmd.pseudoatom(f"y_line_{i + 2 * (grid_size + 1)}", pos=(cx - value, cy - value + i * (2 * value / grid_size), cz + value))
    cmd.pseudoatom(f"y_line_{i + 3 * (grid_size + 1)}", pos=(cx + value, cy - value + i * (2 * value / grid_size), cz + value))

    # Create lines parallel to the z-axis
    cmd.pseudoatom(f"z_line_{i}", pos=(cx - value, cy - value, cz - value + i * (2 * value / grid_size)))
    cmd.pseudoatom(f"z_line_{i + grid_size + 1}", pos=(cx + value, cy - value, cz - value + i * (2 * value / grid_size)))
    cmd.pseudoatom(f"z_line_{i + 2 * (grid_size + 1)}", pos=(cx - value, cy + value, cz - value + i * (2 * value / grid_size)))
    cmd.pseudoatom(f"z_line_{i + 3 * (grid_size + 1)}", pos=(cx + value, cy + value, cz - value + i * (2 * value / grid_size)))
    
# Hide the pseudoatoms used for grid lines
cmd.hide("everything", "x_line_*")
cmd.hide("everything", "y_line_*")
cmd.hide("everything", "z_line_*")

# Show the grid lines
cmd.show("lines", "x_line_*")
cmd.show("lines", "y_line_*")
cmd.show("lines", "z_line_*") 

step = 2*value/grid_size
# Generate & bond the X‑parallel grid‑lines
for jy in range(grid_size+1):
    for kz in range(grid_size+1):
        pts = []
        for i in range(grid_size+1):
            x = cx - value + i*step
            y = cy - value + jy*step
            z = cz - value + kz*step
            name = f"grid_x_{jy}_{kz}_{i}"
            cmd.pseudoatom(name, pos=(x,y,z))
            pts.append(name)
        # bond each to the next
        for a,b in zip(pts, pts[1:]):
            cmd.bond(a, b)
            
# Generate & bond the Y‑parallel grid‑lines
for ix in range(grid_size+1):
    for kz in range(grid_size+1):
        pts = []
        for j in range(grid_size+1):
            x = cx - value + ix*step
            y = cy - value + j*step
            z = cz - value + kz*step
            name = f"grid_y_{ix}_{kz}_{j}"
            cmd.pseudoatom(name, pos=(x,y,z))
            pts.append(name)
        # bond each to the next
        for a,b in zip(pts, pts[1:]):
            cmd.bond(a, b)
            
# Generate & bond the Z‑parallel grid‑lines
for ix in range(grid_size+1):
    for jy in range(grid_size+1):
        pts = []
        for k in range(grid_size+1):
            x = cx - value + ix*step
            y = cy - value + jy*step
            z = cz - value + k*step
            name = f"grid_z_{ix}_{jy}_{k}"
            cmd.pseudoatom(name, pos=(x,y,z))
            pts.append(name)
        # bond each to the next
        for a,b in zip(pts, pts[1:]):
            cmd.bond(a, b)
            
# Now hide all the spheres and show the bonds as lines:
cmd.hide("spheres", "grid_*")
cmd.show("lines", "grid_*")
cmd.set("line_width", 1.0)
cmd.set("line_color", "yellow")