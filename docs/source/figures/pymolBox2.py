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
    cmd.pseudoatom(f"m_corner{i+1}", pos=corner)

# hide pseudoatoms
cmd.hide("everything", "m_corner*")

# List of tuples representing pairs of corners to connect
edges = [
    (1, 2), (1, 3), (1, 5), (2, 4), (2, 6),
    (3, 4), (3, 7), (4, 8), (5, 6), (5, 7),
    (6, 8), (7, 8)
]

for edge in edges:
    cmd.distance(f"m_edge{edge[0]}_{edge[1]}", f"m_corner{edge[0]}", f"m_corner{edge[1]}")

# Customize the appearance
cmd.hide("labels")
cmd.set("dash_width", 50.0, "m_edge*")
cmd.set("dash_color", "black", "m_edge*")


# 2) Grid parameters
grid_size = 5
step = 2 * value / grid_size

# 3) Make & store pseudoatoms
grid = [[[None]*(grid_size+1) for _ in range(grid_size+1)]
            for __ in range(grid_size+1)]
for i in range(grid_size+1):
    x = cx - value + i*step
    for j in range(grid_size+1):
        y = cy - value + j*step
        for k in range(grid_size+1):
            z = cz - value + k*step
            name = f"g_{i}_{j}_{k}"
            cmd.pseudoatom(name, pos=(x,y,z))
            grid[i][j][k] = name

# 4) Create dashed “distance” objects between neighbors
#    along X, Y and Z
for i in range(grid_size+1):
    for j in range(grid_size+1):
        for k in range(grid_size+1):
            this = grid[i][j][k]
            if i < grid_size:
                cmd.distance(f"d_x_{i}_{j}_{k}", this, grid[i+1][j][k])
            if j < grid_size:
                cmd.distance(f"d_y_{i}_{j}_{k}", this, grid[i][j+1][k])
            if k < grid_size:
                cmd.distance(f"d_z_{i}_{j}_{k}", this, grid[i][j][k+1])

# 5) Tidy up
cmd.hide("labels", "d_*")               # no distance labels
cmd.hide("spheres", "g_*")              # hide all pseudo‑atoms
cmd.set("dash_width", 0.1, "d_")              # adjust thickness
#cmd.set("dash_gap",   0.2)              # adjust dash spacing
cmd.set("dash_color", "red", "d_*")  # color all d_ distance objs

# 6) (Optional) raise the dash line depth bias 
#    so the grid is always visible in front of your molecule:
cmd.set("depth_cue", 0)


