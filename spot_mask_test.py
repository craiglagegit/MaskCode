from gdsCAD import *

# Create a Cell and add the box
cell=core.Cell('TOP')

LL = (0.0,0.0)
spot_spacing = 200.0
spot_radius = 15.0
Num_rows = 5
Num_cols = 5
layer = 1
# Create an array of spots

for i in range(Num_cols):
    for j in range(Num_rows):
        spot_origin = (LL[0] + i * spot_spacing, LL[1] + j * spot_spacing) 
        spot=shapes.Disk(spot_origin, spot_radius, layer=layer)
        cell.add(spot)

# Create a layout and add the cell
layout = core.Layout('LIBRARY')
layout.add(cell)

# Save the layout
layout.save('spot_test.gds')

