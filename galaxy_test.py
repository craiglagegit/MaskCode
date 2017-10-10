import sys
from pylab import *
from gdsCAD import *

# Create a Cell and add the box
cell=core.Cell('TOP')

LL = (0.0,0.0)
# We're going to create a brightness profile by randomly dithering these pixels
pixel_size = 1.0

sigmax = 20.0
sigmay = 10.0
layer = 1

for n in range(200):
    #Use Box-Muller algorithm to generate two Gaussian random numbers
    rsq = 1000.0 
    while (rsq >= 1.0 or rsq == 0.0):
        v1 = 2.0 * rand() - 1.0
        v2 = 2.0 * rand() - 1.0
        rsq = v1 * v1 + v2 * v2
    fac = sqrt(-2.0 * log(rsq) / rsq)
    x = sigmax * v1 * fac
    y = sigmay * v2 * fac
    x = round(x / pixel_size) * pixel_size
    y = round(y / pixel_size) * pixel_size    
    pixel=shapes.Rectangle((x-pixel_size/2.0, y-pixel_size/2.0), (x+pixel_size/2.0, y+pixel_size/2.0), layer=layer)
    cell.add(pixel)

# Create a layout and add the cell
layout = core.Layout('LIBRARY')
layout.add(cell)

# Save the layout
layout.save('galaxy_test.gds')

