import sys
from pylab import *
import gdsCAD as gds
from astropy.modeling.models import Sersic2D, Gaussian2D

#*********************SUBROUTINES************************

def GaussianGalaxy1(cell_name, center, sigmax, sigmay, angle, num_squares, square_size=1.0, layer=1): 
    # We're going to create a brightness profile by randomly dithering
    # num_squares squares of size square_size
    # Create a Cell and add the box
    cell=gds.core.Cell(cell_name)

    for n in range(num_squares):
        # Use Box-Muller algorithm to generate two Gaussian random numbers
        rsq = 1000.0 
        while (rsq >= 1.0 or rsq == 0.0):
            v1 = 2.0 * rand() - 1.0
            v2 = 2.0 * rand() - 1.0
            rsq = v1 * v1 + v2 * v2
        fac = sqrt(-2.0 * log(rsq) / rsq)
        x0 = sigmax * v1 * fac
        y0 = sigmay * v2 * fac
        x = center[0] + x0 * cos(angle) + y0 * sin(angle)
        y = center[1] - x0 * sin(angle) + y0 * cos(angle)
        x = round(x / square_size) * square_size
        y = round(y / square_size) * square_size    
        pixel=gds.shapes.Rectangle((x-square_size/2.0, y-square_size/2.0), (x+square_size/2.0, y+square_size/2.0), layer=layer)
        cell.add(pixel)
    return cell


def GaussianGalaxy(cell_name, center, sigmax, sigmay, angle, num_squares, square_size=1.0, layer=1): 
    # We're going to create a brightness profile by randomly dithering
    # num_squares squares of size square_size
    # Create a Cell and add the box
    cell=gds.core.Cell(cell_name)

    for n in range(num_squares):
        # Here we use rejection sampling to generate the 2D Gaussian Profile
        # This is slower, but more general.
        Reject = True
        while Reject:
            x0 = - 3.0 * sigmax + 6.0 * sigmax * rand()
            y0 = - 3.0 * sigmay + 6.0 * sigmay * rand()            
            if Gaussian2D.evaluate(x0, y0, 1.0, 0.0, 0.0, sigmax, sigmay, 0.0) > rand():
                Reject = False
        x = center[0] + x0 * cos(angle) + y0 * sin(angle)
        y = center[1] - x0 * sin(angle) + y0 * cos(angle)
        x = round(x / square_size) * square_size
        y = round(y / square_size) * square_size    
        pixel=gds.shapes.Rectangle((x-square_size/2.0, y-square_size/2.0), (x+square_size/2.0, y+square_size/2.0), layer=layer)
        cell.add(pixel)
    return cell

def SersicGalaxy(cell_name, center, r_eff, sersic_n, ellip, angle, num_squares, square_size=1.0, layer=1): 
    # We're going to create a brightness profile by randomly dithering
    # num_squares squares of size square_size
    # Create a Cell and add the box
    cell=gds.core.Cell(cell_name)

    for n in range(num_squares):
        # Here we use rejection sampling to generate the 2D Sersic profile
        Reject = True
        while Reject:
            x0 = - 6.0 * r_eff + 12.0 * r_eff * rand()
            y0 = - 6.0 * r_eff + 12.0 * r_eff * rand()            
            if Sersic2D.evaluate(x0, y0, 1.0, r_eff, sersic_n, 0.0, 0.0, ellip, 0.0) > rand():
                Reject = False
        x = center[0] + x0 * cos(angle) + y0 * sin(angle)
        y = center[1] - x0 * sin(angle) + y0 * cos(angle)
        x = round(x / square_size) * square_size
        y = round(y / square_size) * square_size    
        pixel=gds.shapes.Rectangle((x-square_size/2.0, y-square_size/2.0), (x+square_size/2.0, y+square_size/2.0), layer=layer)
        cell.add(pixel)
    return cell

def Star(cell_name, center, size, layer=1):
    # Stars are round images of variable size 
    cell=gds.core.Cell(cell_name)
    spot=gds.shapes.Disk(center, size, layer=layer)
    cell.add(spot)
    return cell

#*********************MAIN PROGRAM************************
# Create a layout and add the cell
layout = gds.core.Layout('LIBRARY')
topcell=gds.core.Cell('topcell')
cell1 = GaussianGalaxy1('cell1', (200.0,0.0), 20.0, 10.0, 45.0, 800)
topcell.add(cell1)
cell2 = GaussianGalaxy('cell2', (400.0,0.0), 10.0, 20.0, 0.0, 500)
topcell.add(cell2)
cell3 = SersicGalaxy('cell3', (-100.0,0.0), 20.0, 2.0, 0.3, 30.0, 1500)
topcell.add(cell3)
cell4 = SersicGalaxy('cell4', (100.0,300.0), 20.0, 10.0, 0.3, 120.0, 800)
topcell.add(cell4)

cell5 = Star('cell5', (120.0,40.0), 6.0)
topcell.add(cell5)

cell6 = Star('cell6', (70.0,-40.0), 3.0)
topcell.add(cell6)

cell7 = Star('cell7', (-60.0,150.0), 4.0)
topcell.add(cell7)

layout.add(topcell)
# Save the layout
layout.save('scene_test.gds')

