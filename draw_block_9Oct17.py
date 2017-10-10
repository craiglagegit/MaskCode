import sys, datetime
from pylab import *
from random import shuffle
import gdsCAD as gds
from astropy.modeling.models import Sersic2D, Gaussian2D

#*********************SUBROUTINES************************

class BoundingCell:
    def __init__(self, cell_name, LL, UR, square_size=1.0, PixelSizeX=10.0, PixelSizeY=10.0, line_layer=2):
        # This defines an array which contains CCD pixels,
        # and gds 'squares'.  Each grid point in this array
        # represents one gds square
        self.dx=square_size
        self.PixelSizeX = PixelSizeX
        self.PixelSizeY = PixelSizeY
        self.dy=self.dx
        self.xmin = round(LL[0] / square_size) * square_size
        self.ymin = round(LL[1] / square_size) * square_size        
        self.xmax = round(UR[0] / square_size) * square_size
        self.ymax = round(UR[1] / square_size) * square_size        
        self.nx=int(round((self.xmax - self.xmin) / self.dx))
        self.ny=int(round((self.ymax - self.ymin) / self.dy))        
        self.x=linspace(self.xmin+self.dx/2,self.xmax-self.dx/2,self.nx)
        self.y=linspace(self.ymin+self.dy/2,self.ymax-self.dy/2,self.ny)
        self.data=zeros([self.nx,self.ny], dtype=bool)
        self.cell=gds.core.Cell(cell_name)

        # We show the CCD pixels with a layer that won't print
        for i in range(self.nx):
            if int(round((self.x[i] - self.dx))) % int(round(PixelSizeX)) == 0:
                line = gds.shapes.Box((self.x[i]-self.dx/2,self.ymin),(self.x[i]-self.dx/2,self.ymax), square_size/2.0, layer=line_layer)
                self.cell.add(line)
        line = gds.shapes.Box((self.x[0]-self.dx/2,self.ymin),(self.x[0]-self.dx/2,self.ymax), square_size/2.0, layer=line_layer)
        self.cell.add(line)
        line = gds.shapes.Box((self.x[self.nx-1]+self.dx/2,self.ymin),(self.x[self.nx-1]+self.dx/2,self.ymax), square_size/2.0, layer=line_layer)
        self.cell.add(line)

        for j in range(self.ny):
            if int(round((self.y[j] - self.dy))) % int(round(PixelSizeY)) == 0:
                line = gds.shapes.Box((self.xmin,self.y[j]-self.dy/2),(self.xmax,self.y[j]-self.dy/2), square_size/2.0, layer=line_layer)
                self.cell.add(line)
        line = gds.shapes.Box((self.xmin,self.y[0]-self.dy/2),(self.xmax,self.y[0]-self.dy/2), square_size/2.0, layer=line_layer)
        self.cell.add(line)
        line = gds.shapes.Box((self.xmin,self.y[self.ny-1]+self.dy/2),(self.xmax,self.y[self.ny-1]+self.dy/2), square_size/2.0, layer=line_layer)
        self.cell.add(line)

        
def SersicGalaxy(boundingcell, cell_name, center, r_eff, sersic_n, ellip, angle, m, square_size=1.0, layer=1): 
    # We're going to create a brightness profile by randomly dithering
    # num_squares squares of size square_size
    # magnitudes are calculated assuming a 10" exposure at 100% light intensity
    # Create a Cell and add the box
    text_step = 2.0 * boundingcell.PixelSizeX 
    num_squares = int(pow(10.0, 0.4 * (27.9 - m)))    
    cell=gds.core.Cell(cell_name)

    nxmin = max(0, int((center[0] - 10.0 * r_eff - boundingcell.xmin) / boundingcell.dx))
    nxmax = min(boundingcell.nx-1, int((center[0] + 10.0 * r_eff - boundingcell.xmin) / boundingcell.dx))    
    nymin = max(0, int((center[1] - 10.0 * r_eff - boundingcell.ymin) / boundingcell.dy))
    nymax = min(boundingcell.ny-1, int((center[1] + 10.0 * r_eff - boundingcell.ymin) / boundingcell.dy))    
    
    for n in range(num_squares):
        # Here we use rejection sampling to generate the 2D Sersic profile
        Reject = True
        while Reject:
            AlreadyFilled = True
            while AlreadyFilled:
                nx = randint(nxmin, nxmax)
                ny = randint(nymin, nymax)            
                AlreadyFilled = boundingcell.data[nx,ny]
            x = boundingcell.x[nx]
            y = boundingcell.y[ny]
            r = sqrt((x-center[0])**2 + (y-center[1])**2)
            #print n, center, nx, ny, x, y, r, Sersic2D.evaluate(x, y, 1.0, r_eff, sersic_n, center[0], center[1], ellip, angle)
            if Sersic2D.evaluate(x, y, 1.0, r_eff, sersic_n, center[0], center[1], ellip, angle) > rand():
                Reject = False
        boundingcell.data[nx,ny] = True
        pixel=gds.shapes.Rectangle((x-boundingcell.dx/2.0, y-boundingcell.dy/2.0), (x+boundingcell.dx/2.0, y+boundingcell.dy/2.0), layer=layer)
        cell.add(pixel)

    # Now we annotate it with non-printing text.
    name = gds.core.Text('Galaxy', (center[0] - 2.0 * text_step, center[1] - 2.0 * text_step), layer=3, magnification=0.002)
    cell.add(name)
    mag = gds.core.Text('Mag = %.2f'%m, (center[0] - 2.0 * text_step, center[1] - 3.0 * text_step), layer=3, magnification=0.002)
    cell.add(mag)
    reff = gds.core.Text('r_eff = %.2f arcseconds'%(r_eff / boundingcell.PixelSizeX * 0.2), (center[0] - 2.0 * text_step, center[1] - 4.0 * text_step), layer=3, magnification=0.002)
    cell.add(reff)
    sersicn = gds.core.Text('n = %.2f'%sersic_n, (center[0] - 2.0 * text_step, center[1] - 5.0 * text_step), layer=3, magnification=0.002)
    cell.add(sersicn)
    ell = gds.core.Text('ellip = %.2f'%ellip, (center[0] - 2.0 * text_step, center[1] - 6.0 * text_step), layer=3, magnification=0.002)
    cell.add(ell)
    rot = gds.core.Text('angle = %.2f degrees'%(angle*180.0/pi), (center[0] - 2.0 * text_step, center[1] - 7.0 * text_step), layer=3, magnification=0.002)
    cell.add(rot)
    #print "Galaxy, m = %f, n = %d"%(m, num_squares)
    return cell

def Star(boundingcell, cell_name, center, m, layer=1):
    # Stars are approximately round images
    # We use the same algorithm to calculate how many squares per star
    # magnitudes are calculated assuming a 10" exposure at 100% light intensity    
    text_step = 2.0 * boundingcell.PixelSizeX 
    num_squares = int(pow(10.0, 0.4 * (27.9 - m)))    
    cell=gds.core.Cell(cell_name)
    r0 = boundingcell.dx/2.0
    n = 0
    # ramp up r until we get the right number of squares
    while n < num_squares:
        r0 += 0.01
        nxmin = max(0, int((center[0] - r0 - boundingcell.xmin) / boundingcell.dx))
        nxmax = min( boundingcell.nx-1, int((center[0] + r0 - boundingcell.xmin) / boundingcell.dx))
        nymin = max(0, int((center[1] - r0 - boundingcell.ymin) / boundingcell.dy))
        nymax = min(boundingcell.ny-1, int((center[1] + r0 - boundingcell.ymin) / boundingcell.dy))    
        n = 0
        for i in range(nxmin, nxmax+1):
            for j in range(nymin, nymax+1):
                x = boundingcell.x[i]
                y = boundingcell.y[j]
                r = sqrt((x-center[0])**2 + (y-center[1])**2)
                if r < r0:
                    n += 1
    # Now populate the actual cell with the desired number of squares, as closely as we can.
    n = 0
    for i in range(nxmin, nxmax+1):
        for j in range(nymin, nymax+1):
            x = boundingcell.x[i]
            y = boundingcell.y[j]
            r = sqrt((x-center[0])**2 + (y-center[1])**2)
            if r < r0:
                pixel=gds.shapes.Rectangle((x-boundingcell.dx/2.0, y-boundingcell.dy/2.0), (x+boundingcell.dx/2.0, y+boundingcell.dy/2.0), layer=layer)
                cell.add(pixel)
                n += 1

    # Now adjust the magnitude to reflect the actual number of squares
    m = 27.9 - 2.5 * log10(n)
    # Now we annotate it with non-printing text.                    
    name = gds.core.Text('Star', (center[0] - 2.0 * text_step, center[1] - 2.0 * text_step), layer=3, magnification=0.002)
    cell.add(name)
    mag = gds.core.Text('Mag = %.2f'%m, (center[0] - 2.0 * text_step, center[1] - 3.0 * text_step), layer=3, magnification=0.002)
    cell.add(mag)
    #print "Star, m = %f, n = %d, actual squares = %d"%(m, num_squares, n)    
    return cell
    

#*********************MAIN PROGRAM************************

# This defines a rough distribution of star magnitudes
star_mags = list(linspace(20.0,21.0,5))\
          + list(linspace(21.0,22.0,10))\
          + list(linspace(22.0,23.0,20))\
          + list(linspace(23.0,24.0,40))\
          + list(linspace(24.0,25.0,80))          

# First, define the object list
Nx = 30
Ny = 30
NumStars = 300
NumGals = Nx * Ny - NumStars
ObjectSpacing = 400.0 # in microns
# Create the bounding cell
boundingcell = BoundingCell('boundingcell', (0.0,0.0), (Nx * ObjectSpacing, Ny * ObjectSpacing))
objects = []
shuffled_index = range(NumStars + NumGals) # shuffling the indices randomizes the order
shuffle(shuffled_index)

for i in range(NumStars):
    star_mag = star_mags[int((len(star_mags)-1) * rand())]
    objects.append(['Star', star_mag])

for i in range(NumGals):
    # This defines the galaxy characteristics we are drawing
    gal_mag = 20.0 + rand() * 5.0
    r_eff = 15.0 + rand() * 10.0
    sersic_n = 1.0 + rand() * 3.0
    ellip = rand() * 0.8
    angle = rand() * pi
    objects.append(['Galaxy',r_eff,sersic_n,ellip,angle,gal_mag])

# Create a layout and add the cell
layout = gds.core.Layout('LIBRARY')

# Now create the objects from the list
now = datetime.datetime.now()
timestamp = "%4d%02d%02d%02d%02d%02d"%(now.year,now.month,now.day,now.hour,now.minute,now.second)
catfile = open("star_galaxy_catalog_%s.txt"%timestamp,'w')
headerline = "ID\tS_ID\tType\tX\t\tY\t\tMag\tR_eff\tN\tEllip\tAngle(degrees)\t\n"
catfile.write(headerline)
StarIndex = 0
GalIndex = 0
for i in range(Nx):
    x = ObjectSpacing * (i + 0.5)
    for j in range(Ny):
        y = ObjectSpacing * (j + 0.5)
        index = i * Ny + j
        obj = objects[shuffled_index[index]] # This randomizes the object placement
        if obj[0] == 'Star':
            cell = Star(boundingcell, 'starcell_%d'%StarIndex, (x,y), obj[1])
            boundingcell.cell.add(cell)
            catline = "%d\t%d\tStar\t%.1f\t\t%.1f\t\t%.3f\n"%(index,shuffled_index[index],x,y,obj[1])
            StarIndex += 1
        else:
            [r_eff,sersic_n,ellip,angle,gal_mag] = obj[1:] 
            cell = SersicGalaxy(boundingcell, 'galcell_%d'%GalIndex, (x,y), r_eff, sersic_n, ellip, angle, gal_mag)
            boundingcell.cell.add(cell)
            catline = "%d\t%d\tGalaxy\t%.1f\t\t%.1f\t\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t\n"%(index,shuffled_index[index],x,y,gal_mag,r_eff,sersic_n,ellip,angle*180.0/pi)
            GalIndex += 1
        catfile.write(catline)
    print "%d stars and %d galaxies done"%(StarIndex,GalIndex)
catfile.close()
# Create the top cell, within which we will place multiple copies of the bounding cell
topcell=gds.core.Cell('topcell')

# Now step this cell in a 3x3 array
Mx = 3
My = 3
for i in range(Mx):
    x= Nx * ObjectSpacing * i
    for j in range(My):
        y = Ny * ObjectSpacing * j
        topcell.add(boundingcell.cell, origin=(x,y))
layout.add(topcell)
# Save the layout
layout.save('star_galaxy_mask_%s.gds'%timestamp)

