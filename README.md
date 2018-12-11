# MaskCode - Craig Lage - UC Davis - 10-Oct-17
Working with the LSST Optical Simulator at  UC Davis, we desired to build
a set of masks containing simulated stars and galaxies that can be
projected onto the CCDs we are testing.  I wrote some Python code that placed
sub resolution spots in such a way as to simulate galaxy light profiles.
The code generates .gds files that can be sent to the mask vendor.
I am posting it here in case others find it useful.  The files spot_mask_test.py,
galaxy_test.py, draw_scene.py, and draw_block.py are tests of increasing
complexity I used to develop the code.

The file draw_block_9Oct17.py is the actual code I used to draw a mask containing
2700 simulated stars and 5400 simulated galaxies drawn from distributions.  The
galaxies have Sersic light profiles of various size, brightness, Sersic index,
ellipticity, and orientation.  A catalog of the objects is generated concurrently
with the .gds file.  We will be ordering a mask of this type soon.

Acknowledgements:
This research was supported by DOE grant DE-SC0009999 and NSF/AURA grant N56981C.