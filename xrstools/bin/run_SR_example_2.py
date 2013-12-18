from xrstools import xrs_read, theory, extraction
from pylab import *

t = xrs_read.read_id20('/home/christoph/data/ch3914/orig/raman',energycolumn='energy',monitorcolumn='kap4dio')
t.make_posscan_image(592,'sty','/home/christoph/data/ch3914/scratch/test.dat')
