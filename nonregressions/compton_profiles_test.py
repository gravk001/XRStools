from xrstools import xrs_ComptonProfiles as CP
import numpy as np

from pylab import *
ion()

filename = '/home/christop/Dropbox/XRSTools/xrstools/data/ComptonProfiles.dat'

cp = CP.AtomProfile('Si',filename)


