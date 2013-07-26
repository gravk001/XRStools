import numpy as np
from scipy import interpolate
from xrstools import xrs_read, theory, extraction
from pylab import *
ion()

# try loading old Si data
counters=[1, 2, 3, 4, 5, 6, 7, 8, 9]
data = np.loadtxt('xrstools/things/diamond_data_E016keV_tth123.dat')

tth = np.array([123-13.0, 123-6.5, 123-6.5, 123+0.0, 123+0.0, 123+6.5, 123+6.5, 123+13.0, 123+13.0])+8.0

# initiate data
eloss   = data[:,0]
signals = data[:,1:]
errors  = np.sqrt(np.absolute(data[:,1:]))
E0      = 16

# initiate id20 instance
dia = xrs_read.read_id20('xrstools/things/licl_test_files/raman')

dia.eloss   = eloss
dia.signals = np.absolute(signals)
dia.errors  = errors
dia.E0      = E0
dia.tth     = tth

# initiate cprofiles instance
hf = theory.HFspectrum(dia,['C'],[1.0],correctasym=[[0.0]])

# initiate instance of extraction class for background removal:
extr = extraction.extraction(dia,hf,prenormrange=[282,np.inf])

plot(extr.eloss,extr.J,extr.eloss,extr.signals)

