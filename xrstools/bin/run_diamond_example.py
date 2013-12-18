import numpy as np
from scipy import interpolate
from xrstools import xrs_read, theory, extraction
from pylab import *
ion()

# try loading old Si data
counters=[1, 2, 3, 4, 5, 6, 7, 8, 9]
data = np.loadtxt('xrstools/things/diamond_data_E016keV_tth123.dat')

# remove some glitches in the data
glitchindex1 = np.where(data[:,0]==1147.0999999999999)[0][0]
data = np.delete(data, glitchindex1, 0)
glitchindex2 = np.where(data[:,0]==397.10000000000002)[0][0]
data = np.delete(data, glitchindex2, 0)

meantth = 129
tth = np.array([meantth-13.0, meantth-6.5, meantth-6.5, meantth+0.0, meantth+0.0, meantth+6.5, meantth+6.5, meantth+13.0, meantth+13.0])

# initiate data
eloss   = data[:,0]
signals = data[:,1:]
errors  = np.sqrt(np.absolute(data[:,1:]))
E0      = 15.8

# initiate id20 instance (from a random Spec-file, since the init checks if the specified file actually exists)
dia = xrs_read.read_id16('xrstools/things/licl_test_files/raman')

dia.eloss   = eloss
dia.signals = np.absolute(signals)
dia.errors  = errors
dia.E0      = E0
dia.tth     = tth

# initiate cprofiles instance
hf = theory.HFspectrum(dia,['C'],[1.0],correctasym=[[0.0]])

# initiate instance of extraction class for background removal:
#reload(extraction)
extr = extraction.extraction(dia,hf,prenormrange=[282,np.inf])

extr.extractval_test(range(9),linrange1=[340,470],linrange2=[1500,2500])

extr.getallvalprof(7,smoothgval=50.0,stoploop=False)
extr.remvalenceprof_test(range(9),eoffset=20.0)




