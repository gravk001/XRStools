import numpy as np
from scipy import interpolate
from xrstools import xrs_read, theory, extraction
from pylab import *

# try loading old Si data
counters=[3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19]
data = []
for counter in counters:
	try:
		nr = '%02d' % counter
		fn = '/home/christoph/data/APS_old/Si/fig_raw_si_' + nr + '.dat'
		data.append(np.loadtxt(fn))
	except:
		print 'Det ' + nr + ' not found!'

tth = np.array([27, 36, 54, 63, 72, 81, 90, 99, 108, 117, 126, 135, 144, 153, 171])+8

# initiate data
eloss   = data[0][:,0]
signals = np.zeros((len(eloss),len(counters)))
errors  = np.zeros(np.shape(signals))
E0      = 9.893 # np.zeros((len(counters),1))+9.893

# spline onto first detector
for n in range(len(counters)):
	f = interpolate.interp1d(data[n][:,0],data[n][:,1], bounds_error=False, fill_value=0.0)
	signals[:,n] = f(eloss)
	f = interpolate.interp1d(data[n][:,0],data[n][:,2], bounds_error=False, fill_value=0.0)
	errors[:,n]  = f(eloss)

# initiate id20 instance
si = xrs_read.read_id20('xrstools/things/licl_test_files/raman')

si.eloss = eloss
si.signals = signals
si.errors  = errors
si.E0 = E0
si.tth = tth

#pylab.plot(si.eloss,si.signals,'.')
#pylab.show()

# initiate cprofiles instance
hf = theory.HFspectrum(si,['Si'],[1.0],correctasym=[[1.4]])

# initiate instance of extraction class for background removal:
extr = extraction.extraction(si,hf)

extr.removecoreppearson(0,[35,99],[310,1050])

# high q valence extraction
extr.extractval(14,linrange1=[40,95],linrange2=[850,1200],mirror=True)
extr.getallvalprof(14,smoothgval=5.0,stoploop=False)
extr.remvalenceprof(14)

# quick valence removal
extr.remquickval(13,[700.0,1500.0],[98.0,115.0],10.0)

pylab.plot(extr.eloss,extr.signals[:,0],extr.eloss,extr.J[:,0],extr.eloss,extr.C[:,0])
pylab.show()



