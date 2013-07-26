from xrstools import xrs_read, theory, extraction
import pylab

licl = xrs_read.read_id20('xrstools/things/licl_test_files/raman')
licl.loadelastic(90)
licl.getautorois(90)
licl.loadloop([165,170],5)
licl.loadlong(157)

licl.getrawdata()
licl.getspectrum()
licl.geteloss()
licl.gettths(35)
licl.getqvals()

for n in range(len(licl.signals[0,:])):
	pylab.plot(licl.eloss,licl.signals[:,n])

pylab.axis([-50.0,620.0,0.0,0.005])
pylab.show()

hf = theory.HFspectrum(licl,['H2O'],[1.0],correctasym=[[0.0,0.0]])
extr = extraction.extraction(licl,hf)

extr.removeconstpcore(range(9),[527.0,534.5],[545.0,588.0],weights=[1,1],stoploop=False)
extr.averageqs(range(9))

pylab.plot(extr.eloss,extr.sqwav)
pylab.show()



