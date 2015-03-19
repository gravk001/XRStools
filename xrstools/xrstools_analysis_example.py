from xrstools import xrs_read, theory, extraction, roiSelectionWidget, xrs_rois
from pylab import *
ion()
import numpy as np
from scipy import interpolate, signal, integrate, constants, optimize, ndimage

from  PyQt4 import Qt, QtCore


lowq  = range(12)
lowq.extend(range(36,60))
medq  = range(12,24)
highq = range(24,36)
highq.extend(range(60,72))

##############################################################################
# H2o example
##############################################################################

repertorio = "/scisoft/users/mirone/WORKS/Christoph/for_alessandro"
# repertorio = "/home/alex/WORKS/Christoph/for_alessandro"
# repertorio = "/home/christoph/data/ihr_feb15/"

h2o = xrs_read.read_id20(repertorio + '/hydra',monitorcolumn='kapraman')

image4roi =  h2o.SumDirect( [623] )

app=Qt.QApplication([])
w4r = roiSelectionWidget.mainwindow()
w4r.showImage( image4roi , xrs_rois.get_geo_informations(image4roi.shape) )
w4r.show()
app.exec_()

masks = w4r.getMasks()
print masks

print image4roi.shape

import pickle
f = open(repertorio +'/rois/as_2m_72_big.txt','rb')
roi_obj = pickle.load(f)

print roi_obj.red_rois

# f = open(  repertorio +'/rois/as_2m_72_Dict.pick','wb'    )
# pickle.dump(roi_obj.red_rois, f)
# f.close()


raise "OK"


h2o.load_rois(repertorio +'/rois/as_2m_72_big.txt')
print "OK "
# h2o.getautorois_eachdet([623],thresfrac=10)
# h2o.get_zoom_rois([623])
# h2o.save_rois('/home/christoph/data/ihr_feb15/rois/h2o_72_big.txt')
h2o.loadelasticdirect([623])


h2o.loadloopdirect([625,629,633,637,641],1)
h2o.loadlongdirect(624)
print " OK " 
# h2o.getrawdata()
h2o.getspectrum()
h2o.geteloss()
h2o.gettths(rvd=-41,rvu=85,rvb=121.8,rhl=41.0,rhr=41.0,rhb=121.8,order=[0,1,2,3,4,5])


# O K edge
hf   = theory.HFspectrum(h2o,['O'],[1.0],correctasym=[[0.0,0.0,0.0]])
extr1 = extraction.extraction(h2o,hf)

extr1.analyzerAverage(lowq,errorweighing=False)
#extr1.removePearsonAv2([480.0,533.0],[550.0,555.0],scale=3.6,hfcoreshift=0.0)
extr1.removeLinearAv([520.0,532.0],scale=0.052)
extr1.savetxtsqwav(repertorio+'/as_abs/low_q/h2o_OK.txt', emin=500.0, emax=700.0)

x= raw_input()

