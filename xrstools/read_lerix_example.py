from XRStools import xrs_read, theory, extraction, xrs_utilities
import numpy as np
import matplotlib.pyplot as plt
path = "/home/jovyan/work/xrstools_testing/CNT_PP/"
myq = [3,6,9] #chosen analysers of 19

# Create read_lerix class at chosen path
cnt = xrs_read.read_lerix(path)

#read elastics with chosen Analysers, choosing analysers at this stage calculates E0 and fwhm, should be done before reading nixs
# changing analyzers='all' you can see that some analysers don't regiser enough counts to be used to calculate E0 or fwhm (min:200, max:5000 cps)
cnt.load_elastics(analyzers=myq)
# read all the NIXS and analysers
cnt.load_nixs(scans='all')
cnt.load_wides(scans='all')

# Shape of ndarray containing the averaged signals
cnt.signals.shape
# Can then average the signals from good analysers (myq)
good_ans = cnt.signals[:,myq].mean(axis=1)

# Example of the scans class which holds all the information about each scan.
cnt.scans
cnt.scans['wide0001'].signals[,:]
# motors currently holds some random info which might be useful later on.
cnt.scans['nixs0001'].motors
cnt.E0
cnt.tth = np.array(range(9,180,9))
cnt.update_cenom('all')
cnt.cenom
# Plotting the good analysers
import matplotlib.pyplot as plt
plt.plot(cnt.eloss, good_ans)

element = ['C']
!pwd
filename = '/home/jovyan/work/xrstools_testing/XRStools/data/ComptonProfiles.dat'
import os
os.path.isfile(filename)
cnt.tth = np.array([35,36,37])
test = xrs_utilities.makeprofile_comp('C')
hf   = theory.HFspectrum(cnt,['C'],[1.0])






cnt.signals[0,:].shape
cnt.cenom
