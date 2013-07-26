from xrstools import calctools
from pylab import *

tmao = calctools.erkale('/home/christoph/programs/erkale-runs/tmao-6M/clusters_4A/results/cluster_mol','/trans_four-1.50.dat',0,980,20,stepformat=3)

tmao.cut_rawspecs(500,600)
tmao.broaden_lin()
tmao.sum_specs()
tmao.plot_spec()

urea = calctools.erkale('/home/christoph/programs/erkale-runs/urea-6M/clusters_4A/results/cluster_mol','/trans_four-1.50.dat',0,980,20,stepformat=3)
urea.cut_rawspecs(500,600)
urea.broaden_lin()

ion()
plot(urea.energy,urea.sqw,tmao.energy,tmao.sqw)
