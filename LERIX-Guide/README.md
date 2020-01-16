Author: L.Higgins
Date: 16 Jan 20

-------------------------------------------------------------------------------
Introduction
-------------------------------------------------------------------------------
LERIX (Low Energy Resolution Inelastic X-ray) spectrometer is the instrument in operation at [beamline 20-ID](https://www.aps.anl.gov/Spectroscopy/Beamlines/20-ID) of the Advanced Photon Source (APS). The data collected from X-ray Raman scattering spectroscopy (XRSS) experiments at the APS has a very different format (Labview prodcing ASCII text docs) to the more data-efficient ESRF "spec" format. Therefore, in order to use the feature-rich XRStools software, a small add-on (read_lerix) was made. This guide should help those who want to analyze their data using this tool. It is is development software, made by a very amateur coder in python - it will be buggy. Any help from those with more skill is welcome.

-------------------------------------------------------------------------------
Getting Started - an example
-------------------------------------------------------------------------------
Included in the LERIX-guide directory is some test data of a graphite pellet on a spinner from 20-ID. The experiment was performed using the Si(311) mono at 9890 eV with XRS data collected for the C K edge.

If XRStools has been installed properly then on opening an iPython command line the following modules should be able to be called:

```python
    from XRStools import xrs_read xrs_utilities xrs_extraction
    import numpy as np
    import matplotlib.pyplot as plt
    import os
```

If at this stage there are errors such as 'conda: module not found', try and install the missing modules. The next step would be to rebuild XRStools `python setup.py install`.

If successful you should see the following output:

```
>>>>>>>>  use_PyMca  True
>>>>>>>>  use_PyMca  True
/opt/conda/envs/christoph6/lib/python2.7/site-packages/PyMca5/PyMcaGui/plotting/MaskImageWidget.pyc
Could not load PyTango
>>>>>>>>  use_PyMca  True

############################# Welcome to XRStools #############################
# If you are using this software, please cite the following work:             #
# Ch.J. Sahle, A. Mirone, J. Niskanen, J. Inkinen, M. Krisch, and S. Huotari: #
# "Planning, performing, and analyzing X-ray Raman scattering experiments."   #
# Journal of Synchrotron Radiation 22, No. 2 (2015): 400-409.                 #
###############################################################################
```

At this point we set the directory of the test data and read the experiment:

```python
path = '../XRStools/LERIX-guide/graphite/' #change this to the test data
graphite = xrs_read.read_lerix(exp_dir=path)
```
