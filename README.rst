XRStools
===========================================

Main development website: https://

|Build Status| |Appveyor Status|

Bla Bla

References
----------

* Ch.J. Sahle, A. Mirone, J. Niskanen, J. Inkinen, M. Krisch, and S. Huotari:Journal of Synchrotron Radiation 22, No. 2 (2015): 400-409.

Installation
------------

With Conda
........

Conda installation is buggy but works on Windows/Linux & Macintosh OS. Tested with conda version: 4.7.12

Step 0 - Install conda -- [miniconda comes with less bulk](https://conda.io/en/latest/miniconda.html)
Step 1 - Install git and clone the XRStools repository:

    conda install git

    git clone git@github.com:LJRH/XRStools.git

    git checkout LERIX-working #this step for APS LERIX support

Step 2 - Change directory to the cloned XRStools repo and create a new conda virtual environment called "xrstools":

    cd <...>/XRStools/

    conda create -n xrstools python=2.7 --file requirements.txt

Step 3 - Install other required libraries including praxes's version of pymca (newer versions don't work yet):

    conda activate xrstools

    conda install -yc praxes pymca

    conda install -c conda-forge silx

Step 4 - Intitiate XRStools (this might generate errors, ignore these.)

    python setup.py install [--user]

Step 5 - Open iPython and check that you are able to load important modules:

    ipython
    
    from XRStools import xrs_read xrs_utilities xrs_extraction

With PIP
........

As most Python packages, pyFAI is available via PIP::

   pip install xrstools [--user]

Provide the *--user* to perform an installation local to your user.
Under UNIX, you may have to run the command via *sudo* to gain root access an
perform a system wide installation.



Documentation
-------------

Documentation can be build using this command and Sphinx (installed on your computer)::

    python setup.py build build_doc
