#!/usr/bin/python
# Filename: xrs_ComptonProfiles.py

#/*##########################################################################
#
# The XRStools software package for XRS spectroscopy
#
# Copyright (c) 2013-2014 European Synchrotron Radiation Facility
#
# This file is part of the XRStools XRS spectroscopy package developed at
# the ESRF by the DEC and Software group.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
__author__ = "Christoph J. Sahle - ESRF"
__contact__ = "christoph.sahle@esrf.fr"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"

import numpy as np
import xrs_utilities
import xrs_fileIO

from scipy import interpolate, integrate, constants

class SqwPredict:
    """Class to build a S(q,w) prediction based on HF Compton Profiles.

    Attributes:
    -----------
    sampleStr (list of strings): one string per compound (e.g. ['C','SiO2'])
    concentrations (list of floats): relative compositional weight for each compound
class ComptonProfile:
    """Class to construct and handle Hartree-Fock atomic Compton Profiles.

    Attributes:
    -----------
      element (string): Element symbol as in the periodic table.
      elementNr (int) : Number of the element as in the periodic table.
      shells (list)   :
      edges (list)    :
      C (np.array)    :
      J (np.array)    :
      V (np.array)    :
      CperShell (dict. of np.arrays):
      JperShell (dict. of np.arrays):
      VperShell (dict. of np.arrays):
    """
    def __init__(self, element, filename):
        self.filename  = filename
        self.element   = element
        self.elementNr = xrs_utilities.element(z)

        self.shells    = []
        self.edges     = []
        self.C         = []
        self.J         = []
        self.V         = []
        self.CperShell = {}
        self.JperShell = {}
        self.VperShell = {}

class ComptonProfiles:
    """Class for multiple HF Compton profiles.

    This class should hold one or more instances of the ComptonProfile class
    and have methods to return profiles from single atoms, single shells, all
    atoms. It should be able to apply corrections etc. on those... 

    Attributes:
    -----------
      element (string): Element symbol as in the periodic table.
      elementNr (int) : Number of the element as in the periodic table.
      shells (list)   :
      edges (list)    :
      C (np.array)    :
      J (np.array)    :
      V (np.array)    :
      CperShell (dict. of np.arrays):
      JperShell (dict. of np.arrays):
      VperShell (dict. of np.arrays):
    """
    def __init__(self, element):
        self.element   = element
        self.elementNr = xrs_utilities.element(z)
        self.shells    = []
        self.edges     = []
        self.C         = []
        self.J         = []
        self.V         = []
        self.CperShell = {}
        self.JperShell = {}
        self.VperShell = {}


def trapz_weights(x):
    dx = np.diff(x)
    w = np.empty(x.shape)
    w[1:-1] = (dx[1:] + dx[:-1])/2.
    w[0] = dx[0] / 2.
    w[-1] = dx[-1] / 2.
    return w

def PzProfile(element,filename):
    """Returnes tabulated HF Compton profiles.

    Reads in tabulated HF Compton profiles from the Biggs paper,
    interpolates them, and normalizes them to the # of electrons 
    in the shell.

    Args:
    -----
      element (string):  element symbol (e.g. 'Si', 'Al', etc.)
      filename (string): absolute path and filename to tabulated profiles

    Returns:
    --------
      CP_profile (np.array): Matrix of the Compton profile
        1. column: pz-scale
        2. ... n. columns: Compton profile of nth shell
      binding_energy (list): binding energies of shells
      occupation_num (list): number of electrons in the according shells
    """
    # load Biggs data, mirror at pz = 0.0
    CP_tab, occupation_num, binding_energies, shell_names = xrs_fileIO.readbiggsdata(filename,element)
    pz_tab = np.append(-1.0*np.flipud(CP_tab[1::,0]),CP_tab[:,0])
    CP_tab = np.append(np.flipud(CP_tab[1::,:]),CP_tab,axis=0)
    CP_tab[:,0] = pz_tab

    # pad with zeros as large pos. and large neg. pz for nicer spline
    CP_tab = np.append(np.zeros((1,CP_tab.shape[1])),CP_tab,axis=0)
    CP_tab = np.append(CP_tab,np.zeros((1,CP_tab.shape[1])),axis=0)
    CP_tab[0,0]  = -10000.0
    CP_tab[-1,0] = 10000.0

    # interpolate
    pz_scale  = np.arange(-100.0,100.0,0.01)
    CP_profile = np.zeros((len(pz_scale),len(binding_energies)+1))
    CP_profile[:,0] = pz_scale
    for n in range(len(binding_energies)):
        interp_func = interpolate.pchip(CP_tab[:,0], CP_tab[:,n+2])
        CP_profile[:,n+1]  = interp_func(pz_scale)

    # normalize to one electron, multiply by number of electrons
    for n in range(len(binding_energies)):
        norm              = np.trapz(CP_profile[:,n+1],CP_profile[:,0])
        CP_profile[:,n+1] = CP_profile[:,n+1]/norm*long(occupation_num[n])
	binding_energies    = [float(energy) for energy in binding_energies]
	occupation_num    = [float(value) for value in occupation_num]

    return CP_profile, binding_energies, occupation_num, shell_names

def elossProfile(element,filename,E0,tth,correctasym=None,valence_cutoff=20.0):
    """Returns HF Compton profiles on energy loss scale.

    Uses the PzProfile function to read read in Biggs HF profiles
    and converts them onto energy loss scale. The profiles are cut
    at the respective electron binding energies and are normalized
    to the f-sum rule (i.e. S(q,w) is in units of [1/eV]).

    Args:
    -----
    element (string): element symbol.
    filename (string): absolute path and filename to tabulated Compton profiles.
    E0 (float): analyzer energy in [keV].
    tth (float): scattering angle two theta in [deg].
    correctasym (np.array): vector of scaling factors to be applied.
    valence_cutoff (float): energy value below which edges are considered as valence

    Returns:
    --------
    enScale (np.array): energy loss scale in [eV]
    J_total (np.array): total S(q,w) in [1/eV]
    C_total (np.array): core contribution to S(q,w) in [1/eV]
    V_total (np.array): valence contribution to S(q,w) in [1/eV], the valence is defined by valence_cutoff
    q (np.array): momentum transfer in [a.u]
    J_shell (dict of np.arrays): dictionary of contributions for each shell, the key are defines as in Biggs table.
    C_shell (dict of np.arrays): same as J_shell for core contribution
    V_shell (dict of np.arrays): same as J_shell for valence contribution

    """
    # read in the Biggs data
    CP_profile, binding_energies, occupation_num, shell_names = PzProfile(element,filename)

    # convert pz to energy loss scale
    enScale = ((np.flipud(xrs_utilities.pz2e1(E0,CP_profile[:,0],tth))-E0)*1.0e3)

    # define the momentum transfer
    q = xrs_utilities.momtrans_au(enScale/1000.0 + E0, E0, tth)

    # calculate asymmetry after Holm and Ribberfors for filles 1s and 2p shells
    # if correctasym == True
    asymmetry = np.flipud(xrs_utilities.HRcorrect(CP_profile, occupation_num, q))
    if correctasym:
        CP_profile[:,1:4] = CP_profile[:1:4] + asymmetry * correctasym

    # discard profiles, q, and enScale for energy losses smaller than zero
    HF_profile = CP_profile[np.nonzero(enScale.T>=0.0)[0],:]
    q          = q[np.nonzero(enScale.T>=0)[0]]
    enScale    = enScale[np.nonzero(enScale.T>=0)[0]]
    HF_profile[:,0] = enScale 

    # discard profiles for energy losses  below according binding energies
    for n in range(len(binding_energies)):
        HF_profile[np.where(enScale<binding_energies[n]),n+1] = 0 

    # convert J(pz) to S(q,w) via J(pz)=N_electrons*hartree*q*S(q,w) and
    # normalize using the f-sum rule (sum(S(q,w)*w)=f)
    # 1. convert to a.u.
    hartree  = 1.0/constants.physical_constants['electron volt-hartree relationship'][0]
    enScaleH = enScale/hartree # eloss in a.u.

    # 2. normalize to one then multiply by N_el*q**2.0/2.0
    for n in range(len(binding_energies)):
        HF_profile[:,n+1] = HF_profile[:,n+1]/(integrate.trapz(np.multiply(HF_profile[:,n+1],enScaleH),enScaleH))
        HF_profile[:,n+1] = np.multiply(HF_profile[:,n+1],(q**2.0)/2.0)*occupation_num[n]

    # 3. convert back to [1/eV] and sum up
    J_total = np.zeros((len(enScale)))
    V_total = np.zeros((len(enScale)))
    for n in range(len(binding_energies)):
        if binding_energies[n] < enScale[-1]:
            J_total += HF_profile[:,n+1]/hartree
            if binding_energies[n] < valence_cutoff:
                V_total += HF_profile[:,n+1]/hartree
    C_total = J_total - V_total

    # make dictionaries for individual shells
    J_shell = {}
    C_shell = {}
    V_shell = {}
    counter = 1
    for name in shell_names:
        if 'Shell' in name:
            if binding_energies[counter-1] < enScale[-1]:
                J_shell[name] = HF_profile[:,counter]
                if binding_energies[counter-1] < valence_cutoff:
                    V_shell[name] = HF_profile[:,counter]
                else:
                    V_shell[name] = np.zeros_like(HF_profile[:,counter])
            C_shell[name] = J_shell[name] - V_shell[name]
            counter += 1

    return enScale, J_total, C_total, V_total, q, J_shell, C_shell, V_shell
