#!/usr/bin/python
# Filename: xrs_extraction.py

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

import os
import numpy as np
import xrs_utilities, xrs_fileIO, math_functions, xrs_ComptonProfiles

installation_dir = os.path.dirname(os.path.abspath(__file__))

debug = 1

if not debug:
	HFCP_PATH = os.path.join(installation_dir,'../../../../share/XRStools/data/ComptonProfiles.dat')
	LOGTABLE_PATH = os.path.join(installation_dir,'../../../../share/XRStools/data/logtable.dat')
else:
	HFCP_PATH     = '/home/christoph/sources/XRStools/data/ComptonProfiles.dat'
	LOGTABLE_PATH = '/home/christoph/sources/XRStools/data/logtable.dat'


class HF_dataset:
	"""
	**dataset**
	A class to hold all information from HF Compton profiles necessary to subtract background from the experiment.
	"""
	def __init__(self, data, formulas, stoich_weights, edges):
		self.formulas       = formulas
		self.stoich_weights = stoich_weights
		self.edges          = edges #e.g. {'Li':['K','L23'], 'O':'K'}
		self.E0             = data.E0
		self.cenom          = data.cenom
		self.tth            = data.tth
		self.eloss          = data.eloss
		self.HFProfile      = xrs_ComptonProfiles.HFProfile(formulas, stoich_weights, HFCP_PATH)
		self.HFProfile.get_elossProfiles(self.E0,self.tth)

		# interpolate total HF profiles onto experimental eloss scale
		self.J_total   = np.zeros((len(self.eloss),len(self.tth)))
		self.C_total   = np.zeros((len(self.eloss),len(self.tth)))
		self.V_total   = np.zeros((len(self.eloss),len(self.tth)))
		for ii in range(len(self.tth)):
			self.J_total[:,ii] = np.interp(self.eloss, self.HFProfile.eloss,self.HFProfile.J_total[:,ii])
			self.C_total[:,ii] = np.interp(self.eloss, self.HFProfile.eloss,self.HFProfile.C_total[:,ii])
			self.V_total[:,ii] = np.interp(self.eloss, self.HFProfile.eloss,self.HFProfile.V_total[:,ii])

		# initialize double dicts for {'element1':{'edge1','edge2',...}, 'element2':{'edge1','edge2',...} }
		self.C_edges = {}
		for key in self.edges:
			self.C_edges[key] = {}
			for edge in self.edges[key]:
				self.C_edges[key][edge] = np.zeros((len(self.eloss),len(self.tth)))

		# interpolate the core profiles for the desired elements
		for key in self.edges:
			for edge in self.edges[key]:
				edge_keyword = xrs_ComptonProfiles.mapShellNames(edge,xrs_utilities.element(key))
				for formula in self.formulas:
					if key in formula:
						# cp core-edge profile
						atom_profile = self.HFProfile.FormulaProfiles[formula].AtomProfiles[key]
						for jj in range(len(self.tth)):
							self.C_edges[key][edge][:,jj] = np.interp(self.eloss,atom_profile.eloss,atom_profile.CperShell[edge_keyword][:,jj])
					else:
						print('Could not find ' + key + ' in any of the provided formulas.')

class valence_CP:
	"""
	**valence_CP**
	Class to organize information about extracted experimental valence Compton profiles.
	"""
	def __init__(self)
		self.pzscale        = np.flipud(np.arange(-10,10,0.05)) # definition according to Huotari et al, JPCS 62 (2001) 2205
		self.valencepz      = np.zeros((len(self.pzscale),len(self.signals[0,:])))
		self.valasymmetrypz = np.zeros((len(self.pzscale),len(self.signals[0,:])))
		self.valence        = np.zeros((len(self.eloss),len(self.signals[0,:])))
		self.valasymmetry   = np.zeros((len(self.eloss),len(self.signals[0,:])))

class edge_extraction:
	"""
	**edge_extraction**
	Class to destill core edge spectra from x-ray Raman scattering experiments.
	"""
	def __init__(self,experimental_data, formulas, stoich_weights, edges ,prenormrange=[5,np.inf]):
		# input
		self.eloss   = experimental_data.eloss
		self.signals = experimental_data.signals
		self.errors  = experimental_data.errors
		self.E0      = experimental_data.E0
		self.tth     = experimental_data.tth
		self.prenormrange = prenormrange
		self.HF_dataset = HF_dataset(experimental_data, formulas, stoich_weights, edges)

		# output
		self.background     = np.zeros(np.shape(data.signals))
		self.sqw            = np.zeros(np.shape(data.signals))
		self.sqwav          = np.zeros(np.shape(data.eloss))
		self.sqwaverr       = np.zeros(np.shape(data.eloss))
		self.valence_CP     = valence_CP()

		# some variables for averaging rawdata over analyzers/whole chambers
		self.avsignals  = np.array([])
		self.averrors   = np.array([])
		self.avC        = np.array([])
		self.avqvals    = np.array([])

		# rough normalization over range given by prenormrange
		if prenormrange:
			for n in range(len(self.tth)):
				HFnorm  = np.trapz(self.HF_dataset.J_total[:,n], self.eloss)
				inds    = np.where(np.logical_and(self.eloss >= prenormrange[0], self.eloss <= prenormrange[1]))[0]
				EXPnorm = np.trapz(self.signals[inds,n], self.eloss[inds])
				self.signals[:,n] = self.signals[:,n]/EXPnorm*HFnorm
				self.errors[:,n]  = self.errors[:,n]/EXPnorm*HFnorm

	def analyzerAverage(self,roi_numbers,errorweighing=True):
		"""
		**analyzerAverage**
		Averages signals from several crystals before background subtraction.

		Args:
        -----
		roi_numbers : list, str
			list of ROI numbers to average over of keyword for analyzer chamber
			(e.g. 'VD','VU','VB','HR','HL','HB')
		errorweighing : boolean (True by default)
			keyword if error weighing should be used for the averaging or not

		"""	
		if isinstance(roi_numbers,list):
			columns = roi_numbers
		elif: isinstance(roi_numbers,str): 
			columns = map_chamber_names(roi_numbers)
		else:
			print('Unsupported type for keyword \'roi_numbers\'.')
			return

		self.avsignals = np.zeros_like(self.eloss)
		self.averror   = np.zeros_like(self.eloss)
		self.avC       = np.zeros_like(self.eloss)
		self.avqvals   = np.zeros_like(self.eloss)

		# build matricies to sum over
		av      = np.zeros((len(self.eloss),len(columns)))
		averr   = np.zeros((len(self.eloss),len(columns)))
		avcvals = np.zeros((len(self.eloss),len(columns)))
		avqvals = np.zeros((len(self.eloss),len(columns)))

		for column,ii in zip(columns,range(len(columns))):
			# find data points with error = 0.0 and replace error by 1.0
			inds = np.where(self.errors[:,column] == 0.0)[0]
			self.errors[inds,column] = 1.0
			av[:,ii]      = self.signals[:,column]]
			averr[:,ii]   = self.errors[:,column]
			avqvals[:,ii] = self.qvals[:,column]

		# sum things up
		if errorweighing:
			self.avsignals = np.sum( av/averr**2.0 ,axis=1)/( np.sum(1.0/averr**2.0,axis=1))
			self.averrors  = np.sqrt( 1.0/np.sum(1.0/(averr)**2.0,axis=1) )
		else: 
			self.avsignals = np.sum(av,axis=1)
			self.averrors  = np.sqrt(np.sum(np.absolute(averr)**2.0,axis=1)) # check this again

		self.avqvals = np.mean(avqvals,axis=1)

	def removeCorePearsonAv2(self,range1,range2,element,edge,weights=[2,1],guess=None):
		"""
		**removeCorePearsonAv2**
		
		"""
		# find indices for fit-range
		region1 = np.where(np.logical_and(self.eloss >= range1[0], self.eloss <= range1[1]))[0]
		region2 = np.where(np.logical_and(self.eloss >= range2[0], self.eloss <= range2[1]))[0]
		region  = np.append(region1*weights[0],region2*weights[1])

		# find indices for guessing start values
		guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]
		if not guess: 
			guess = []
			ind   = self.avsignals[guessregion].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
			guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
			guess.append(guess[0]*1.0) # once the position of the peason maximum
			guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
			guess.append(self.avsignals[guessregion][ind]) # Peak intensity
			guess.append(0.0) # linear slope
			guess.append(0.0) # no background
			guess.append(1.0) # scaling factor for exp. data

		# manage some plotting things
		plt.ion()
		plt.cla()

		# get the HF core spectrum
		HF_core = self.HF_dataset.C_edges[element][edge]

		bnds = ((0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None))
		cons    =  ({'type': 'ineq', 'fun': lambda x:  x[2]  },
					{'type': 'ineq', 'fun': lambda x:  x[3]  },
					{'type': 'ineq', 'fun': lambda x:  x[6]  },
					{'type': 'eq',   'fun': lambda x:   np.sum( (x[6]*self.avsignals[region2] - 
																pearson7_zeroback(self.eloss[region2],x[0:4]) - 
																np.polyval(x[4:6],self.eloss[region2]) - 
																self.avC[region2])**2.0 )  },

)
res_ca350_50h = optimize.minimize(fitfctn, (0.1,0.1), method='SLSQP', bounds=bnds,constraints=cons).x

fitfct  = lambda a: np.sum( (a[6]*self.avsignals[region1] - pearson7_zeroback(self.eloss[region1],a[0:4]) - np.polyval(a[4:6],self.eloss[region1]) - self.avC[region1])**2.0 ) + np.sum( (a[6]*self.avsignals[region2] - pearson7_zeroback(self.eloss[region2],a[0:4]) - np.polyval(a[4:6],self.eloss[region2]) - self.avC[region2])**2.0 )


		# some sensible boundary conditions for the fit:

		c2 = lambda a: a[2] # shape should not be negative
		c3 = lambda a: a[3] # peak intensity should not be negative
		c4 = lambda a: np.absolute(5e-1 - a[4]) # slope for linear background should be small
		c5 = lambda a: a[3] - a[5] # offset for linear should be smaller than maximum of pearson
		c6 = lambda a: a[6] # scaling factor for the data should not be negative
		c7 = lambda a: np.sum( (a[6]*self.avsignals[region2] - pearson7_zeroback(self.eloss[region2],a[0:4]) - self.avC[region2] - np.polyval(a[4:6],self.eloss[region2]))**2.0 )

		fitfct  = lambda a: np.sum( (a[6]*self.avsignals[region1] - pearson7_zeroback(self.eloss[region1],a[0:4]) - np.polyval(a[4:6],self.eloss[region1]) - self.avC[region1])**2.0 ) + np.sum( (a[6]*self.avsignals[region2] - pearson7_zeroback(self.eloss[region2],a[0:4]) - np.polyval(a[4:6],self.eloss[region2]) - self.avC[region2])**2.0 )
		cons = []#[c1, c2, c3, c4, c5, c6, c7]
		res     = optimize.minimize(fitfct,guess)
		yres    = pearson7_zeroback(self.eloss,res[0:4]) + np.polyval(res[4:6],self.eloss)

		plt.plot(self.eloss,self.avsignals*res[6],self.eloss,yres+self.avC,self.eloss,self.avsignals*res[6]-yres,self.eloss,self.avC)
		plt.legend(('scaled data','pearson + linear + core','data - (pearson + linear)','core'))
		plt.draw()

		self.sqwav = self.avsignals*res[6] - yres
		self.sqwaverr = self.averrors*res[6]




def map_chamber_names(name):
	"""
	**map_chamber_names**
	Maps names of chambers to range of ROI numbers.
	"""
	chamberNames = {'VD':range(0,12),
					'VU':range(12,24),
					'VB':range(24,36),
					'HR':range(36,48),
					'HL':range(48,60),
					'HB':range(60,72)}

	return chamberNames[name.upper()]





