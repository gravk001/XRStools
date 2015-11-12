#!/usr/bin/python
# Filename: xrs_scans.py

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
from itertools import groupby
import xrs_utilities, math_functions

class scan:
	"""
	Container class, holding information of single scans performed with 2D detectors. 
	"""
	def __init__(self,edf_arrays,scannumber,energy_scale,monitor_signal,counters,motor_positions,specfile_data,scantype='generic'):
		# rawdata
		self.edfmats  = np.array(edf_arrays)          # numpy array of all 2D images that belong to the scan
		self.number   = scannumber                    # number under which this scan can be found in the SPEC file 
		self.scantype = scantype                      # keyword, later used to group scans (add similar scans, etc.)
		self.energy   = np.array(energy_scale)        # numpy array of the energy axis
		self.monitor  = np.array(monitor_signal)      # numpy array of the monitor signal
		# some things maybe not imediately necessary 
		self.counters = counters                      # names of all counters that appear in the SPEC file for this scan (maybe unnecessary)
		self.motors   = motor_positions               # all motor positions as found in the SPEC file header for this scan ( " )
		self.specdata = np.array(specfile_data)       # all data that is also in the SPEC file for this scan
		# data (to be filled after defining rois)
		self.eloss    = []     # numpy array of the energy loss scale
		self.signals  = []     # numpy array of signals extracted from the ROIs
		self.errors   = []     # numpy array of all Poisson errors
		self.cenom    = []     # list with center of masses (used if scan is an elastic line scan)
		# pixel wise things
		self.signals_pw = []   
		self.cenom_pw   = []
		self.signals_pw_interp = []

	def applyrois(self,indices,scaling=None):
		"""
		Sums up intensities found in the ROIs of each detector image
		and stores it into the self.signals attribute.
		roi_object = instance of the 'rois' class redining the ROIs
		scaling    = numpy array of numerical scaling factors (has to be one for each ROIs)
		"""
		data = np.zeros((len(self.edfmats),len(indices)))
		for n in range(len(indices)): # each roi
			for m in range(len(self.edfmats)): # each energy point along the scan
				for l in range(len(indices[n])): # each pixel on the detector
					data[m,n] += self.edfmats[m,indices[n][l][0],indices[n][l][1]]
		self.signals = np.array(data)
		if np.any(scaling):
			assert len(scaling) == len(indices) # make sure, there is one scaling factor for each roi
			for ii in range(len(indices)):
				self.signals[:,ii] *= scaling[ii]

	def applyrois_pw(self,indices,scaling=None):
		"""
		Pixel-wise reading of the ROI's pixels into a list of arrays. I.e. each n-pixel ROI will have n Spectra, saved in a 2D array.
		Parameters
		----------
		indices : list
			List of indices (attribute of the xrs_rois class).
		scaling : list of flaots, optional
			Python list of scaling factors (one per ROI defined) to be applied to all pixels of that ROI.

		"""
		data = [] # list of 2D arrays (energy vs. intensity for each pixel inside a single ROI) 
		for n in range(len(indices)): # each ROI
			roidata = np.zeros((len(self.edfmats),len(indices[n]))) # 2D np array energy vs pixels in current roi
			for m in range(len(self.edfmats)): # each energy point along the scan
				for l in range(len(indices[n])): # each pixel on the detector
					roidata[m,l] = self.edfmats[m,indices[n][l][0],indices[n][l][1]]
			data.append(roidata) # list which contains one array (energy point, pixel) per ROI
		self.signals_pw = data

	def get_eloss_pw(self):
		"""
		Finds the center of mass for each pixel in each ROI, sets the energy loss scale in 
		and interpolates the signals to a common energy loss scale. Finds the resolution (FWHM) for
		each pixel.

		"""
		if not self.scantype == 'elastic': # return, if scantype is not elastic
			print 'This method is meant for elastic line scans only!'
			return
		if not self.signals_pw: # return, if there is no data
			print 'Please use the applyrois_pw function first!'
			return
		else:
			cenom_pw      = []
			resolution_pw = []
			for roiind in range(len(self.signals_pw)): # each ROI
				oneroi_cenom      = []
				oneroi_resolution = []
				for pixelind in range(len(self.signals_pw[roiind])): # each pixel in the ROI
					oneroi_cenom.append(xrs_utilities.find_center_of_mass(self.energy,self.signals_pw[roiind][:,pixelind]))
					try: 
						FWHM,x0 = xrs_utilities.fwhm((self.energy - oneroi_cenom[roiind][pixelind])*1e3,self.signals_pw[roiind][:,pixelind])
						oneroi_resolution.append(FWHM)
					except:
						oneroi_resolution.append(0.0)
				cenom_pw.append(oneroi_cenom)
				resolution_pw.append(oneroi_resolution)
			# define mean of first ROI as 'master' energy loss scale for all Pixels
			self.eloss = (self.energy - np.mean(cenom_pw[0]))*1e3 # average energy loss in eV from first ROI
			# !!!! mayb it is better to just save the CENOM values and to the interpolation later with the stiched 
			# !!!! whole spectrum... 
			# define eloss-scale for each ROI and interpolate onto the 'master' eloss-scale
			#for roiind in range(len(self.signals_pw)): # each ROI
			#	oneroi_signals = np.zeros_like(self.signals_pw[roiind])
			#	for pixelind in range(len(self.signals_pw[roiind])):
			#		x = (self.energy-self.cenom_pw[roiind][pixelind])*1e3
			#		y =
			#		f = interp1d(x,y, bounds_error=False,fill_value=0.0)
			#		self.signals[:,n] = f(self.eloss)
			#self.signals_pw_interp = []

	def get_type(self):
		return self.scantype

	def get_scannumber(self):
		return self.number

	def get_shape(self):
		if not np.any(self.signals):
			print 'please apply the ROIs first.'
			return
		else:
			return np.shape(self.signals)

	def get_numofrois(self):
		if not self.signals.any():
			print 'please apply the ROIs first.'
			return
		else:
			return np.shape(self.signals)[1]

class scangroup:
	"""
	Container class holding information from a group of scans.
	"""
	def __init__(self,energy,signals,errors,grouptype='generic'):
		self.energy    = energy     # common energy scale
		self.eloss     = []         # common energy loss scale
		self.signals   = signals    # summed up signals
		self.errors    = errors     # Poisson errors
		self.grouptype = grouptype  # keyword common to all scans
		self.signals_orig = signals # keep a copy of uninterpolated data

	def get_type(self):
		return self.grouptype

	def get_meanenergy(self):
		return np.mean(self.energy)

	def get_estart(self):
		return self.energy[0]

	def get_eend(self):
		return self.energy[-1]

	def get_meanegridspacing(self):
		return np.mean(np.diff(self.energy))

	def get_maxediff(self):
		return (self.energy[-1]-self.energy[0])

class offDiaDataSet:
	""" **offDiaDataSet**
	Class to hold information from an off-diagonal dataset.
	"""
	def __init__(self):
		self.RCmonitor    = np.array([])
		self.signalMatrix = np.array([])
		self.motorMatrix  = np.array([])
		self.energy       = np.array([])
		self.ROI_number   = 0
		self.G_vector     = np.array([])
		self.q0           = np.array([])
		self.qh           = np.array([])
		self.k0           = np.array([])
		self.kh           = np.array([])
		self.kprime       = np.array([])

def findgroups(scans):
	"""
	this groups together instances of the scan class based on their  "scantype" attribute and returns ordered scans
	"""	
	allscannames = []
	for scan in scans:
		print scan
		allscannames.append(scan)
	allscannames.sort() # 
	allscans = []
	for scan in allscannames:
		 allscans.append(scans[scan])
	allscans = sorted(allscans,key=lambda x:x.get_type())
	rawgroups = []
	results = groupby(allscans,lambda x:x.get_type())
	print 'the following scangroups were found:'
	for key,thegroup in results:
		print key
		rawgroups.append(list(thegroup))
	return rawgroups

def makegroup(groupofscans,grouptype=None):
	"""
	takes a group of scans, sums up the signals and monitors, estimates poisson errors, and returns an instance of the scangroup 		class (turns several instances of the "scan" class into an instance of the "scangroup" class)
	"""
	if not grouptype:
		grouptype = groupofscans[0].get_type() # the type of the sum of scans will be the same as the first from the list	
	theenergy   = groupofscans[0].energy # all scans are splined onto energy grid of the first scan
	thesignals  = np.zeros(groupofscans[0].get_shape())
	theerrors   = np.zeros(groupofscans[0].get_shape())
	themonitors = np.zeros(np.shape(groupofscans[0].monitor))
	for scan in groupofscans:
		f  = interpolate.interp1d(scan.energy,scan.monitor, bounds_error=False, fill_value=0.0)
		moni = f(theenergy)
		themonitors += moni
		for n in range(thesignals.shape[-1]):
			f = interpolate.interp1d(scan.energy,scan.signals[:,n], bounds_error=False, fill_value=0.0)
			signal = f(theenergy)
			thesignals[:,n] += signal
	for n in range(thesignals.shape[-1]):
		theerrors[:,n]  = np.sqrt(thesignals[:,n])
	# and normalize	
	for n in range(thesignals.shape[-1]):
		thesignals[:,n] = thesignals[:,n]/themonitors
		theerrors[:,n]  = theerrors[:,n]/themonitors

	group = scangroup(theenergy,thesignals,theerrors,grouptype)
	return group

def makegroup_nointerp(groupofscans,grouptype=None):
	"""
	takes a group of scans, sums up the signals and monitors, estimates poisson errors, and returns an instance of the scangroup class (turns several instances of the "scan" class into an instance of the "scangroup" class), same as makegroup but withouth interpolation to account for encoder differences... may need to add some "linspace" function in case the energy scale is not monotoneous... 
	"""
	if not grouptype:
		grouptype = groupofscans[0].get_type() # the type of the sum of scans will be the same as the first from the list	
	theenergy   = groupofscans[0].energy
	thesignals  = np.zeros(groupofscans[0].get_shape())
	theerrors   = np.zeros(groupofscans[0].get_shape())
	themonitors = np.zeros(np.shape(groupofscans[0].monitor))

	for scan in groupofscans:
		themonitors += scan.monitor
		for n in range(thesignals.shape[-1]):
			thesignals[:,n] += scan.signals[:,n]
	for n in range(thesignals.shape[-1]):
		theerrors[:,n]  = np.sqrt(thesignals[:,n])
	# and normalize
	for n in range(thesignals.shape[-1]):
		thesignals[:,n] = thesignals[:,n]/themonitors
		theerrors[:,n]  = theerrors[:,n]/themonitors
	group = scangroup(theenergy,thesignals,theerrors,grouptype)
	return group

def append2Scan_right(group1,group2,inds=None,grouptype='spectrum'):
	"""
	append two instancees of the scangroup class, return instance of scangroup
	append group2[inds] to the right (higher energies) of group1
	if inds is not None, only append rows indicated by inds to the first group 
	"""
	assert isinstance(group1,scangroup) and isinstance(group2,scangroup)
	energy  = group1.energy
	signals = group1.signals
	errors  = group1.errors
	gtype   = group1.grouptype
	if not inds:
		energy  = np.append(energy,group2.energy)
		signals = np.append(signals,np.squeeze(group2.signals),0)
		errors  = np.append(errors,np.squeeze(group2.errors),0)
		return scangroup(energy,signals,errors,grouptype=gtype)
	else:
		energy  = np.append(energy,group2.energy[inds])
		signals = np.append(signals,np.squeeze(group2.signals[inds,:]),0)
		errors  = np.append(errors,np.squeeze(group2.errors[inds,:]),0)
		return scangroup(energy,signals,errors,grouptype=gtype)

def append2Scan_left(group1,group2,inds=None,grouptype='spectrum'):
	"""
	append two instancees of the scangroup class, return instance of scangroup
	append group1[inds] to the left (lower energies) of group2
	if inds is not None, only append rows indicated by inds to the first group 
	"""
	assert isinstance(group1,scangroup) and isinstance(group2,scangroup)
	if not inds:
		energy  = group1.energy
		signals = group1.signals
		errors  = group1.errors
		gtype   = group1.grouptype
		energy  = np.append(energy,group2.energy)
		signals = np.append(signals,np.squeeze(group2.signals),0)
		errors  = np.append(errors,np.squeeze(group2.errors),0)
		return scangroup(energy,signals,errors,grouptype=gtype)
	else:
		energy  = group1.energy[inds]
		signals = group1.signals[inds,:]
		errors  = group1.errors[inds,:]
		gtype   = group1.grouptype
		energy  = np.append(energy,group2.energy)
		signals = np.append(signals,np.squeeze(group2.signals),0)
		errors  = np.append(errors,np.squeeze(group2.errors),0)
		return scangroup(energy,signals,errors,grouptype=gtype)

def insertScan(group1,group2,grouptype='spectrum'):
	"""
	inserts group2 into group1
	NOTE! there is a numpy insert function, maybe it would be better to use that one!
	"""
	# find indices for below and above group2
	lowinds  = np.where(group1.energy<group2.get_estart())
	highinds = np.where(group1.energy>group2.get_eend())

	energy  = np.append(group1.energy[lowinds],np.append(group2.energy,group1.energy[highinds]))
	signals = np.append(np.squeeze(group1.signals[lowinds,:]),np.append(group2.signals,np.squeeze(group1.signals[highinds,:]),0),0)
	errors  = np.append(np.squeeze(group1.errors[lowinds,:]),np.append(group2.errors,np.squeeze(group1.errors[highinds,:]),0),0)

	return scangroup(energy,signals,errors,grouptype)

def catScansLong(groups,include_elastic):
	"""
	takes a longscan and inserts other backgroundscans (scans that have 'long' in their name) and other scans and inserts them into the long scan. 
	"""
	# the long scan
	spectrum = groups['long']

	# groups that don't have 'long' in the grouptype	
	allgroups  = []
	for group in groups:
		if not 'long' in group:
			allgroups.append(groups[group])
	allgroups.sort(key = lambda x:x.get_estart())

	# groups that have 'long' in the grouptype	
	longgroups = []
	for group in groups:
		if 'long' in group and group != 'long':
			longgroups.append(groups[group])
	longgroups.sort(key = lambda x:x.get_estart())

	# if there are other longscans: insert those first into the long scan
	for group in longgroups:
		spectrum = insertScan(spectrum,group)

	# insert other scans into the long scan
	for group in allgroups:
		spectrum = insertScan(spectrum,group)

	# cut off the elastic line if present in the groups
	if 'elastic' in groups and not include_elastic:
		inds = np.where(spectrum.energy > groups['elastic'].get_eend())[0]
		return spectrum.energy[inds], spectrum.signals[inds,:], spectrum.errors[inds,:]
	else:
		return spectrum.energy, spectrum.signals, spectrum.errors

def catScans(groups,include_elastic):
	"""
	concatenate all scans in groups, return the appended energy, signals, and errors
	"""
	# sort the groups by their start-energy
	allgroups  = []
	for group in groups:
		allgroups.append(groups[group])
	allgroups.sort(key = lambda x:x.get_estart())

	# assign first group to the spectrum (which is an instance of the scangroup class as well)
	spectrum = scangroup(allgroups[0].energy,allgroups[0].signals,allgroups[0].errors,grouptype='spectrum')

	# go through all other groups and append them to the right of the spectrum
	if len(allgroups)>1: # check if there are groups to append
		for group in allgroups[1:]:
			spectrum = append2Scan_right(spectrum,group)

	# cut off the elastic line if present in the groups
	if 'elastic' in groups and not include_elastic:
		inds = np.where(spectrum.energy > groups['elastic'].get_eend())[0]
		return spectrum.energy[inds], spectrum.signals[inds,:], spectrum.errors[inds,:]
	else:
		return spectrum.energy, spectrum.signals, spectrum.errors

def appendScans(groups,include_elastic):
    """
    try including different background scans... 
    append groups of scans ordered by their first energy value.
    long scans are inserted into gaps that at greater than two times the grid of the finer scans
    """
    # find all types of groups	
    grouptypes = [key for key in groups.keys()]

    if 'long' in grouptypes:
        print " going to refine "
        return catScansLong(groups,include_elastic)
    else: 
        return catScans(groups,include_elastic)

def catXESScans(groups):
	"""
	Concatenate all scans in groups, return the appended energy, signals, and errors.
	This needs to be a bit smarter to also work for scans that are scanned from small to large energy...
	"""
	# sort the groups by their end-energy (is the smallest in XES/energy2-scans)
	allgroups  = []
	for group in groups:
		allgroups.append(groups[group])
	allgroups.sort(key = lambda x:x.get_eend())

	# find lowest energy, highest energy, smallest increment, define grid
	eStart = allgroups[-1].energy[0]
	eStop  = allgroups[0].energy[-1]
	stepSizes = []
	for group in allgroups:
		stepSizes.append(np.diff(group.energy)[0])
	eStep   = np.amin(stepSizes)

	energy  = np.arange(eStop,eStart,-eStep)
	signals = np.zeros((len(energy),1))
	errors  = np.zeros((len(energy),1))

	# interpolate all groups onto grid, put into a matrix
	for group in allgroups:
		interp_signals = np.interp(energy,np.flipud(group.energy),np.flipud(np.squeeze(group.signals)), left=0.0, right=0.0)
		signals        = np.append(signals, interp_signals.reshape((len(energy),1)),axis=1)
		interp_errors  = np.interp(energy,np.flipud(group.energy),np.flipud(np.squeeze(group.errors)), left=0.0, right=0.0)
		errors         = np.append(errors,interp_errors.reshape((len(energy),1)),axis=1)

	# sum up and weigh by number of available non-zero points
	sum_signals = np.zeros_like(energy)
	sum_errors  = np.zeros_like(energy)
	for ii in range(len(energy)):
		nParts = len(np.where(signals[ii,:]>0.0)[0])
		sum_signals[ii] = np.sum(signals[ii,:])/nParts
		sum_errors[ii]  = np.sqrt(np.sum(errors[ii,:]**2.0))/nParts

	spectrum = scangroup(energy,sum_signals,sum_errors,grouptype='spectrum')
	return spectrum.energy, spectrum.signals, spectrum.errors

def appendXESScans(groups):
	"""
	try including different background scans... 
	append groups of scans ordered by their first energy value.
	long scans are inserted into gaps that at greater than two times the grid of the finer scans
	"""
	# find all types of groups	
	grouptypes = [key for key in groups.keys()]
	return catXESScans(groups)

def create_sum_image(scans,scannumbers):
	"""
	Returns a summed image from all scans with numbers 'scannumbers'.
	scans       = dictionary of objects from the scan-class
	scannumbers = single scannumber, or list of scannumbers from which an image should be constructed
	"""
	# make 'scannumbers' iterable (even if it is just an integer)
	numbers = []
	if not isinstance(scannumbers,list):
		numbers.append(scannumbers)
	else:
		numbers = scannumbers

	key = 'Scan%03d' % numbers[0]
	image = np.zeros_like(scans[key].edfmats[0,:,:])
	for number in numbers:
		key = 'Scan%03d' % number
		for ii in range(scans[key].edfmats.shape[0]):
			image += scans[key].edfmats[ii,:,:]

	return image


def create_diff_image(scans,scannumbers,energy_keV):
	"""
	Returns a summed image from all scans with numbers 'scannumbers'.
	scans       = dictionary of objects from the scan-class
	scannumbers = single scannumber, or list of scannumbers from which an image should be constructed
	"""
	# make 'scannumbers' iterable (even if it is just an integer)
	numbers = []
	if not isinstance(scannumbers,list):
		numbers.append(scannumbers)
	else:
		numbers = scannumbers

	key = 'Scan%03d' % numbers[0]
	below_image = np.zeros_like(scans[key].edfmats[0,:,:])
	above_image = np.zeros_like(scans[key].edfmats[0,:,:])

	# find indices below and above 'energy'
	below_inds = scans[key].energy < energy_keV
	above_inds = scans[key].energy > energy_keV
	for number in numbers:
		key = 'Scan%03d' % number
		for ii in below_inds:
			below_image += scans[key].edfmats[ii,:,:]
		for ii in above_inds:
			above_image += scans[key].edfmats[ii,:,:]

	return (above_image - below_image)

def findRCscans(scans):
	""" **findRCscans**
	Returns a list of scans with name RC.
	"""
	RCscans = []
	for key in scans:
		if scans[key].get_type() == 'RC':
			RCscans.append(scans[key])
	return RCscans



