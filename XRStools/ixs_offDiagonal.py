#!/usr/bin/python
# Filename: ixs_offDiagonal.py

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


#from helpers import *
import xrs_rois, xrs_scans, xrs_utilities, math_functions, xrs_fileIO, roifinder_and_gui

from numpy import array
import scipy.io
import traceback
import sys
import os
import numpy as np
import array as arr
import pickle
from itertools import groupby
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from pylab import *
from scipy import signal
from scipy.ndimage import measurements
from itertools import groupby

import matplotlib.pyplot as plt
import warnings

# try to import the fast PyMCA parsers
try:
    import PyMca5.PyMcaIO.EdfFile as EdfIO
    import PyMca5.PyMcaIO.specfilewrapper as SpecIO
    use_PyMca = True
except:
    use_PyMca = False

print " >>>>>>>>  use_PyMca " , use_PyMca
__metaclass__ = type # new style classes

class offDiagonal:
	""" **offDiagonal**
	Class for reading scans from off-diagonal IXS experiments on the high-resolution setup at ID20.

	Arguments:
	----------
	absFilename (string): Absolute path and filename of the SPEC-file.
	scanMotor (string): Mnemonic of the motor that is scanned.
	monitorName (string): Mnemonic of the counter used for normalization.
	edfName (string): EDF-file base file name (default is None, i.e. same as SPEC-file).
	armLength (float): Legth (in m) of the spectrometer arm used (either 1.0 or 2.0).
	"""
	def __init__(self,absFilename,scanMotor='srz',monitorName='kaprixs',edfName=None,armLength=1.0):

		# SPEC-file name issues
		if absFilename:
			if not os.path.isfile(absFilename):
				raise Exception('IOError! No such file, please check filename.')
			self.path          = os.path.split(absFilename)[0] + '/'
			self.filename = os.path.split(absFilename)[1]
			if not edfName:
				self.edfName = os.path.split(absFilename)[1]
			else:
				self.edfName = edfName
		self.EDF_PREFIX    = 'edf/'
		self.EDF_POSTFIX   = '.edf'

		# scan numbers and dictionaries
		self.scans         = {}
		self.scanNumbers   = []

		# detector info
		self.DET_PIXEL_NUMx = 0
		self.DET_PIXEL_NUMy = 0
		self.DET_PIXEL_NUM  = 0

		# analyzer info
		self.armLength      = armLength # 1.0 m or 2.0 m arm
		self.crystBendR     = armLength # analyzer crystal bending radius
		if self.armLength == 1.0:
			self.TTH_OFFSETS1   = [] # [off-set direction 1, off-set direction 2]
		elif self.armLength == 2.0:
			self.TTH_OFFSETS1   = [] # [off-set direction 1, off-set direction 2]
		else:
			warnings.warn('Only 1.0 m and 2.0 m spectrometer arms available!')

		# ROI object from the roifinder class
		self.roi_obj = []

		# which column in the SPEC file to be used for the energy and monitor
		self.scanMotor     = scanMotor.lower()
		self.moniColumn    = monitorName.lower()

	def readscan(self,scanNumber):
		""" **readscan**
		Reads the SPEC- and EDF-files for a single scan.

		Arguments:
		----------
		scanNumber (int): number of scan to read.

		Returns:
		--------
		data (np.ndarray): Matrix holding data from SPEC-file.
		motors (list): List of motors as in header of SPEC-file.
		counters (dictionary): Dictionary holding same data as in data but as dictionary.
		edfmats (np.ndarray): 3D Matrix with the EDF-files.
		"""

		# load SPEC-file
		print('Parsing EDF- and SPEC-files of Scan No. %s' % scanNumber)
		fn = self.path + self.filename
		if use_PyMca == True:
			data, motors, counters = xrs_fileIO.PyMcaSpecRead(fn,scanNumber)
			edfmats = xrs_fileIO.ReadEdfImages_PyMca(counters['ccdno'], self.path, self.EDF_PREFIX, self.edfName, self.EDF_POSTFIX)
		else:
			data, motors, counters = xrs_fileIO.SpecRead(fn,scanNumber)
			edfmats = xrs_fileIO.ReadEdfImages_my(counters['ccdno'], self.path, self.EDF_PREFIX, self.edfName, self.EDF_POSTFIX)

		# add the scannumber to self.scannumbers, if not already present
		if not scanNumber in self.scanNumbers:
			self.scanNumbers.append(scanNumber)
		return data, motors, counters, edfmats


	def loadscan(self,scanNumbers,scanType='generic',direct=False):
		""" **loadscan**
		Loads one or multiple scans. 

		Arguments:
		----------
		scanNumbers (int or list of ints): Scan number or list of scan numbers to be loaded.
		scanType (string): String describing the scan for later automatic stitching/interpolation. Few special types exist: elastic, long.
		direct (boolean): Keyword if EDF-files should be deleted or kept (default).
        """
		# make sure scannumbers are iterable (list)
		if not isinstance(scanNumbers,list):
			scannums = []
			scannums.append(scanNumbers)
		else:
			scannums = scanNumbers 

		# load all scans from the list
		for number in scannums:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number)
			# create an instance of "scan" class for every scan
			onescan = xrs_scans.scan(edfmats,number,counters[self.scanMotor],counters[self.moniColumn],counters,motors,data,scanType)

			# applyrois and delete EDF-files if direct==True
			if direct:
				onescan.applyrois(self.roi_obj.indices,scaling=scaling)
				print('Deleting EDF-files of Scan No. %03d' % number)
				onescan.edfmats = [] # delete the edfmats
				self.scans[scanname] = onescan
			else:
				self.scans[scanname] = onescan

	def loadRockingCurve(self,scanNumbers,energyCoor=[9,3],direct=False):
		""" **loadRockingCurve**
		Load one or more rocking curves.

		Arguments:
		----------
		scanNumbers (int or list of ints): Scan number or list of scan numbers of rocking curve scans to be loaded.
		scanType (string): String describing the scan for later automatic stitching/interpolation. Few special types exist: elastic, long.
		energyCoor (list): Indices to find the energy during the scan based on the SPEC-file header.
		direct (boolean): Keyword if EDF-files should be deleted or kept (default).
		"""
		# make sure scanNumbers are iterable (list)
		if not isinstance(scanNumbers,list):
			scannums = []
			scannums.append(scanNumbers)
		else:
			scannums = scanNumbers 

		for number in scannums:
			self.loadscan(number,scanType='RC',direct=direct)
			scanname = 'Scan%03d' % number
			energy   = self.scans[scanname].motors[energyCoor[0]][energyCoor[1]]
			self.scans[scanname].offdia_energy = energy

	def stitchRockingCurves(self):
		""" **stitchRockingCurves**
		Go through all rocking curves and stitch them together to a 3D matrix.
		"""
		RcScans = xrs_scans.findRCscans(self.scans)
		sorted(RcScans,key=lambda x:x.offdia_energy)


	def set_roiObj(self,roiobj):
		""" **set_roiObj**
		Assigns an object of the roi_obj class to this class.
		"""
		self.roi_obj = roiobj

	def getrawdata(self):
		""" **getrawdata**
		Iterates through all instances of the scan class and calls it's applyrois method
		to sum up over all rois.
		"""
		if not np.any(self.roi_obj.indices):
			print 'Please define some ROIs first.'
			return
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ("integrating "+scan)
				self.scans[scan].applyrois(self.roi_obj.indices)

	def getrawdata_pixelwise(self):
		"""
		Goes through all instances of the scan class and calls it's applyrois_pw method
		to extract intensities for all pixels in each ROI.
		"""
		if not np.any(self.roi_obj.indices):
			print 'Please define some ROIs first.'
			return
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ("integrating pixelwise "+scan)
				self.scans[scan].applyrois_pw(self.roi_obj.indices)
     
        def SumDirect(self,scannumbers):
			""" **SumDirect**
			"""
			sum = None
			for number in scannumbers:
				data, motors, counters, edfmats = self.readscan(number)
				if sum is None:
					sum = np.zeros(edfmats[0].shape ,"f") 
					sum[:] += edfmats.sum(axis=0)
			return sum

	def deletescan(self,scannumbers):
		""" **deletescan**
		Deletes scans from the class.
		INPUT:
		scannumbers = integer or list of integers (SPEC scan numbers) to delete
		"""
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers
		for number in numbers:
			scanname = 'Scan%03d' % number
			del(self.scans[scanname])
			self.scannumbers.remove(number)

	def save_raw_data(self,filename):
		data = np.zeros((len(self.eloss),len(self.signals[0,:])))
		data[:,0]   = self.eloss
		data[:,1::] = self.signals
		np.savetxt(filename,data)

