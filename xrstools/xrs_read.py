#!/usr/bin/python
# Filename: xrs_read.py

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
import xrs_rois, xrs_scans, xrs_utilities, math_functions

from numpy import array
import scipy.io

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

import matplotlib.pyplot as plt

__metaclass__ = type # new style classes

class read_id20:
	"""
	Main class for handling raw data from XRS experiments on ESRF's ID20. This class
    is used to read scans from SPEC files and the according EDF-files, it provides access
    to all tools from the xrs_rois module for defining ROIs, it can be used to integrate
    scans, sum them up, stitch them together, and define the energy loss scale.
    INPUT:
    absfilename   = path and filename of the SPEC-file
    energycolumn  = name (string) of the counter for the energy as defined in the SPEC session (counter mnemonic)
    monitorcolumn = name (string) of the counter for the monitor signals as defined in the SPEC session (counter mnemonic)
    edfName       = name/prefix (string) of the EDF-files (default is the same as the SPEC-file name)
    single_image  = boolean switch, 'True' (default) if all 6 detectors are merged in a single image,
                    'False' if two detector images per point exist.
	"""
	def __init__(self,absfilename,energycolumn='energy',monitorcolumn='kap4dio',edfName=None,single_image=True):
		self.scans         = {} # was a dictionary before
		if not os.path.isfile(absfilename):
			raise Exception('IOError! No such file, please check filename.')
		self.path          = os.path.split(absfilename)[0] + '/'
		self.filename = os.path.split(absfilename)[1]
		if not edfName:
			self.edfName = os.path.split(absfilename)[1]
		else:
			self.edfName = edfName
		self.single_image  = single_image
		self.scannumbers   = []
		self.EDF_PREFIXh   = 'edf/h_'
		self.EDF_PREFIXv   = 'edf/v_'
		self.EDF_PREFIX    = 'edf/'
		self.EDF_POSTFIX   = '.edf'
		self.DET_PIXEL_NUMx = 768
		self.DET_PIXEL_NUMy = 512
		self.DET_PIXEL_NUM  = 256
		# which column in the SPEC file to be used for the energy and monitor
		self.encolumn      = energycolumn.lower()
		self.monicolumn    = monitorcolumn.lower()
		# here are the attributes of the old rawdata class
		self.eloss    = []	# common eloss scale for all analyzers
		self.energy   = []	# common energy scale for all analyzers
		self.signals  = []	# signals for all analyzers
		self.errors   = []	# poisson errors
		self.qvalues  = []	# for all analyzers
		self.groups   = {}  # dictionary of groups (such as 2 'elastic', or 5 'edge1', etc.)
		self.cenom    = []  # list of center of masses of the elastic lines
		self.E0       = []  # energy value, mean value of self.cenom
		self.tth      = []  # list of scattering angles (one for each ROI)
		self.VDtth    = []
		self.VUtth    = []
		self.VBtth    = []
		self.HRtth    = []
		self.HLtth    = []
		self.HBtth    = []
		self.resolution   = []  # list of FWHM of the elastic lines for each analyzer
		self.signals_orig = []  # signals for all analyzers before interpolation
		self.errors_orig  = []  # errors for all analyzers before interpolation
		# TTH offsets from center of V and H modules
		# tth offsets in one direction (horizontal for V-boxes, vertical for H-boxes)
		self.TTH_OFFSETS1   = np.array([5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0])
		# tth offsets in the other direction (horizontal for H-boxes, vertical for V-boxes)
		self.TTH_OFFSETS2   = np.array([-9.71, -9.75, -9.71, -3.24, -3.25, -3.24, 3.24, 3.25, 3.24, 9.71, 9.75, 9.71]) 
		# input
		self.rois  = []	# object of a container class from the helpers module (old)
		self.roi_obj = [] # an instance of the roi_object class from the xrs_rois module (new)

	def readscan(self,scannumber,fromtofile=False):
		"""
		Returns the data, motors, counter-names, and edf-files from the SPEC file defined when
        the xrs_read object was initiated.
		There should be an alternative that uses the PyMca module if installed.
        INPUT:
        scannumber = number of the scan to be loaded
        fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
		"""
		# first see if scan can be loaded from npz-file
		if fromtofile:
			try:
				print 'Trying to load scan from file.'
				scanname = 'Scan%03d' % scannumber
				sub_path = os.path.split(self.path[:-1])[0]
				fname    = sub_path + '/scans/' + scanname + '.npz'
				scan     = np.load(fname)
				data     = list(scan['data'])
				motors   = list(scan['motors'])
				counters = scan['counters'].item()
				edfmats  = scan['edfmats']
				return data, motors, counters, edfmats
			except:
				print 'Failed loading scan from file, will read EDF- and SPEC-file.'
				pass

		# proceed with parsing edf- and SPEC-files, if scan could not be loaded from zip-archive
		if not fromtofile:
			print 'Parsing EDF- and SPEC-files of scan No. %s' % scannumber

		# load SPEC-file
		fn = self.path + self.filename
		data, motors, counters = xrs_utilities.specread(fn,scannumber)

		if not self.single_image:
			# initiate arrays for the edf-files
			edfmatsh = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy/2,self.DET_PIXEL_NUMx)))
			edfmatsv = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy/2,self.DET_PIXEL_NUMx)))
			edfmats  = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy,self.DET_PIXEL_NUMx)))
			# load edf-files
			for m in range(len(counters['ccdno'])):
				ccdnumber = counters['ccdno'][m]
				edfnameh   = self.path + self.EDF_PREFIXh + self.edfName + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
				edfnamev   = self.path + self.EDF_PREFIXv + self.edfName + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
				edfmatsh[m,:,:] = (xrs_utilities.edfread_test(edfnameh))
				edfmatsv[m,:,:] = (xrs_utilities.edfread_test(edfnamev))

				edfmats[m,0:self.DET_PIXEL_NUMy/2,:] = edfmatsv[m,:,:]
				edfmats[m,self.DET_PIXEL_NUMy/2:,:]  = edfmatsh[m,:,:]
		else:
			# initiate arrays for the edf-files
			edfmats  = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy,self.DET_PIXEL_NUMx)))
			# load edf-files
			for m in range(len(counters['ccdno'])):
				ccdnumber = counters['ccdno'][m]
				edfname   = self.path + self.EDF_PREFIX + self.edfName + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
				edfmats[m,:,:] = (xrs_utilities.edfread_test(edfname))

		# add the scannumber to self.scannumbers, if not already present
		if not scannumber in self.scannumbers:
			self.scannumbers.extend([scannumber])

		# store scan in numpy zip-archive, if desired
		if fromtofile:
			scanname = 'Scan%03d' % scannumber
			sub_path = os.path.split(self.path[:-1])[0]
			fname = sub_path + '/scans/' + scanname
			if os.path.exists(sub_path):
				print 'trying to save file in numpy-archive.'
				np.savez(fname, data=data, motors=motors, counters=counters, edfmats=edfmats)
			else:
				print 'please create a folder ' + fname + ' first!'
				pass

		return data, motors, counters, edfmats #<type 'list'> <type 'list'> <type 'dict'> <type 'numpy.ndarray'>

	def loadscan(self,scannumbers,scantype='generic',fromtofile=False):
		"""
		Loads the files belonging to scan No. "scannumber" and puts it into an instance
		of the xrs_scan-class 'scan'. The default scantype is 'generic', later the scans
		will be grouped (and added) based on the scantype.
        INPUT:
        scannumbers = integer or list of scannumbers that should be loaded
        scantype    = string describing the scan to be loaded (e.g. 'edge1' or 'K-edge') 
        fromtofile  = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
        """
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers 
		for number in scannums:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number,fromtofile)
			# can assign some things here already (even if maybe redundant)
			monitor   = counters[self.monicolumn]
			monoangle = 1 # have to check for this later ... !!! counters['pmonoa']
			energy    = counters[self.encolumn]
			# create an instance of "scan" class for every scan
			onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
			# assign one dictionary entry to each scan 
			self.scans[scanname] = onescan

	def loadloop(self,begnums,numofregions,fromtofile=False):
		"""
		Loads a whole loop of scans based on their starting scannumbers and the number of single
		scans in the loop.
        INPUT:
		begnums      = list of scannumbers of the first scans of each loop (is a list)
		numofregions = number of scans in each loop (integer)
		"""
		typenames = []
		for n in range(numofregions):
			typenames.append('edge'+str(n+1))
		numbers = []
		for n in range(len(begnums)):
			for m in range(numofregions):
				numbers.append((begnums[n]+m))
		typenames = typenames*len(begnums)
		for n in range(len(typenames)):
			thenumber = []
			thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
			self.loadscan(thenumber,typenames[n],fromtofile)

	def loadelastic(self,scann,fromtofile=False):
		"""
		Loads a scan using the loadscan function and sets the scantype attribute to 'elastic'.
        I.e. shorthand for 'obj.loadscan(scannumber,type='elastic')'.
        INPUT:
        scann      = integer or list of integers
        fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)        
		"""
		self.loadscan(scann,'elastic',fromtofile)
	
	def loadlong(self,scann,fromtofile=False):
		"""
		Loads a scan using the loadscan function and sets the scantype attribute to 'long'.
        I.e. shorthand for 'obj.loadscan(scannumber,type='long')'.
        INPUT:
        scann      = integer or list of integers
        fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
		"""
		self.loadscan(scann,'long',fromtofile)

        
	def get_zoom_rois(self,scannumbers,logscaling=True,colormap='jet',interpolation='nearest'):
		"""
		Define ROIs by zooming into an image constructed from the sum of all edf-files 
		in 'scannumbers'
		scannumbers = either single scannumber or list of scannumbers
		logscaling  = set to 'True' if images is to be shown on log-scale (default is True)
		colormap    = string to define the colormap which is to be used for display (anything 
		              supported by matplotlib, 'jet' by default)
		interpolation = interpolation scheme to be used for displaying the image (anything
		                supported by matplotlib, 'nearest' by default)
		"""
		image = xrs_rois.create_sum_image(self.scans,scannumbers)
		roi_finder_obj = xrs_rois.roi_finder()
		roi_finder_obj.get_zoom_rois(image,logscaling=logscaling,colormap=colormap,interpolation=interpolation)
		self.roi_obj = roi_finder_obj.roi_obj

	def get_linear_rois(self,scannumbers,logscaling=True,height=5,colormap='jet',interpolation='nearest'):
		"""
		Define ROIs by zooming into an image constructed from the sum of all edf-files 
		in 'scannumbers'
		scannumbers = either single scannumber or list of scannumbers
		logscaling  = set to 'True' if images is to be shown on log-scale (default is True)
		colormap    = string to define the colormap which is to be used for display (anything 
		              supported by matplotlib, 'jet' by default)
		interpolation = interpolation scheme to be used for displaying the image (anything
		                supported by matplotlib, 'nearest' by default)
		"""
		image = xrs_rois.create_sum_image(self.scans,scannumbers)
		roi_finder_obj = xrs_rois.roi_finder()
		roi_finder_obj.get_linear_rois(image,logscaling=logscaling,height=height,colormap=colormap,interpolation=interpolation)
		self.roi_obj = roi_finder_obj.roi_obj

	def get_auto_rois(self,scannumbers,kernel_size=5,threshold=100.0,logscaling=True,colormap='jet',interpolation='bilinear'):
		"""
		Define ROIs automatically using median filtering and a variable threshold.
		scannumbers   = either single scannumber or list of scannumbers
		kernel_size   = used kernel size for the median filter (must be an odd integer)
		logscaling    = set to 'True' if images is to be shown on log-scale (default is True)
		colormap      = string to define the colormap which is to be used for display (anything 
		                supported by matplotlib, 'jet' by default)
		interpolation = interpolation scheme to be used for displaying the image (anything
		                supported by matplotlib, 'nearest' by default)
		"""
		# check that kernel_size is odd
		if not kernel_size % 2 == 1:
			print 'The \'kernal_size\' must be an odd number.'
			return
		image = xrs_rois.create_sum_image(self.scans,scannumbers)
		roi_finder_obj = xrs_rois.roi_finder()
		roi_finder_obj.get_auto_rois(image,kernel_size=kernel_size,threshold=threshold,logscaling=logscaling,colormap=colormap,interpolation=interpolation)
		self.roi_obj = roi_finder_obj.roi_obj


	def get_auto_rois_eachdet(self,scannumbers,kernel_size=5,threshold=100.0,logscaling=True,colormap='jet',interpolation='bilinear'):
		"""
		Define ROIs automatically using median filtering and a variable threshold for each detector
		separately.
		scannumbers   = either single scannumber or list of scannumbers
		kernel_size   = used kernel size for the median filter (must be an odd integer)
		logscaling    = set to 'True' if images is to be shown on log-scale (default is True)
		colormap      = string to define the colormap which is to be used for display (anything 
		                supported by matplotlib, 'jet' by default)
		interpolation = interpolation scheme to be used for displaying the image (anything
		                supported by matplotlib, 'nearest' by default)
		"""
		# check that kernel_size is odd
		if not kernel_size % 2 == 1:
			print 'The \'kernal_size\' must be an odd number.'
			return

		# create a big image
		image = xrs_rois.create_sum_image(self.scans,scannumbers)

		# break down the image into 256x256 pixel images
		det_images, offsets = xrs_rois.break_down_det_image(image,self.DET_PIXEL_NUM)

		# create one roi_object per sub-image
		temp_objs = []
		for ii in range(det_images.shape[0]):
			temp = xrs_rois.roi_finder()
			temp.get_auto_rois(det_images[ii,:,:],kernel_size=kernel_size,threshold=threshold,logscaling=logscaling,colormap=colormap,interpolation=interpolation)
			temp_objs.append(temp)

		# merge all roi_objects into one
		merged_obj   = xrs_rois.merge_roi_objects_by_matrix(temp_objs,image.shape,offsets,self.DET_PIXEL_NUM)
		self.roi_obj = merged_obj

	def get_polygon_rois(self,scannumbers,logscaling=True,colormap='jet',interpolation='nearest'):
		"""
		Define a polygon shaped ROI from an image constructed from 
		the sum of all edf-files in 'scannumbers'
		image_shape = tuple with shape of the current image (i.e. (256,256))
		scannumbers = either single scannumber or list of scannumbers
		logscaling  = set to 'True' if images is to be shown on log-scale (default is True)
		colormap    = string to define the colormap which is to be used for display (anything 
		              supported by matplotlib, 'jet' by default)
		interpolation = interpolation scheme to be used for displaying the image (anything
		                supported by matplotlib, 'nearest' by default)
		"""
		image = xrs_rois.create_sum_image(self.scans,scannumbers)
		roi_finder_obj = xrs_rois.roi_finder()
		roi_finder_obj.get_polygon_rois(image,logscaling=logscaling,colormap=colormap,interpolation=interpolation)
		self.roi_obj = roi_finder_obj.roi_obj

	def get_polygon_rois_eachdet(self,scannumbers,logscaling=True,colormap='jet',interpolation='nearest'):
		"""
		Define a polygon shaped ROI from an image constructed from 
		the sum of all edf-files in 'scannumbers'
		image_shape = tuple with shape of the current image (i.e. (256,256))
		scannumbers = either single scannumber or list of scannumbers
		logscaling  = set to 'True' if images is to be shown on log-scale (default is True)
		colormap    = string to define the colormap which is to be used for display (anything 
		              supported by matplotlib, 'jet' by default)
		interpolation = interpolation scheme to be used for displaying the image (anything
		                supported by matplotlib, 'nearest' by default)
		"""

		# create a big image
		image = xrs_rois.create_sum_image(self.scans,scannumbers)

		# break down the image into 256x256 pixel images
		det_images, offsets = xrs_rois.break_down_det_image(image,self.DET_PIXEL_NUM)

		# create one roi_object per sub-image
		temp_objs = []
		for modind in range(det_images.shape[0]):
			temp = xrs_rois.roi_finder()
			temp.get_polygon_rois( det_images[modind,:,:],modind,logscaling=logscaling,colormap=colormap,interpolation=interpolation)
			temp_objs.append(temp)

		# merge all roi_objects into one
		merged_obj   = xrs_rois.merge_roi_objects_by_matrix(temp_objs,image.shape,offsets,self.DET_PIXEL_NUM)
		self.roi_obj = merged_obj

	def save_rois(self,filename):
		"""
		Saves the defined ROIs in a pick
		"""
		import pickle
		f = open(filename, 'wb')
		pickle.dump(self.roi_obj, f,protocol=-1)
		f.close()

	def loadrois(self,filename):
		"""
		loads a file written with the saverois-function using pickle
		filename = absolute path to file and filename
		"""
		import pickle
		f = open(filename,'rb')
		self.rois = pickle.load(f)
		self.roi_obj = xrs_rois.roi_object()
		self.roi_obj.indices = xrs_rois.swap_indices_old_rois(self.rois.inds)
		f.close()

	def load_rois(self,filename):
		"""
		Loads ROIs from a file written with the save_rois-function.
		filename = absolute path to file and filename
		"""
		import pickle
		f = open(filename,'rb')
		self.roi_obj = pickle.load(f)
		f.close()

	def set_roiObj(self,roiobj):
		self.roi_obj = roiobj
		self.roi_obj.load_rois_fromMasksDict(masksDict)

	def show_rois(self):
		"""
		Produces an image of all selected ROIs and their according number.
		"""
		# check if ROIs are defined
		if not self.roi_obj:
			print 'Please define some ROIs first.'
			return
		else:
			self.roi_obj.show_rois()		

	def findroisColumns(self,scannumbers,whichroi,logscaling=False):
		"""
		Constructs a waterfall plot from the columns in a scan, i.e. energy vs. pixel number along the ROI
		scannumbers = scannumber or list of scannumbers from which to construct the plot
		whichroi    = integer (starting from 0) from which ROI to use for constructing the plot
		"""
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers 

		if not self.roi_obj.indices:
			'Please define some zoom ROIs first.'
			return
		if not self.roi_obj.kind == 'zoom':
			'Currently this feature only works for ROIs of type \'zoom\'.'
			return

		xinds     = np.unique(self.roi_obj.x_indices[whichroi])
		yinds     = np.unique(self.roi_obj.y_indices[whichroi])
		scanname  = 'Scan%03d' % scannums[0]
		edfmats   = np.zeros_like(self.scans[scanname].edfmats)
		energy    = self.scans[scanname].energy
		waterfall = np.zeros((len(yinds),len(energy)))

		for scannum in scannums:
			scanname = 'Scan%03d' % scannum
			scanedf  = self.scans[scanname].edfmats
			scanmonitor = self.scans[scanname].monitor
			for ii in range(len(energy)):
				edfmats[ii,:,:] += scanedf[ii,:,:]/scanmonitor[ii]

		for ii in range(len(energy)):
			for jj in range(len(yinds)):
				waterfall[jj,ii] = np.sum(edfmats[ii,xinds,yinds[jj]])

		plt.figure()
		for ii in range(len(yinds)):
			plt.plot(waterfall[ii,:])

		fig = plt.figure()
		ax = fig.add_subplot(111)

		if logscaling:
			ax.imshow(np.log(np.transpose(waterfall)), interpolation='nearest')
		else:
			ax.imshow(np.transpose(waterfall), interpolation='nearest')

		ax.set_aspect('auto')
		plt.xlabel('ROI pixel')
		plt.ylabel('energy point')
		plt.show()

	def orderrois(self,arrangement='vertical',missing=None):
		"""
		order the rois in an order provided such that e.g. autorois have the correct order 
		"""
		if not self.roi_obj:
			print 'Please select some ROIs first!'
			return

		# get some detector infos
		all_det_names = ['VD','VU','VB','HR','HL','HB']
		all_dets      = []
		for name in all_det_names:
			all_dets.append(xrs_utilities.maxipix_det(name,arrangement))

		# find the ROIs for each detector
		# go through all detectors and ROIs defined and see which one has centers in which detector
		for det in all_dets:
			det_x_min = det.get_pixel_range()[0]
			det_x_max = det.get_pixel_range()[1]
			det_y_min = det.get_pixel_range()[2]
			det_y_max = det.get_pixel_range()[3]
			for ii in range(len(self.roi_obj.indices)):
				x_mean = np.mean(self.roi_obj.x_indices[ii])
				y_mean = np.mean(self.roi_obj.y_indices[ii])
				if x_mean >= det_x_min and x_mean <= det_x_max and y_mean >= det_y_min and y_mean <= det_y_max:
					det.roi_indices.append(self.roi_obj.indices[ii])
					det.roi_x_indices.append(self.roi_obj.x_indices[ii])
					det.roi_y_indices.append(self.roi_obj.y_indices[ii])
					det.roi_x_means.append(x_mean)
					det.roi_y_means.append(y_mean)
	
		# count, check with missing, break if something is wrong
		for det in all_dets:
			if not len(det.roi_indices) == 12:
				print 'WARNING! Module ' + det.name + ' only has ' + '%d' % len(det.roi_indices) + ' ROIs defined, your numbering will be messed up!'
		# rearrange according to 'arrangement' keyword
		if arrangement == 'vertical':
			verticalIndex = [0,3,6,9,1,4,7,10,2,5,8,11]	
			for det in all_dets:
				# order from low to high y-mean
				roi_coords   = np.array(det.roi_indices)
				roi_y_means  = np.array(det.roi_x_means) 
				sorting_inds = roi_y_means.argsort()
				roi_coords_increasing = [roi_coords[i] for i in sorting_inds]
				if det.get_det_name() in ['VD','VU','VB']:
					# sort from high to low y-center
					det.roi_indices = [roi_coords_increasing[i] for i in verticalIndex[::-1]]
				elif det.get_det_name() in ['HR','HL','HB']:
					# sort from low to high y-center
					det.roi_indices = [roi_coords_increasing[i] for i in verticalIndex]
				else:
					print 'Sorry, no such module.'

		# reassign all roi_obj variables... 
		allrois = []
		for det in all_dets:
			allrois.extend(det.roi_indices)

		self.roi_obj.indices = allrois
		self.roi_obj.roi_matrix     = xrs_rois.convert_inds_to_matrix(self.roi_obj.indices,(512,768))
		self.roi_obj.red_rois       = xrs_rois.convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.x_indices      = xrs_rois.convert_inds_to_xinds(self.roi_obj.indices)
		self.roi_obj.y_indices      = xrs_rois.convert_inds_to_yinds(self.roi_obj.indices)
		self.roi_obj.masks          = xrs_rois.convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def getrawdata(self):
		"""
		Goes through all instances of the scan class and calls it's applyrois method
		to sum up over all rois.
		"""
		if not np.any(self.roi_obj.indices):
			print 'Please define some ROIs first.'
			return
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ("integrating "+scan)
				self.scans[scan].applyrois(self.roi_obj.indices)


        def SumDirect(self,scannumbers):
                sum=None
		for number in scannumbers:
			data, motors, counters, edfmats = self.readscan(number)
                        if sum is None:
                                sum = np.zeros(edfmats[0].shape ,"f") 
                        sum[:] += edfmats.sum(axis=0) 
                return sum

	def loadscandirect(self,scannumbers,scantype='generic',fromtofile=False,scaling=None):
		"""
		Loads a scan without saving the edf files in matrices.
		scannumbers = integer or list of integers defining the scannumbers from the SPEC file
		scantype    = string describing the scan to be loaded (e.g. 'edge1' or 'K-edge') 
        fromtofile  = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
		scaling     = list of scaling factors to be applied, one for each ROI defined
		"""
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers 
		# check if there are ROIs defined		
		if not self.roi_obj:
			print 'Please define some ROIs first'
			return
		for number in scannums:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number,fromtofile)
			# can assign some things here already (even if maybe redundant)
			monitor   = counters[self.monicolumn]
			monoangle = 1 # counters['pmonoa'] # this still needs checking
			energy    = counters[self.encolumn]
			# create an instance of "scan" class for every scan

			onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)

			onescan.applyrois(self.roi_obj.indices,scaling=scaling)

			print 'Deleting EDF-files of Scan No. %03d' % number
			onescan.edfmats = [] # delete the edfmats
			self.scans[scanname] = onescan

	def loadloopdirect(self,begnums,numofregions,fromtofile=False,scaling=None):
		"""
		Loads a whole loop of scans based on their starting scannumbers and the number of single
        INPUT:
		begnums      = list of scannumbers of the first scans of each loop (is a list)
		numofregions = number of scans in each loop (integer)
        fromtofile   = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
		scaling      = list of scaling factors to be applied, one for each ROI defined
		"""
		typenames = []
		for n in range(numofregions):
			typenames.append('edge'+str(n+1))
		numbers = []
		for n in range(len(begnums)):
			for m in range(numofregions):
				numbers.append((begnums[n]+m))
		typenames = typenames*len(begnums)
		for n in range(len(typenames)):
			thenumber = []
			thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
			self.loadscandirect(thenumber,typenames[n],fromtofile,scaling=scaling)

	def loadelasticdirect(self,scann,fromtofile=False):
		"""
		Loads a scan using the loadscan function and sets the scantype attribute to 'elastic'.
        I.e. shorthand for 'obj.loadscan(scannumber,type='elastic')'.
        INPUT:
        scann      = integer or list of integers
        fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)        
		"""
		self.loadscandirect(scann,'elastic',fromtofile)

	def loadlongdirect(self,scann,fromtofile=False,scaling=None):
		"""
		Loads a scan using the loadscan function and sets the scantype attribute to 'long'.
        I.e. shorthand for 'obj.loadscan(scannumber,type='long')'.
        INPUT:
        scann      = integer or list of integers
        fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)        
		"""
		self.loadscandirect(scann,'long',fromtofile,scaling=scaling)

	def deletescan(self,scannumbers):
		"""
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

	def getspectrum(self, include_elastic=False):
		"""
		Groups the instances of the scan class by their scantype attribute, 
		adds equal scans (each group of equal scans) and appends them.
		INPUT:
		include_elastic = boolean flag, skips the elastic line if set to 'False' (default)
		"""
		# find the groups 
		allgroups = xrs_scans.findgroups(self.scans)
		for group in allgroups:
			# self.makegroup(group)
			onegroup = xrs_scans.makegroup_nointerp(group)
			self.groups[onegroup.get_type()] = onegroup
		self.energy,self.signals,self.errors = xrs_scans.appendScans(self.groups,include_elastic)
		self.signals_orig = self.signals
		self.errors_orig  = self.errors

	def geteloss(self):
		"""
		Defines the energy loss scale for all ROIs by finding the center of mass for each ROI's elastic line.
		Interpolates the signals and errors onto a commom energy loss scale. Finds the resolution (FWHM) of the
		'elastic' groups.
		"""
		if not 'elastic' in self.groups:
			print 'Please load/integrate at least one elastic scan first!'
			return
		else:
			# reset values, in case this is run several times
			self.cenom = []
			self.resolution = []
			self.E0 = []
			for n in range(len(self.roi_obj.indices)):
				# find the center of mass for each ROI
				self.cenom.append(xrs_utilities.find_center_of_mass(self.groups['elastic'].energy,self.groups['elastic'].signals_orig[:,n]))
				# try finden the FWHM/resolution for each ROI
				try:
					FWHM,x0 = fwhm((self.groups['elastic'].energy - self.cenom[n])*1e3,self.groups['elastic'].signals_orig[:,n])
					self.resolution.append(FWHM)
				# append a zero if the FWHM routine fails
				except:
					self.resolution.append(0.0) # need a more sofisticated way of finding the FWHM
			self.E0 = np.mean(self.cenom)
			# define the first eloss scale as the 'master' scale for all ROIs
			self.eloss = (self.energy - self.cenom[0])*1e3 # energy loss in eV
			# define eloss-scale for each ROI and interpolate onto the 'master' eloss-scale 
			for n in range(len(self.roi_obj.indices)):
				# inserting zeros at beginning and end of the vectors to avoid interpolation errors
				x = (self.energy-self.cenom[n])*1e3
				x = np.insert(x,0,-1e10)
				x = np.insert(x,-1,1e10)
				y = self.signals_orig[:,n]
				y = np.insert(y,0,0)
				y = np.insert(y,-1,0)
				f = interp1d(x,y, bounds_error=False,fill_value=0.0)
				self.signals[:,n] = f(self.eloss)
			# do the same for the errors 
			for n in range(len(self.roi_obj.indices)):
				# inserting zeros at beginning and end of the vectors to avoid interpolation errors
				x = (self.energy-self.cenom[n])*1e3
				x = np.insert(x,0,-1e10)
				x = np.insert(x,-1,1e10)
				y = self.errors_orig[:,n]
				y = np.insert(y,0,0)
				y = np.insert(y,-1,0)
				f = interp1d(x,y, bounds_error=False,fill_value=0.0)
				self.errors[:,n] = f(self.eloss)

	def gettths(self,rvd=0.0,rvu=0.0,rvb=0.0,rhl=0.0,rhr=0.0,rhb=0.0,order=[0,1,2,3,4,5]):
		"""
		Uses the defined TT_OFFSETS of the read_id20 class to set all scattering angles tth
		from the mean angle avtth of the analyzer modules.
		INPUT:
		rhl = mean tth angle of HL module (default is 0.0)
		rhr = mean tth angle of HR module (default is 0.0)
		rhb = mean tth angle of HB module (default is 0.0)
		rvd = mean tth angle of VD module (default is 0.0)
		rvu = mean tth angle of VU module (default is 0.0)
		rvb = mean tth angle of VB module (default is 0.0)
		order = list of integers (0-5) which describes the order of modules in which the 
				ROIs were defined (default is VD, VU, VB, HR, HL, HB; i.e. [0,1,2,3,4,5])
		"""
		# reset all values, just in case mean angles are redefined (values are otherwise appended to existing values)
		self.VDtth = []
		self.VUtth = []
		self.VBtth = []
		self.HRtth = []
		self.HLtth = []
		self.HBtth = []
		self.tth   = []
		# horizontal modules
		# HL (motor name rhl)
		v_angles = self.TTH_OFFSETS1 
		h_angles = self.TTH_OFFSETS2 + rhl
		for n in range(len(h_angles)):
			self.HLtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)
		# HR (motor name rhr)
		v_angles = self.TTH_OFFSETS1 
		h_angles = self.TTH_OFFSETS2 + rhr
		for n in range(len(h_angles)):
			self.HRtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)
		# HB (motor name rhb)
		v_angles = self.TTH_OFFSETS1 
		h_angles = self.TTH_OFFSETS2 + rhb
		for n in range(len(h_angles)):
			self.HBtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)

		# vertical modules
		# VD
		v_angles = self.TTH_OFFSETS2 + rvd
		h_angles = self.TTH_OFFSETS1 
		for n in range(len(h_angles)):
			self.VDtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)
		# VU
		v_angles = self.TTH_OFFSETS2 + rvu
		h_angles = self.TTH_OFFSETS1 
		for n in range(len(h_angles)):
			self.VUtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)
		# VB
		v_angles = self.TTH_OFFSETS2 + rvb
		h_angles = self.TTH_OFFSETS1 
		for n in range(len(h_angles)):
			self.VBtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)

		# list of TTH values
		tth = [self.VDtth, self.VUtth, self.VBtth, self.HRtth, self.HLtth, self.HBtth]
		# list all TTH values in one long list ordered by the 'order'-keyword
		for n in order:
			self.tth.extend(tth[n])

	def getqvals(self,invangstr=False):
		"""
		Calculates q-values from E0 and tth values in either atomic units (defalt) or
		inverse angstroms.
		"""
		theqvals = np.zeros_like(self.signals)
		if invangstr:
			for n in range(len(self.tth)):
				theqvals[:,n] = xrs_utilities.momtrans_inva(self.E0+self.eloss/1e3,self.E0,self.tth[n]) 
		else:
			for n in range(len(self.tth)):
				theqvals[:,n] = xrs_utilities.momtrans_au(self.E0+self.eloss/1e3,self.E0,self.tth[n])
		self.qvalues = theqvals

	def getqvals_energy(self,energy):
		"""
		Returns all q-values at a certain energy loss.
		INPUT:
		energy = energy loss value for which all q-values are stored
		"""
		ind = np.abs(self.eloss - energy).argmin()
		return self.qvalues[ind,:]

	def copy_edf_files(self,scannumbers,destdir):
		"""
		Copies all edf-files from scan with scannumber or scannumbers into directory 'destdir'
		INPUT:
		scannumbers = integer or list of integers defining the scannumbers from the SPEC file
		destdir     = string with absolute path for the destination
		"""
		import shutil
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers
		fn = self.path + self.filename
		if not self.single_image:
			for n in range(len(numbers)):
				data, motors, counters = xrs_utilities.specread(fn,numbers[n])
				for m in range(len(counters['ccdno'])):
					ccdnumber = counters['ccdno'][m]
					edfnameh   = self.path + self.EDF_PREFIXh + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
					edfnamev   = self.path + self.EDF_PREFIXv + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
					shutil.copy2(edfnameh, destdir)
					shutil.copy2(edfnamev, destdir)
		if self.single_image:
			for n in range(len(numbers)):
				data, motors, counters = xrs_utilities.specread(fn,numbers[n])
				for m in range(len(counters['ccdno'])):
					ccdnumber = counters['ccdno'][m]
					edfname   = self.path + self.EDF_PREFIX + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
					shutil.copy2(edfname, destdir)

	def make_posscan_image(self,scannumber,motorname,filename=None):
		"""
		Loads a scan from a sample position scan (x-scan, y-scan, z-scan), lets you choose a zoomroi and constructs a 2D image from this
		INPUT:
		scannumber = number of the scan
		motorname  = string that contains the motor name (must be the same as in the SPEC file)
		filename   = optional parameter with filename to store the image
		"""
		plt.clf()
		# load the scan
		data, motors, counters, edfmats = self.readscan(scannumber)
		# the scan motor
		position = counters[motorname.lower()]
		# define a zoom ROI
		image = xrs_utilities.sumx(edfmats)
		roi_finder_obj = xrs_rois.roi_finder()
		roi_finder_obj.get_zoom_rois(image)
		# construct the image
		roixinds = roi_finder_obj.roi_obj.x_indices[0]
		roiyinds = roi_finder_obj.roi_obj.y_indices[0]
		# go through all edf files of the scan, sum over the height of the roi and stack the resulting lines into a matrix
		axesrange = [0,roiyinds[-1],position[-1],position[0]]
		theimage = (np.sum(edfmats[:,np.amin(roixinds):np.amax(roixinds)+1,np.amin(roiyinds):np.amax(roiyinds)+1],axis=1))
		plt.close()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.imshow(np.log(theimage),extent=axesrange)
		ax.set_aspect('auto')
		plt.xlabel('pixel along the beam')
		ylabelstr = motorname.lower() + ' position [mm]'
		plt.ylabel(ylabelstr)
		plt.show()
		# save the image, if a filename is provided
		if filename:
			f = open(filename, 'wb')
			yrange = np.arange(np.amin(roixinds),np.amax(roixinds)+1)
			theobject = LRimage(theimage,position,yrange)
			pickle.dump(theobject, f,protocol=-1)
			f.close()

	def animation(self,scannumber,logscaling=True,timeout=-1,colormap='jet'):
		"""
		Shows the edf-files of a scan as a 'movie'.
		INPUT:
		scannumber = integer/scannumber
		logscaling = set to 'True' (default) if edf-images are to be shown on logarithmic-scale
		timeout    = time in seconds defining pause between two images, if negative (default)
					 images are renewed by mouse clicks
		colormap   = matplotlib color scheme used in the display
		"""
		if isinstance(scannumber,list):
			if len(scannumber)>1:
				print 'this only works for a single scan, sorry'
				return
			else:
				scannumber = scannumber[0]
		scanname = 'Scan%03d' % scannumber
		edfmats = self.scans[scanname].edfmats
		scanlen  = np.shape(edfmats)[0]
		plt.ion()
		plt.clf()
		for n in range(scanlen):
			plt.clf()
			if logscaling:
				theimage = plt.imshow(np.log(edfmats[n,:,:]))
			else:
				theimage = plt.imshow(edfmats[n,:,:])
			plt.xlabel('detector x-axis [pixel]')
			plt.ylabel('detector y-axis [pixel]')
			if timeout<0:
				titlestring = 'Frame No. %d' % (n+1) + ' of %d' % scanlen + ', press key or mouse botton to continue' 
				plt.title(titlestring)
			else:
				titlestring = 'Frame No. %d' % (n+1) + ' of %d' % scanlen + ', updating every %2.2f ' % timeout + ' seconds'
				plt.title(titlestring)
			theimage.set_cmap(colormap)
			plt.draw()
			plt.waitforbuttonpress(timeout=timeout)

	def animationScans(self,scannumbers,logscaling=True,timeout=-1,colormap='jet'):
		"""
		shows a sum of all edf-files from a scan for one or several scans
		"""
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers

		for number in numbers:
			scanname = 'Scan%03d' % number
			edfmats  = self.scans[scanname].edfmats
			edfsum   = sumx(edfmats)
			if logscaling:
				theimage = plt.imshow(np.log(edfsum))
			else:
				theimage = plt.imshow(edfsum)
			plt.xlabel('detector x-axis [pixel]')
			plt.ylabel('detector y-axis [pixel]')
			if timeout<0:
				titlestring = 'Scan No. %3d' % (number) + ', press key or mouse botton to continue' 
				plt.title(titlestring)
			else:
				titlestring = 'Scan No. %3d' % (number) + ', updating every %2.2f ' % timeout + ' seconds'
				plt.title(titlestring)
			theimage.set_cmap(colormap)
			plt.draw()
			plt.waitforbuttonpress(timeout=timeout)

	def showscans(self,scannumbers,columns=None):
		"""
		plots spectra from a scan or series of scans
		"""
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers
		plt.ion()
		plt.clf()
		for number in numbers:
			scanname = 'Scan%03d' % number
			# show all individual ROIs if no columns are specified
			if not columns:
				y = np.array(self.scans[scanname].signals)
				for n in range(len(self.roi_obj.indices)):
					y[:,n] = y[:,n]/self.scans[scanname].monitor
				try:
					titlestring = 'all analyzer signals of scan No. ' + scanname
					plt.title(titlestring)
					plt.plot(self.scans[scanname].eloss,y)
					plt.xlabel('energy loss [eV]')
					plt.ylabel('signal [arb. units]')
				except:
					titlestring = 'all analyzer signals of scan No. ' + scanname
					plt.title(titlestring)
					plt.plot(self.scans[scanname].energy,y)
					plt.xlabel('primary energy [keV]')
					plt.ylabel('signal [arb. units]')
			# sum over the columns if they are provided
			if columns:
				if isinstance(columns,list):
					y = np.zeros_like(self.scans[scanname].energy)
					for n in columns:
						y += self.scans[scanname].signals[:,n]
					y = y/self.scans[scanname].monitor
				else:
					y = self.scans[scanname].signals[:,columns]
				try:
					titlestring = 'all analyzer signals of scan No. ' + scanname
					plt.title(titlestring)
					plt.plot(self.scans[scanname].eloss,y)
					plt.xlabel('energy loss [eV]')
					plt.ylabel('signal [arb. units]')
				except:
					titlestring = 'all analyzer signals of scan No. ' + scanname
					plt.title(titlestring)
					plt.plot(self.scans[scanname].energy,y)
					plt.xlabel('primary energy [keV]')
					plt.ylabel('signal [arb. units]')

	def showloops(self,scannumbers,looplength):
		"""
		plots the results one loop at a time
		"""
		pass

	def plotspectrum(self,columns=None):
		"""
		make a plot of the raw spectrum of the entire instance, sum over selected analyzers if porvided
		"""
		if len(self.eloss) == 0:
			print 'please provide elastic line scans first and run the geteloss() method'
			return
		if not columns:
			plt.plot(self.eloss,self.signals)
			plt.title('all analyzers vs. energy loss')
			plt.xlabel('energy loss [eV]')
			plt.ylabel('normalized signal [arb. units]')
		else:
			plt.plot(self.eloss,np.sum(self.signals[:,columns],axis=1))
			plt.title('all analyzers vs. energy loss')
			plt.xlabel('energy loss [eV]')
			plt.ylabel('normalized signal [arb. units]')

	def printlength(self,scannumbers):
		"""
		Prints the number of energy points in a scan or a number of scans.
		INPUT:
		scannumbers = integer or list of integers
		"""
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers

		for i in numbers:
			name = 'Scan%03d' % i
			print 'length of scan %03d ' %i + ' is ' + str(len(self.scans[name].energy))

	def removeBackgroundRoi(self,backroinum,estart=None,estop=None):
		if not estart:
			estart = self.eloss[0]
		if not estop:
			estop = self.elsos[-1]
		if not self.tth:
			print 'Please define the scattering angles first using the gettth method!'
			return

		for ii in range(len(self.tth)):
			if ii != backroinum:
				inds       = np.where(np.logical_and(self.eloss>=estart, mpc96.eloss<=estop))
				expnorm    = np.trapz(self.signals[inds,ii],self.eloss[inds])
				backnorm   = np.trapz(self.signals[inds,backroinum],self.eloss[inds])
				background = self.signals[:,backroinum]/backnorm*expnorm
				self.signals[:,ii] -= background # subtract background roi

class read_id16:
	"""
	a class for loading data form SPEC files and the according edf files from the old id16
	"""
	def __init__(self,absfilename,energycolumn='energy_cc',monitorcolumn='monitor'):
		self.scans         = {} # was a dictionary before
		if not os.path.isfile(absfilename):
			raise Exception('IOError! No such file, please check filename.')
		self.path          = os.path.split(absfilename)[0] + '/'
		self.filename      = os.path.split(absfilename)[1]
		self.scannumbers   = []
		self.EDF_PREFIX    = 'edf/'
		self.EDF_POSTFIX   = '.edf'
		self.DET_PIXEL_NUM = 256
		self.TTH_OFFSETS   = np.array([-13.0,-6.5,-6.5,0.0,0.0,6.5,6.5,13.0,13.0])
		# which column in the SPEC file to be used for the energy and monitor
		self.encolumn      = energycolumn.lower()
		self.monicolumn    = monitorcolumn.lower()
		# here are the attributes of the old rawdata class
		self.eloss    = []	# common eloss scale for all analyzers
		self.energy   = []	# common energy scale for all analyzers
		self.signals  = []	# signals for all analyzers
		self.errors   = []	# poisson errors
		self.qvalues  = []	# for all analyzers
		self.groups   = {}  # groups (such as 2 'elastic', or 5 'edge1', etc.)
		self.cenom    = []
		self.E0       = []
		self.tth      = []
		# input
		self.rois     = []	# can be set once the rois are defined to integrate scans directly
		# intermediates
		self.rawsignals = []	# structure of rawdata (energies, raw-signals, monitor, etc.)
		self.rawenergy  = []
		self.rawmonitor = []

	def readscan(self,scannumber):
		"""
		returns the data, motors, counter-names, and edf-files from the calss's specfile of "scannumber"
		should use PyMca's routines to read the Spec- and edf-files
		"""
		if self.path: # if a path was provided upon class construction, use this one
			path = self.path
		else: # if no path was provided
			print 'please provide a path to the SPEC-file'
			return
		if not self.filename:
			print 'please provide a SPEC-file name!'
			return
		else:
			fn = path + self.filename
		#try: #try loading variables from file
		print ("parsing edf- and SPEC-files of scan No. %s" % scannumber)
		data, motors, counters = xrs_utilities.specread(fn,scannumber)
		edfmats = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUM,self.DET_PIXEL_NUM)))
		for m in range(len(counters['ccdno'])):
			ccdnumber = counters['ccdno'][m]+1
			edfname   = path + self.EDF_PREFIX + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
			edfmats[m,:,:] = xrs_utilities.edfread(edfname)
		self.scannumbers.extend([scannumber]) # keep track of the scnanumbers
		return data, motors, counters, edfmats

	def loadscan(self,scannumbers,scantype='generic'):
		"""
		loads the files belonging to scan No "scannumber" and puts it into an instance
		of the container-class scan. the default scantype is 'generic', later the scans
		will be grouped (then added) based on the scantype
		"""
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers 
		for number in scannums:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number)
			# can assign some things here already (even if maybe redundant)
			monitor   = counters[self.monicolumn]
			monoangle = counters['pmonoa']
			energy    = counters[self.encolumn]
			# create an instance of "scan" class for every scan
			onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
			# assign one dictionary entry to each scan 
			self.scans[scanname] = onescan

	def loadloop(self,begnums,numofregions):
		"""
		loads a whole loop of scans based on their starting scannumbers and the number of single
		scans in a single loop
		begnums      = scannumbers of the first scans of each loop (is a list)
		numofregions = number of scans in each loop (integer)
		"""
		typenames = []
		for n in range(numofregions):
			typenames.append('edge'+str(n+1))
		numbers = []
		for n in range(len(begnums)):
			for m in range(numofregions):
				numbers.append((begnums[n]+m))
		typenames = typenames*len(begnums)
		for n in range(len(typenames)):
			thenumber = []
			thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
			# self = self.loadscan(thenumber,typenames[n])
			self.loadscan(thenumber,typenames[n])

	def loadelastic(self,scann):
		"""
		loads a scan using the loadscan function and sets the scantype attribute to 'elastic'
		"""
		self.loadscan(scann,'elastic')
	
	def loadlong(self,scann):
		"""
		loads a scan using the loadscan function and sets the scantype attribute to 'long'
		"""
		self.loadscan(scann,'long')

	def createrois(self,scannumbers):
		"""
		create rois object from this class and scannumbers
		"""
		rois_object = rois(self.scans,scannumbers)
		return rois_object

	def getautorois(self,scannumbers,kernel_size=5,colormap='jet'):
		"""
		define ROIs automatically using median filtering and a variable threshold
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getautorois(kernel_size,colormap)
		self.rois = rois_object.rois
		# feed an instance of the new roi_object class with these results for backward compatibility
		new_obj = xrs_rois.roi_object()
		new_obj.indices = rois_object.rois.inds
		self.roi_obj = new_obj

	def getlinrois(self,scannumbers,numrois=9,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getlinrois(numrois,logscaling,height,colormap)
		self.rois = rois_object.rois
		# feed an instance of the new roi_object class with these results for backward compatibility
		new_obj = xrs_rois.roi_object()
		new_obj.indices = rois_object.rois.inds
		self.roi_obj = new_obj

	def getzoomrois(self,scannumbers,numrois=9,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on a region of interest
		"""
		rois_object = self.createrois(scannumbers)
		rois_object.getzoomrois(numrois,logscaling,colormap)
		self.rois = rois_object.rois
		# feed an instance of the new roi_object class with these results for backward compatibility
		new_obj = xrs_rois.roi_object()
		new_obj.indices = rois_object.rois.inds
		self.roi_obj = new_obj

	def getlinedgerois(self,scannumbers,energyval,numrois=9,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points on a matrix, which is a difference of edf-files above and 
		below an energy value
		"""
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers
		# search for correct index in the first of the given scans
		scanname = 'Scan%03d' % scannums[0]
		index = self.scans[scanname].energy.flat[np.abs(self.scans[scanname].energy - energyval).argmin()]
		rois_object = self.createrois(scannumbers)
		rois_object.getlinedgerois(index,numrois,logscaling,height,colormap)
		self.rois = rois_object.rois
		# feed an instance of the new roi_object class with these results for backward compatibility
		new_obj = xrs_rois.roi_object()
		new_obj.indices = rois_object.rois.inds
		self.roi_obj = new_obj

	def getzoomedgerois(self,scannumbers,energyval,numrois=9,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on a matrix, which is a difference of edf-files above and 
		below an energy value
		"""
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers
		# search for correct index in the first of the given scans
		scanname = 'Scan%03d' % scannums[0]
		index = self.scans[scanname].energy.flat[np.abs(self.scans[scanname].energy - energyval).argmin()]
		rois_object = self.createrois(scannumbers)
		rois_object.getlinedgerois(index,numrois,logscaling,colormap)
		self.rois = rois_object.rois
		# feed an instance of the new roi_object class with these results for backward compatibility
		new_obj = xrs_rois.roi_object()
		new_obj.indices = rois_object.rois.inds
		self.roi_obj = new_obj

	def getrawdata(self):
		"""
		goes through all instances of the scan class and calls it's applyrois method
		to sum up over all rois
		"""
		if not self.rois:
			print 'please define some ROIs first.'
			pass
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ("integrating "+scan)
				self.scans[scan].applyrois(self.roi_obj.indices)

	def loadscandirect(self,scannumbers,rois,scantype='generic'):
		"""
		loads a scan without saving the edf files in matrices 
		needs ROIs as input (this needs testing)
		"""
		for number in scannumbers:
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number)
			# can assign some things here already (even if maybe redundant)
			monitor   = counters[self.monicolumn]
			monoangle = counters['pmonoa']
			energy    = counters[self.encolumn]
			# create an instance of "scan" class for every scan
			onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
			onescan.applyrois(self.roi_obj.indices)
			onescan.edfmats = [] # delete the edfmats
			self.scans[scanname] = onescan
		# if there are no rois yet, set the rois	
		if not self.rois:
			self.rois = rois

	def loadloopdirect(self,begnums,numofregions,rois):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'edge1', 'edge2', ...
		"""
		typenames = []
		for n in range(numofregions):
			typenames.append('edge'+str(n+1))
		numbers = []
		for n in range(len(begnums)):
			for m in range(numofregions):
				numbers.append((begnums[n]+m))
		typenames = typenames*len(begnums)
		for n in range(len(typenames)):
			thenumber = []
			thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
			self.loadscandirect(thenumber,typenames[n])

	def loadelasticdirect(self,scann,rois):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'elastic'
		"""
		self.loadscandirect(scann,rois,'elastic')
	
	def loadlongdirect(self,scann,rois):
		"""
		loads a scan without saving the edf files in matrices and sets
		the scantype attribute to 'long'
		"""
		self.loadscandirect(scann,rois,'long')

	def getspectrum(self):
		"""
		groups the instances of the scan class by their scantype attribute, 
		adds equal scans (each group of equal scans) and appends them 
		"""
		# find the groups 
		allgroups = findgroups(self.scans)
		for group in allgroups:
			# self.makegroup(group)
			onegroup = makegroup(group)
			self.groups[onegroup.gettype()] = onegroup
		self.energy,self.signals,self.errors = appendscans(self.groups)

	def geteloss(self):
		"""
		finds the center of mass of each roi and interpolates the signals onto a common grid
		"""	
		if not 'elastic' in self.groups:
			print 'please load/integrate at least one elastic scan first!'
		else:
			for n in range(len(self.rois)):
				self.cenom.append(trapz(self.groups['elastic'].signals[:,n]*self.groups['elastic'].energy,x=self.groups['elastic'].energy)/trapz(self.groups['elastic'].signals[:,n],x=self.groups['elastic'].energy))
			self.E0 = np.mean(self.cenom)
			self.eloss = (self.energy - self.cenom[0])*1e3 # energy loss in eV
			# interpolate all signals onto the energy scale of the first analyzer
			#if self.signals.any():
			for n in range(len(self.rois)):
				f = interp1d((self.energy-self.cenom[n])*1e3,self.signals[:,n], bounds_error=False,fill_value=0.0)
				self.signals[:,n] = f(self.eloss)

	def gettths(self,avtth):
		"""
		uses the defined TT_OFFSETS of the read_id20 class to set all scattering angles tth
		from the mean angle avtth of the analyzer module
		"""
		self.tth = np.absolute(avtth + self.TTH_OFFSETS)

	def getqvals(self,invangstr=False):
		"""
		calculates q-values from E0 and tth values in either atomic units (defalt) or
		inverse angstroms 
		"""
		theqvals = np.zeros_like(self.signals)
		if invangstr:
			for n in range(len(self.tth)):
				theqvals[:,n] = momtrans_inva(self.E0+self.eloss/1e3,self.E0,self.tth[n]) 
		else:
			for n in range(len(self.tth)):
				theqvals[:,n] = momtrans_au(self.E0+self.eloss/1e3,self.E0,self.tth[n])
		self.qvalues = theqvals

	def getqvals_energy(self,energy):
		"""
		returns all q-values at a certain energy loss
		"""
		ind = np.abs(self.eloss - energy).argmin()
		return self.qvalues[ind,:]

	def copy_edf_files(self,scannumbers,destdir):
		"""
		copies all edf-files from scan with scannumber or scannumbers into directory 'destdir'
		"""
		import shutil
		numbers = []
		if not isinstance(scannumbers,list):
			numbers.append(scannumbers)
		else:
			numbers = scannumbers
		fn = self.path + self.filename
		for n in range(len(numbers)):
			data, motors, counters = xrs_utilities.specread(fn,numbers[n])
			for m in range(len(counters['ccdno'])):
				ccdnumber = counters['ccdno'][m]+1
				edfname   = self.path + self.EDF_PREFIX + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
				shutil.copy2(edfname, destdir)



# horizontal zoom selector
"""
ax = subplot(111)
ax.plot(x,y)

def onselect(vmin, vmax):
    print vmin, vmax
span = SpanSelector(ax, onselect, 'horizontal')
"""

# example for rectangle selector from widgets
"""
from matplotlib.widgets import RectangleSelector

fig, ax = plt.subplots()
x = np.random.normal(size=1000)
y = np.random.normal(size=1000)
c = np.zeros((1000, 3))
c[:, 2] = 1  # set to blue
points = ax.scatter(x, y, s=20, c=c)

def selector_callback(eclick, erelease):
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    global c
    c[(x >= min(x1, x2)) & (x <= max(x1, x2))
      & (y >= min(y1, y2)) & (y <= max(y1, y2))] = [1, 0, 0]
    points.set_facecolors(c)
    fig.canvas.draw()
    
    
selector = RectangleSelector(ax, selector_callback,
                             drawtype='box', useblit=True,
                             button=[1,3], # don't use middle button
                             minspanx=5, minspany=5,
                             spancoords='pixels')
"""
# another one
"""
from matplotlib.widgets import  RectangleSelector
from pylab import *

def onselect(eclick, erelease):
  'eclick and erelease are matplotlib events at press and release'
  print ' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata)
  print ' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata)
  print ' used button   : ', eclick.button

def toggle_selector(event):
    print ' Key pressed.'
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print ' RectangleSelector deactivated.'
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print ' RectangleSelector activated.'
        toggle_selector.RS.set_active(True)

x = arange(100)/(99.0)
y = sin(x)
fig = figure
ax = subplot(111)
ax.plot(x,y)

toggle_selector.RS = RectangleSelector(ax, onselect, drawtype='line')
connect('key_press_event', toggle_selector)
show()

"""

class read_p01:
	"""
	a class for loading data fio, DSC, and binary files from PetraIII beamline P01 (experimental)
	"""
	def __init__(self,prefix):
		self.scans         = {}
		self.prefix        = prefix
		self.scannumbers   = []
		self.DET_PIXEL_NUM = 256
		self.TTH_OFFSETS   = np.array([])
		self.eloss    = []	# common eloss scale for all analyzers
		self.energy   = []	# common energy scale for all analyzers
		self.signals  = []	# signals for all analyzers
		self.errors   = []	# poisson errors
		self.qvalues  = []	# for all analyzers
		self.groups   = {}  # groups (such as 2 'elastic', or 5 'edge1', etc.)
		self.cenom    = []
		self.E0       = []
		self.tth      = []
		# input
		self.rois     = []		# can be set once the rois are defined to integrate scans
		# intermediates
		self.rawsignals = []	# structure of rawdata (energies, raw-signals, monitor, etc.)
		self.rawenergy  = []
		self.rawmonitor = []
	
	def loadscan(self,scannumbers,reps=None,usemonoa=False,scantype='generic'):
		monicol  = 6
		encol    = 0
		anglecol = 0
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers
		if not reps:
			for number in scannums:
				scanname = 'Scan%03d' % number
				fiodata, mats, mats1, mats2 = readp01scan(self.prefix,number)
				# can assign some things here already (even if maybe redundant)
				print np.shape(fiodata)				
				monitor   = fiodata[:,monicol]
				monoangle = fiodata[:,anglecol]				
				if usemonoa:
					# need to flipup energy scale since monoangle is scanned from small to large angles
					inds    = np.argsort(energy_monoangle(fiodata[:,anglecol]))
					fiodata = fiodata[inds,:]						
					energy  = energy_monoangle(fiodata[:,anglecol])
					mats    = mats[inds,:,:]
					mats1   = mats1[inds,:,:]
					mats2   = mats2[inds,:,:]
					monitor = monitor[inds]
					monoangle= monoangle[inds]
				else:
					energy = fiodata[:,encol]/1.0e3
				
				# create an instance of "scan" class for every scan
				onescan = scan(mats,number,energy,monoangle,monitor,[],[],fiodata,scantype)
				# assign one dictionary entry to each scan 
				self.scans[scanname] = onescan
		else:
			for rep in reps:
				scanname = 'Scan%03d' % scannums[0] + 'r%01d' % rep
				fiodata, mats, mats1, mats2 = readp01scan_rep(self.prefix,scannums[0],rep)
				# can assign some things here already (even if maybe redundant)
				monitor   = fiodata[:,monicol]
				monoangle = fiodata[:,anglecol]				
				if usemonoa:
					# need to flipup energy scale since monoangle is scanned from small to large angles
					inds    = np.argsort(energy_monoangle(fiodata[:,anglecol]))
					fiodata = fiodata[inds,:]						
					energy  = energy_monoangle(fiodata[:,anglecol])
					mats    = mats[inds,:,:]
					mats1   = mats1[inds,:,:]
					mats2   = mats2[inds,:,:]
					monitor = monitor[inds]
					monoangle= monoangle[inds]
				else:
					energy = fiodata[:,encol]/1.0e3

				# create an instance of "scan" class for every scan
				onescan = scan(mats,scannums[0],energy,monoangle,monitor,[],[],fiodata,scantype)
				# assign one dictionary entry to each scan 
				self.scans[scanname] = onescan

	def loadelastic(self,scann,reps=None):
		self = self.loadscan(scann,reps,scantype='elastic')
	
	def loadedge(self,scann,reps=None,usemonoa=True):
		self = self.loadscan(scann,reps,usemonoa,scantype='edge')
	
	def loadlong(self,scann,reps=None):
		self = self.loadscan(scann,reps,scantype='long')

	def createrois(self,scannumbers):
		"""
		create rois object from this class and scannumbers
		"""
		rois_object = rois(self.scans,scannumbers)
		return rois_object		

	def getzoomrois(self,scannumbers,numrois=12,logscaling=True,colormap='jet'):
		rois_object = self.createrois(scannumbers)
		rois_object.getzoomrois(numrois,logscaling,colormap)
		self.rois = rois_object.rois

	def getautorois(self,scannumbers,kernel_size=5,colormap='jet'):
		rois_object = self.createrois(scannumbers)
		rois_object.getautorois(kernel_size,colormap)
		self.rois = rois_object.rois

	def getlinrois(self,scannumbers,numrois=9,logscaling=True,height=5,colormap='jet'):
		rois_object = self.createrois(scannumbers)
		rois_object.getlinrois(numrois,logscaling,height,colormap)
		self.rois = rois_object.rois

	def getrawdata(self):
		if not self.rois:
			print 'please define some ROIs first.'
			pass
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ("integrating "+scan)
				self.scans[scan].applyrois(self.rois)

	def getspectrum(self):
		# find the groups 
		allgroups = findgroups(self.scans)
		for group in allgroups:
			# self.makegroup(group)
			onegroup = makegroup(group)
			self.groups[onegroup.gettype()] = onegroup
		self.energy,self.signals,self.errors = appendscans(self.groups)

	def plotoneq(self,scannumbers,whichq,norm=True,reps=False):
		if not rep:
			for number in scannumbers:
				scanname = 'Scan%03d' % number
				if norm:
					plot(self.scans[scanname].energy,self.scans[scanname].signals[:,whichq]/self.scans[scanname].monitor)
				else:
					plot(self.scans[scanname].energy,self.scans[scanname].signals[:,whichq])
			show()
		else:
			for number in scannumbers:
				for rep in reps:
					scanname = 'Scan%03d' % number + 'r%01d' % rep 
					if norm:
						plot(self.scans[scanname].energy,self.scans[scanname].signals[:,whichq]/self.scans[scanname].monitor)
					else:
						plot(self.scans[scanname].energy,self.scans[scanname].signals[:,whichq])
			show()

