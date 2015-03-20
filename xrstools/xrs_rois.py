#!/usr/bin/python
# Filename: xrs_rois.py

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
import matplotlib.pyplot as plt

from matplotlib.path import Path
from xrs_utilities import *
from math_functions import *
from matplotlib.widgets import Cursor, Button
from scipy.ndimage import measurements
from scipy import signal

class container:
	"""
	Random container class to hold values
	"""
	def __init__(self):
		pass




sl={}
sl[0  ] = slice(0  ,256)
sl[256] = slice(256,512)
sl[512] = slice(512,768)

V147 = [1,4,7,10,2,5,8,11,3,6,9,12]
V1296 = [12,9,6,3,11,8,5,2,10,7,4,1]
V1074 = [10,7,4,1,11,8,5,2,12,9,6,3     ]
V369  = [3,6,9,12,2,5,8,11,1,4,7,10             ]
geo_informations = {(256,768): { "DET_PIXEL_NUM":256, "geo":[256,768], "nofrois":36,
                                 "subnames": ["RD" ,"LU","B"],
                                 "subgeos" : [(sl[0] ,sl[0]),
                                              (sl[0],sl[256]),
                                              (sl[0],sl[512])]   ,
                                 "analyser_nIDs": {"LU":{"3x4":V1074,"Vertical": V147},  
                                                   "RD":{"3x4":V147,"Vertical": V147},
                                                   "B": {"3x4":V147,"Vertical": V147}  
                                                   }
                                 },
                    (512,768) : { "DET_PIXEL_NUM":256, "geo":[512,768],"nofrois":72,
                                  "subnames":["HR" ,"HL","HB","VD" , "VU","VB", ] ,
                                  "subgeos"  :[(sl[0],sl[0]     ),
                                               (sl[0],sl[256]   ),
                                               (sl[0],sl[512]   ),
                                               (sl[256],sl[0]   ),
                                               (sl[256],sl[256] ),
                                               (sl[256],sl[512] )],
                                  "analyser_nIDs": {"VU":{"3x4":V1074,"Vertical": V147},
                                                    "VD":{"3x4":V147 ,"Vertical": V147},
                                                    "VB":{"3x4":V147 ,"Vertical": V147},
                                                    "HR":{"3x4":V369 ,"Vertical":V1296},
                                                    "HL":{"3x4":V1296,"Vertical":V1296},
                                                    "HB":{"3x4":V1296,"Vertical":V1296},
                                                    }
                                  },
                    (256,256):{"DET_PIXEL_NUM":256, "geo":[256,256],"nofrois":1,"subnames":["DETECTOR"],"subgeos" : [(sl[0] ,sl[0])],
                               "analyser_nIDs": {"DETECTOR":{"3x4":V147},"Vertical": V147}
                       }
            }

def get_geo_informations(shape):
    return geo_informations[shape]



class roi_object:
	"""
	Container class to hold all relevant information about given ROIs.
	"""
	def __init__(self):
                self.roi_matrix     = np.array([]) # single matrix of zeros, ones, twos, ... , n's (where n is the number of ROIs defined)
		self.red_rois       = {}           # dictionary, one entry for each ROI, each ROI has an origin and a rectangular box of ones and zeros defining the ROI
		self.indices        = [] # list of list of tuples (one list of tuples for each ROI)
		self.number_of_rois = 0            # number of ROIs defined
		self.kind           = [] # keyword (e.g. 'zoom', 'line', 'auto', etc.), certain features (esp. in imaging) are only available for certain kinds of ROIs
		self.x_indices      = [] # list of numpy arrays of x-indices (for each ROI)
		self.y_indices      = [] # list of numpy arrays of y-indices (for each ROI)
		self.masks          = [] # 3D numpy array with slices of zeros and ones (same size as detector image) for each roi

		
	def load_rois_fromMasksDict(self, masksDict, newshape=None, kind="zoom"):
		self.kind=kind
		self.red_rois = masksDict
		if newshape is not None:
			self.roi_matrix = np.zeros(newshape)
		self.roi_matrix = convert_redmatrix_to_matrix( masksDict,self.roi_matrix , offsetX=0, offsetY=0)

		self.masks          = convert_roi_matrix_to_masks(self.roi_matrix)
		self.indices        = convert_matrix_rois_to_inds(self.roi_matrix)
		selfnumber_of_rois = np.amax(self.roi_matrix)
		self.x_indices      = convert_inds_to_xinds(self.indices)
		self.y_indices      = convert_inds_to_yinds(self.indices)


	def show_rois(self):
		"""
		Creates a figure with the defined ROIs as numbered boxes on it.
		"""
		roi_matrix = self.roi_matrix

		# check if there are ROIs defined
		if not np.any(roi_matrix):
			print 'Please select some rois first.'

		# make a figure
		else:
			plt.ioff()
			plt.imshow(roi_matrix)
			plt.xlabel('x-axis [pixel]')
			plt.ylabel('y-axis [pixel]')

			# add a label with the number to the center of each ROI
			for ii in range(int(np.amax(roi_matrix))):
				# find center of the ROI and write the label
				inds    = np.where(roi_matrix[:,:] == ii+1)
				xcenter = np.mean(inds[1])
				ycenter = np.mean(inds[0])
				string  = '%02d' % (ii+1)
				plt.text(xcenter,ycenter,string)
			plt.show()
			plt.ion()

	def get_number_of_rois(self):
		return self.number_of_rois

	def get_indices(self):
		return self.indices

	def get_x_indices(self):
		return self.x_indices

	def get_y_indices(self):
		return self.y_indices

	def get_bounding_boxes(self):
		return self.bounding_boxes

	def get_masks(self):
		return self.masks

class roi_finder:
	"""
	Class to define ROIs from a 2D image.
	"""

	def __init__(self):
		self.roi_obj = roi_object() # empty roi object

	def deleterois(self):
		"""
		Clear the existing ROIs by creating a fresh roi_object.
		"""
		self.roi_obj = roi_object()

	def get_linear_rois(self,input_image,logscaling=True,height=5,colormap='jet',interpolation='nearest'):
		"""
		Define ROIs by clicking two points on a 2D image.
		number_of_rois = integer defining how many ROIs should be determined
		input_object   = 2D array, scan_object, or dictionary of scans to define the ROIs from
		logscaling     = boolean, to determine wether the image is shown on a log-scale (default = True)
		height         = integer defining the height (in pixels) of the ROIs
		"""
		# make sure the matplotlib interactive mode is off
		plt.ioff()

		# clear all existing rois
		self.deleterois()

		# check that the input is a 2d matrix
		if not len(input_image.shape) == 2:
			print 'Please provide a 2D numpy array as input!'
			return

		# calculate the logarithm if 'logscaling' == True
		if logscaling:
			# set all zeros to ones:
			input_image[input_image[:,:] == 0.0] = 1.0
			input_image = np.log(np.abs(input_image))

		# prepare a figure
		fig, ax = plt.subplots()
		plt.subplots_adjust(bottom=0.2)
		cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

		# Initialize suptitle, which will be updated
		titlestring = 'Start by clicking the \'Next\'-button.'
		titleInst=plt.suptitle(titlestring)

		# generate an image to be displayed
		figure_obj = plt.imshow(input_image,interpolation=interpolation)

		# set the colormap for the image
		figure_obj.set_cmap(colormap)

		rois   = []
		class Index:
			ind  = 0
			def next(self, event):
				self.ind += 1
				titlestring = 'Click two points for ROI Nr. %02d, hit \'Next\' to continue, \'Finish\' to end.' % self.ind
				titleInst.set_text(titlestring) # Update title
				one_roi = define_lin_roi(height,input_image.shape)
				for index in one_roi:
					input_image[index[0],index[1]] += 1.0e6
				figure_obj.set_data(input_image)
				plt.hold(True)
				rois.append(one_roi)
			def prev(self, event):
				self.ind -= 1
				try:
					titlestring = 'Click the \'Next\' button again to continue.'
					titleInst.set_text(titlestring) # Update title
					# print titlestring
					for index in rois[-1]:
						input_image[index[0],index[1]] -= 1.0e6
					figure_obj.set_data(input_image)
					plt.hold(True)
					rois.pop()
				except:
					pass
			def close(self, event):
				plt.close()
			def dmy(self, event):
				pass # adding a dummy function for the dummy button

		callback = Index()
		axprev   = plt.axes([0.5, 0.05, 0.1, 0.075])
		axnext   = plt.axes([0.61, 0.05, 0.1, 0.075])
		axclose  = plt.axes([0.72, 0.05, 0.1, 0.075])
		axdmy    = plt.axes([0.001, 0.001, 0.001, 0.001]) # for some reason the first botton disappears when clicked
		bdmy     = Button(axdmy,'')                        # which is why I am including a dummy button here
		bdmy.on_clicked(callback.dmy)                   # this way, the real buttons work
		bnext    = Button(axnext, 'Next')
		bnext.on_clicked(callback.next)
		bprev    = Button(axprev, 'Back')
		bprev.on_clicked(callback.prev)
		bclose   = Button(axclose, 'Finish')
		bclose.on_clicked(callback.close)
		plt.show()

		# assign the defined rois to the roi_object class
		self.roi_obj.roi_matrix     = convert_inds_to_matrix(rois,input_image.shape)
		self.roi_obj.red_rois       = convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = rois 
		self.roi_obj.kind           = 'linear'
		self.roi_obj.x_indices      = convert_inds_to_xinds(rois)
		self.roi_obj.y_indices      = convert_inds_to_yinds(rois)
		self.roi_obj.masks          = convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def get_zoom_rois(self,input_image,logscaling=True,colormap='jet',interpolation='nearest'):
		"""
		Define ROIs by clicking two points on a 2D image.
		number_of_rois = integer defining how many ROIs should be determined
		input_object   = 2D array, scan_object, or dictionary of scans to define the ROIs from
		logscaling     = boolean, to determine wether the image is shown on a log-scale (default = True)
		height         = integer defining the height (in pixels) of the ROIs
		"""
		# make sure the matplotlib interactive mode is off
		plt.ioff()

		# clear all existing rois
		self.deleterois()

		# check that the input is a 2d matrix
		if not len(input_image.shape) == 2:
			print 'please provide a 2D numpy array as input!'
			return                

		# calculate the logarithm if 'logscaling' == True
		if logscaling:
			# set all zeros to ones:
			input_image[input_image[:,:] == 0.0] = 1.0
			input_image = np.log(np.abs(input_image))

		# prepare a figure
		fig, ax = plt.subplots()
		plt.subplots_adjust(bottom=0.2)
		cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

		# Initialize suptitle, which will be updated
		titlestring = 'Start by clicking the \'Next\' button.'
		titleInst=plt.suptitle(titlestring)

		# generate an image to be displayed
		figure_obj = plt.imshow(input_image,interpolation=interpolation)

		# activate the zoom function already
		thismanager = plt.get_current_fig_manager()
		thismanager.toolbar.zoom()

		# set the colormap for the image
		figure_obj.set_cmap(colormap)
		#plt.colorbar()

		# initialize a matrix for the rois (will be filled with areas of ones, twos, etc
		rois = []

		# print info to start: 
		print 'Start by clicking the \'Next\' button.'

		class Index:
			ind          = 0
			initstage    = True
			next_clicked = False
			def next(self, event):
				# for some reason, the first time this is used, it doesn't work, so here is one dummy round
				if self.initstage:
					#self.ind += 1
					self.initstage = False
					plt.sca(ax)
					one_roi = define_zoom_roi(input_image,verbose=True)
					#for index in one_roi:
					#	input_image[index[0],index[1]] *= 1.0
					# reset the matrix to be displayed
					figure_obj.set_data(input_image)
					# reset the zoom
					plt.xlim(0.0,input_image.shape[1])
					plt.ylim(input_image.shape[0],0.0)
					plt.draw()
					titlestring = 'Zoom in to define ROI Nr. %02d, hit \'Next\' to continue.' % (self.ind + 1)
					titleInst.set_text(titlestring) # Update title
				else:
					self.ind += 1
					plt.sca(ax)
					one_roi = define_zoom_roi(input_image)
					for index in one_roi:
						input_image[index[0],index[1]] += 1.0e10
					# reset the matrix to be displayed
					figure_obj.set_data(input_image)
					# reset the zoom
					plt.xlim(0.0,input_image.shape[1])
					plt.ylim(input_image.shape[0],0.0)
					plt.draw()
					rois.append(one_roi)
					titlestring = 'Zoom in to define ROI Nr. %02d, hit \'Next\' to continue, \'Finish\' to end.' % (self.ind + 1)
					titleInst.set_text(titlestring) # Update title
                                                
			def prev(self, event):
				self.ind -= 1
				titlestring = 'Undoing ROI Nr. %02d. Zoom again, click the \'Next\' button to continue.' % (self.ind + 1)
				titleInst.set_text(titlestring) # Update title
				#thedata[roimatrix == self.ind+1]   -= 1.0e6
				#roi_matrix[roimatrix == self.ind+1] = 0.0
				for index in rois[-1]:
					input_image[index[0],index[1]] -= 1.0e10
				figure_obj.set_data(input_image)
				plt.hold(True)
				rois.pop()

			def close(self, event):
				plt.sca(ax)
				one_roi = define_zoom_roi(input_image)
				for index in one_roi:
					input_image[index[0],index[1]] += 1.0e10
				# reset the matrix to be displayed
				figure_obj.set_data(input_image)
				# reset the zoom
				plt.xlim(0.0,input_image.shape[1])
				plt.ylim(input_image.shape[0],0.0)
				plt.draw()
				rois.append(one_roi)
				titlestring = 'Last ROI is Nr. %02d.' % (self.ind + 1)
				titleInst.set_text(titlestring) # Update title
				plt.close()

			def dmy(self, event):
				pass # adding a dummy function for the dummy button

		callback = Index()
		axprev   = plt.axes([0.5, 0.05, 0.1, 0.075])
		axnext   = plt.axes([0.61, 0.05, 0.1, 0.075])
		axclose  = plt.axes([0.72, 0.05, 0.1, 0.075])
		axdmy    = plt.axes([0.001, 0.001, 0.001, 0.001]) # for some reason the first botton disappears when clicked
		bdmy     = Button(axdmy,'')                       # which is why I am including a dummy button here
		bdmy.on_clicked(callback.dmy)                     # this way, the real buttons work
		bnext    = Button(axnext, 'Next')
		bnext.on_clicked(callback.next)
		bprev    = Button(axprev, 'Back')
		bprev.on_clicked(callback.prev)
		bclose   = Button(axclose, 'Finish')
		bclose.on_clicked(callback.close)
		plt.show()

		# assign the defined rois to the roi_object class
		self.roi_obj.roi_matrix     = (convert_inds_to_matrix(rois,input_image.shape))
		self.roi_obj.red_rois       = convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = rois 
		self.roi_obj.kind           = 'zoom'
		self.roi_obj.x_indices      = convert_inds_to_xinds(rois)
		self.roi_obj.y_indices      = convert_inds_to_yinds(rois)
		self.roi_obj.masks          = convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def get_auto_rois(self,input_image,kernel_size=5,threshold=100.0,logscaling=True,colormap='jet',interpolation='bilinear'):
		"""
		Define ROIs by choosing a threshold using a slider bar under the figure. In this function, the entire 
		detector is shown.
		input_image = 2D numpy array with the image to be displayed
		kernal_size = integer defining the median filter window (has to be odd)
		theshold    = initial number defining the upper end value for the slider bar (amax(input_image)/threshold defines this number), can be within GUI
		logscaling  = boolean, if True (default) the logarithm of input_image is displayed
		colormap    = matplotlib color scheme used in the display
		interpolation = matplotlib interpolation scheme used for the display
		"""
		# make sure the matplotlib interactive mode is off
		plt.ioff()

		# clear all existing rois
		self.deleterois()

		# clear existing figure
		# plt.clf()

		# calculate the logarithm if 'logscaling' == True
		if logscaling:
			# set all zeros to ones:
			input_image[input_image[:,:] == 0.0] = 1.0
			input_image = np.log(np.abs(input_image))


		ax = plt.subplot(111) 
		plt.subplots_adjust(left=0.05, bottom=0.2)

		# print out some instructions
		plt.suptitle('Use the slider bar to select ROIs, close the plotting window when satisfied.')

		# initial threshold value
		thres0 = 0.0

		# create a figure object
		figure_obj = plt.imshow(input_image,interpolation=interpolation)
		figure_obj.set_cmap(colormap)

		# prepare the slider bar
		thresxcolor = 'lightgoldenrodyellow'
		thresxamp  = plt.axes([0.2, 0.10, 0.55, 0.03], axisbg=thresxcolor)
		maxthreshold=np.floor(np.amax(input_image)) # maximum of slider
		sthres = plt.Slider(thresxamp, 'Threshold', 0.0, maxthreshold, valinit=thres0)

		textBox=plt.figtext(0.50, 0.065, 'Multiplier: 1.0',verticalalignment='center')

		# define what happens when the slider is touched
		def update(val):
			# parse a threshold from the slider
			thres     = sthres.val*thresMultiplier.factor
			# median filter the image
			newmatrix = signal.medfilt2d(input_image, kernel_size=kernel_size)
			# set pixels below the threshold to zero
			belowthres_indices = newmatrix < thres
			newmatrix[belowthres_indices] = 0
			# identify connected regions (this is already the roi_matrix)
			self.roi_obj.roi_matrix,numfoundrois = measurements.label(newmatrix)
			print str(numfoundrois) + ' ROIs found!'	
			figure_obj.set_data(newmatrix)
			plt.draw()

		# Buttons for changing multiplier for the value of slider
		class thresMultiplierClass:
			factor = 1.0;
			def __new__(cls):
				return self.factor
			def increase(self,event):
				self.factor *=2.0
				textBox.set_text('Multiplier: ' + str(self.factor))
				return self.factor
			def decrease(self,event):
				self.factor /=2.0
				textBox.set_text('Multiplier: ' + str(self.factor))
				return self.factor

		# call the update function when the slider is touched
		sthres.on_changed(update)

		thresMultiplier = thresMultiplierClass()
		axincrease   = plt.axes([0.8, 0.05, 0.05, 0.03])
		axdecrease   = plt.axes([0.7, 0.05, 0.05, 0.03])
		bnincrease    = Button(axincrease, 'x 2')
		bndecrease    = Button(axdecrease, '/ 2')
		bnincrease.on_clicked(thresMultiplier.increase)  # First change threshold
		bnincrease.on_clicked(update)		         # Then update image
		bndecrease.on_clicked(thresMultiplier.decrease)
		bndecrease.on_clicked(update)
		# ADDITION ENDS

		plt.show()

		# assign the defined rois to the roi_object class
		self.roi_obj.red_rois       = convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = convert_matrix_rois_to_inds(self.roi_obj.roi_matrix)
		self.roi_obj.kind           = 'auto'
		self.roi_obj.x_indices      = convert_inds_to_xinds(self.roi_obj.indices)
		self.roi_obj.y_indices      = convert_inds_to_yinds(self.roi_obj.indices)
		self.roi_obj.masks          = convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def get_polygon_rois(self,input_image,modind=-1,logscaling=True,colormap='jet',interpolation='nearest'):
		"""
		Define ROIs by clicking arbitrary number of points on a 2D image:
		LEFT CLICK to define the corner points of polygon, 
		MIDDLE CLICK to finish current ROI and move to the next ROI,
		RIGHT CLICK to cancel the previous point of polygon
		input_object   = 2D array, scan_object, or dictionary of scans to define the ROIs from
		modind	     = integer to identify module, if -1 (default), no module info will be in title (the case of one big image)
		logscaling     = boolean, to determine wether the image is shown on a log-scale (default = True)

		"""
		# make sure the matplotlib interactive mode is off
		plt.ioff()

		# clear all existing rois
		self.deleterois()

		# check that the input is a 2d matrix
		if not len(input_image.shape) == 2:
			print 'Please provide a 2D numpy array as input!'
			return

		# calculate the logarithm if 'logscaling' == True
		if logscaling:
			# set all zeros to ones:
			input_image[input_image[:,:] == 0.0] = 1.0
			input_image = np.log(np.abs(input_image))

		# prepare a figure
		fig, ax = plt.subplots()
		plt.subplots_adjust(bottom=0.2)
		
		moduleNames='VD:','HR:','VU:','HL:','VB:','HB:',''  # for plot title

		# Initialize suptitle, which will be updated
		titlestring = ''
		titleInst=plt.suptitle(titlestring)

		cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

		# generate an image to be displayed
		figure_obj = plt.imshow(input_image,interpolation=interpolation)

		# set the colormap for the image
		figure_obj.set_cmap(colormap)

		rois   = []
		class Index:
			ind  = 1
			def next(self, event):

				titlestring = '%s next ROI is Nr. %02d:\n Left button to new points, middle to finish ROI. Hit \'Finish\' to end with this image.' % (moduleNames[modind], self.ind)
				titleInst.set_text(titlestring) # Update title

				# Try needed, as FINISH button closes the figure and ginput() generates _tkinter.TclError 
				try:
					one_roi = define_polygon_roi(input_image.shape)
					for index in one_roi:
						input_image[index[0],index[1]] += 1.0e6
					figure_obj.set_data(input_image)
					plt.hold(True)
					plt.draw()
					rois.append(one_roi)
					self.ind += 1
					# Begin defining the next ROI right after current
					self.next(self)
				except KeyboardInterrupt:	# to prevent "dead" figures
					plt.close()
					pass		
				except:
					pass
			def prev(self, event):
				self.ind -= 1
				for index in rois[-1]:
					input_image[index[0],index[1]] -= 1.0e6
				figure_obj.set_data(input_image)
				plt.hold(True)
				plt.draw()
				rois.pop()
				self.next(self)

			def close(self, event):
				plt.close()
			def dmy(self, event):
				pass # adding a dummy function for the dummy button

		callback = Index()	
		axprev   = plt.axes([0.5, 0.05, 0.1, 0.075])
		axnext   = plt.axes([0.61, 0.05, 0.1, 0.075])
		axclose  = plt.axes([0.72, 0.05, 0.1, 0.075])
		axdmy    = plt.axes([0.001, 0.001, 0.001, 0.001]) # for some reason the first botton disappears when clicked
		bdmy     = Button(axdmy,'')                        # which is why I am including a dummy button here
		bdmy.on_clicked(callback.dmy)                   # this way, the real buttons work
		bnext    = Button(axnext, 'Next')
		bnext.on_clicked(callback.next)
		bprev    = Button(axprev, 'Back')
		bprev.on_clicked(callback.prev)
		bclose   = Button(axclose, 'Finish')
		bclose.on_clicked(callback.close)

		# START: initiate NEXT button press
		callback.next(self)		

		# assign the defined rois to the roi_object class
		self.roi_obj.roi_matrix     = convert_inds_to_matrix(rois,input_image.shape)
		self.roi_obj.red_rois       = convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		self.roi_obj.indices        = rois 
		self.roi_obj.kind           = 'polygon'
		self.roi_obj.x_indices      = convert_inds_to_xinds(rois)
		self.roi_obj.y_indices      = convert_inds_to_yinds(rois)
		self.roi_obj.masks          = convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
		self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

	def show_rois(self):
		"""
		Creates a figure with the defined ROIs as numbered boxes on it.
		"""
		roi_matrix = self.roi_obj.roi_matrix

		# check if there are ROIs defined
		if not np.any(roi_matrix):
			print 'Please select some rois first.'

		# make a figure
		else:
			plt.imshow(roi_matrix)
			plt.xlabel('x-axis [pixel]')
			plt.ylabel('y-axis [pixel]')

			# add a label with the number to the center of each ROI
			for ii in range(int(np.amax(roi_matrix))):
				# find center of the ROI and write the label
				inds    = np.where(roi_matrix[:,:] == ii+1)
				xcenter = np.mean(inds[1])
				ycenter = np.mean(inds[0])
				string  = '%02d' % (ii+1)
				plt.text(xcenter,ycenter,string)
				plt.show()

def define_lin_roi(height,image_shape,verbose=False):
	"""
	Lets you pick 2 points on a current image and returns a linear ROI of
	height (2*height+1).
	height      = number of pixels that define the height of the ROI
	image_shape = tuple with shape of the current image (i.e. (256,256))
	"""
	endpoints = list(np.round(plt.ginput(2,timeout=-1)))

	# check that selected points are in the image
	for point in endpoints:
		if point[0] < 0.0:
			point = 0
		if point[0] > image_shape[1]:
			point = image_shape[1]
		if point[1] < 0.0:
			point[1] = 0
		if point[1] > image_shape[0]:
			point[1] = image_shape[0]

	# check that point 2 is bigger than point 1 in the x direction
	if endpoints[1][0]< endpoints[0][0]:
		endpoints[0],endpoints[1] = endpoints[1],endpoints[0]

	# print the limits of the rectangle in the shell
	if not verbose:
		print 'The selected points are: ', [[endpoints[0][1],endpoints[0][0]],[endpoints[1][1],endpoints[1][0]]]

	roix = np.arange(endpoints[0][0],endpoints[1][0])
	roiy = [round(num) for num in np.polyval(np.polyfit([endpoints[0][0],endpoints[1][0]],[endpoints[0][1],endpoints[1][1]],1),roix)]
	roi  = []
	height = np.arange(-height,height)
	for n in range(len(roix)):
		for ind in height:
			roi.append((roiy[n]+ind,roix[n]))
	return roi

def define_zoom_roi(input_image,verbose=False):
	"""
	Parses the current figure limits and uses them to define a rectangle of "
	roi_number"s in the matrix given by roi_matrix.
	input_image  = unzoomed figure
	roi_matrix   = current roi matrix which will be altered 
	"""
	# input_image.shape prints  (y-length, x-length)
	# parse the figure limits from the current zoom
	limits = np.round(plt.axis()) # x-min, x-max, y-max, y-min
	# frac_limits = plt.axis() # x-min, x-max, y-max, y-min as floats
	# limits = [np.round(frac_limits[0]), np.ceil(frac_limits[1]), np.ceil(frac_limits[2]), np.floor(frac_limits[3])] # x-min, x-max, y-max, y-min
	# check that selected zoom area is not outside the image
	inds = limits < 0
	limits[inds] = 0
	if limits[1] > input_image.shape[1]:
		limits[1] = input_image.shape[1]
	if limits[2] > input_image.shape[0]:
		limits[2] = input_image.shape[0]

	# sort the limits in ascenging order
	limitsy = limits[2:4] # vertical
	limitsy.sort()
	limitsx = limits[0:2] # horizontal
	limitsx.sort()

	# print the limits of the rectangle in the shell
	if not verbose:
		print 'The selected limits are: ', limitsx, limitsy

	# prepare a n*m matrix with one roi
	T = np.zeros(input_image.shape)
	T[limitsy[0]:limitsy[1],limitsx[0]:limitsx[1]] = 1
	indsy,indsx = np.where(T == 1)

	roi  = []
	for n in range(len(indsx)):
		roi.append((indsy[n],indsx[n]))
	return roi


def convert_redmatrix_to_matrix( masksDict,mask, offsetX=0, offsetY=0):
    for key, (pos,M)  in masksDict.iteritems():
        num=int("".join([c for c in key if c.isdigit()]))
        S = M.shape
        inset =    (slice(offsetY+pos[0]  , offsetY+pos[0]+S[0]   ), slice(  offsetX+pos[1]  , offsetX+pos[1]+S[1] ) )
	print inset
	print num

        mask[  inset   ] =  num+1
    return mask



def convert_inds_to_matrix(ind_rois,image_shape):
    """
    Converts a ROI defined by a list of lists of tuples into a ROI
    that is defined by an array containing zeros, ones, twos, ..., n's, 
    where n is the number of ROIs.
    ind_rois    = list of lists with pairs of pixel indices
    image_shape = touple defining the shape of the matrix for which the ROIs are valid 
    """
    roi_matrix = np.zeros(image_shape)
    counter = 1
    for pixel in ind_rois:
        for xyind in pixel:
            roi_matrix[xyind[0],xyind[1]] = counter
        counter += 1
    return roi_matrix




def convert_matrix_to_redmatrix(matrix_rois, labelformat= 'ROI%02d'):
	"""
	Converts a ROI defined by an array containing zeros, ones, twos, ..., n's, 
	where n is the number of ROIs, into a dictionary with keys 'ROI00',
	'ROI01', ..., 'ROInn'. Each entry of the dictionary is a list containing a
	tuple with the origin and the reduced ROI.
	matrix_roi = numpy array
	"""
	redmatrix = {}
	for ind in range(int(np.amax(matrix_rois))):
		indices   = np.where(matrix_rois == ind+1)
		origin    = (np.amin(indices[0]), np.amin(indices[1]))
		bound_box = matrix_rois[np.amin(indices[0]):np.amax(indices[0]),np.amin(indices[1]):np.amax(indices[1])]
		thekey    =labelformat % ind
		redmatrix[thekey] = [origin,bound_box]
	return redmatrix

def convert_inds_to_xinds(roi_inds):
    """
    Converts ROIs defined in lists of lists of x-y-coordinate tuples into
    a list of x-coordinates only.
    """
    xind_rois = []
    for roi in roi_inds:
        xinds = []
        for pair in roi:
            xinds.append(pair[0])
        xind_rois.append(xinds)
    return xind_rois

def convert_inds_to_yinds(roi_inds):
    """
    Converts ROIs defined in lists of lists of x-y-coordinate tuples into
    a list of y-coordinates only.
    """
    yind_rois = []
    for roi in roi_inds:
        yinds = []
        for pair in roi:
            yinds.append(pair[1])
        yind_rois.append(yinds)
    return yind_rois

def convert_roi_matrix_to_masks(roi_matrix):
    """
    Converts a 2D ROI matrix with zeros, ones, twos, ..., n's (where n is the number of ROIs) to
    a 3D matrix with one slice of zeros and ones per ROI.
    """
    # get the shape
    roi_masks = np.zeros((np.amax(roi_matrix),roi_matrix.shape[0],roi_matrix.shape[1]))
    for ii in range(int(np.amax(roi_matrix))):
        inds = np.where(roi_matrix[:,:] == ii+1)
	for jj in range(len(inds[0])):
		roi_masks[ii,inds[0][jj],inds[1][jj]] = ii+1
    return roi_masks

def convert_matrix_rois_to_inds(roi_matrix):
	"""
	Converts a 2D ROI matrix with zeros, ones, twos, ..., n's (where n is the number of ROIs) to
	a list of lists each of which has tuples with coordinates for each pixel in each roi.
	"""
	rois = []
	number_of_rois = np.amax(roi_matrix)
	for ii in range(int(number_of_rois)):
		inds = np.where(roi_matrix[:,:] == ii+1)
		oneroi = []
		for i in range(len(inds[0])):
			oneroi.append( (inds[0][i],inds[1][i]) )        
		rois.append(oneroi)
	return rois

def test_roifinder(roi_type_str, imagesize = [512,768], scan = None ):
    """
    Runs the roi_finder class on a random image of given type for 
    testing purposes.
    scan[0] = absolute path to a spec file
    scan[1] = energycolumn='energy'
    scan[2] = monitorcolumn='kap4dio'
    scan[3] = scan number from which to take images
    """
    strings = ['zoom','linear','auto']
    if not roi_type_str in strings:
	    print 'Only ' + str(strings) + ' testable, choose one of them!' 
	    return

    # make a random image 
    if not scan:
	    rand_image = np.random.rand(imagesize[0],imagesize[1])
    else: 
	    import xrs_read, xrs_utilities
	    read_obj = xrs_read.read_id20(scan[0],energycolumn=scan[1],monitorcolumn=scan[2])
	    read_obj.loadelastic(scan[3])
	    key = 'Scan%03d' % scan[3]
	    rand_image = xrs_utilities.sumx(read_obj.scans[key].edfmats)

    # create a roi_finder object 
    roi_finder_obj = roi_finder()

    if roi_type_str == 'zoom':
	    roi_finder_obj.get_zoom_rois(rand_image,logscaling=True,colormap='jet',interpolation='nearest')
	    roi_finder_obj.show_rois()

    elif roi_type_str == 'linear':
	    roi_finder_obj.get_linear_rois(rand_image,logscaling=True,height=5,colormap='jet',interpolation='nearest')
	    roi_finder_obj.show_rois()
	    
    elif roi_type_str == 'auto':
	    roi_finder_obj.get_auto_rois(rand_image,kernel_size=5,threshold=1.0,logscaling=True,colormap='jet',interpolation='nearest')
	    roi_finder_obj.show_rois()



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


def break_down_det_image(image,pixel_num):
	"""
	Desomposes a Detector image into subimages. Returns a 3D matrix. 
	"""
	# check that there are integer multiples of pixel_num in the big image
	if not image.shape[0] % pixel_num == 0 or not image.shape[1] % pixel_num == 0:
		print 'There must be an integer number of \'pixel_num\' in the large image.'
		return

	detnum_row  = image.shape[0]/pixel_num
	detnum_col  = image.shape[1]/pixel_num
	num_of_dets = detnum_row * detnum_col
	det_images  = np.zeros((num_of_dets,pixel_num,pixel_num))
	x_ranges = []
	for ii in range(detnum_col):
		x_ranges.append((ii*pixel_num, (ii+1)*pixel_num))

	y_ranges = []
	for ii in range(detnum_row):
		y_ranges.append((ii*pixel_num, (ii+1)*pixel_num))

	counter = 0
	offsets = []
	for i in x_ranges:
		for j in y_ranges:
			det_images[counter,:,:]=image[j[0]:j[1],i[0]:i[1]]
			offsets.append((j[0],i[0]))
			counter += 1

	return det_images, offsets

def shift_roi_indices(indices,shift):
	"""
	Applies a given shift (xshift,yshift) to given indices. \
	indices = list of (x,y)-tuples
	shift   = (xshift,yshift) tuple
	"""
	for ind in indices:
		ind[0] += shift[0]
		ind[1] += shift[1]
	return indices

def merge_roi_objects_by_matrix(list_of_roi_objects,large_image_shape,offsets,pixel_num):
	"""
	Merges several roi_objects into one big one using the roi_matrices.
	"""
	# prepare a big roi matrix
	roi_matrix = np.zeros(large_image_shape)
	counter = 0
	for roi_obj in list_of_roi_objects:
		max_number = np.amax(roi_matrix)
		single_matrix = roi_obj.roi_obj.roi_matrix
		inds = single_matrix != 0
		single_matrix[inds] += max_number
		roi_matrix[offsets[counter][0]:offsets[counter][0]+pixel_num,offsets[counter][1]:offsets[counter][1]+pixel_num] = single_matrix
		#roi_obj.roi_obj.roi_matrix + max_number
		counter += 1

	# make a new object to return
	merged_obj = roi_object()
	merged_obj.roi_matrix     = roi_matrix
	merged_obj.indices        = convert_matrix_rois_to_inds(roi_matrix)
	merged_obj.red_rois       = convert_matrix_to_redmatrix(roi_matrix)
	merged_obj.number_of_rois = np.amax(roi_matrix)
	merged_obj.kind           = list_of_roi_objects[0].roi_obj.kind
	merged_obj.x_indices      = convert_inds_to_xinds(merged_obj.indices)
	merged_obj.y_indices      = convert_inds_to_yinds(merged_obj.indices)
	merged_obj.masks          = convert_roi_matrix_to_masks(roi_matrix)

	return merged_obj

def swap_indices_old_rois(old_indices):
	"""
	Swappes x- and y-indices from indices ROIs.
	"""
	new_indices = []
	for roi in old_indices:
		one_new_roi = []
		for point in roi:
			one_new_roi.append((point[1],point[0]))
		new_indices.append(one_new_roi)
	return new_indices

def define_polygon_roi(image_shape,verbose=False):
	"""
	Define a polygon shaped ROI from a current image by 
	selecting points.
	"""

	tmptuple=plt.ginput(0,timeout=-1,show_clicks=True)
	pnts=np.array(tmptuple)
	if verbose:
		print 'The selected points are:' 
		print pnts
	# Create the polygon path from points clicked
	path=Path(pnts)

	# bounding box: everything outside this is necessarily zero
	bbx = np.arange( np.floor(min( pnts[:,0])), np.ceil(max( pnts[:,0])))
	bby = np.arange( np.floor(min( pnts[:,1])), np.ceil(max( pnts[:,1])))
	
	roi=[]
	for ii in bbx:
		for jj in bby:
			# Test point
			if path.contains_point([ii, jj]):
				roi.append( (jj, ii))  

	return roi



#######
#
# for compatibility with old style rois, here is the old rois class (should be deleted asap)
#
#######

class rois:
	"""
	a class to define ROIs 
	ROIs are saved in a container (self.rois)
	"""
	def __init__(self,scans,scannumbers):
		self.scandata = scans # dictionary with instances of the onescan class
		if not isinstance(scannumbers,list):
			thescannumbers = []
			thescannumbers.append(scannumbers)
		else:
			thescannumbers = scannumbers  
		self.scannums   = thescannumbers
		self.DET_PIXEL_NUM = 256
		self.rois       = container() # single variable that holds all relevant info about the defined rois
		self.rois.inds  = [] # list of lists of touples (for each ROI)
		self.rois.xinds = [] # list of lists of x-indices (for each ROI)
		self.rois.yinds = [] # list of lists of y-indices (for each ROI)
		self.rois.roiNumber = 0 # number of ROIs defined
		self.rois.kind    = [] # what kind of ROI (auto, zoom, linear)
		self.rois.height  = [] # list of hights (in pixels, for each ROI)
		self.rois.width   = [] # list of widths (in pixels, for each ROI)

	def preparemats(self):
		"""
		sums and squeezes all edf-files of a scan/all scans into one matrix 
		"""
		# take shape of the first edf-matrices
		dim = np.shape(self.scandata['Scan%03d' % self.scannums[0]].edfmats)
		edfmatrices = np.zeros((1,dim[1],dim[2]))
		for number in self.scannums:
			scanname = 'Scan%03d' % number
			edfmatrices = np.append(edfmatrices,self.scandata[scanname].edfmats,axis=0)
		return sumx(edfmatrices)
	
	def prepareedgemats(self,index):
		"""
		difference between two summed and squeezed edf-files of a scan from below and above energyval 
		"""
		# take shape of the first edf-matrices 
		dim = np.shape(self.scandata['Scan%03d' % self.scannums[0]].edfmats)
		abovemats = np.zeros((1,dim[1],dim[2]))
		belowmats = np.zeros((1,dim[1],dim[2]))
		for number in self.scannums:
			scanname = 'Scan%03d' % number
			abovemats = np.append(abovemats,self.scandata[scanname].edfmats[index:,:,:],axis=0)
			belowmats = np.append(belowmats,self.scandata[scanname].edfmats[:index,:,:],axis=0)
		return np.absolute(np.squeeze(np.sum(abovemats,axis=0))-np.squeeze(np.sum(belowmats,axis=0)))

	def getlinrois(self,numrois,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points
		"""
		# clear all existing rois
		self.deleterois()

		plt.clf()

		fig = plt.figure(figsize=(8, 6))
		ax = fig.add_subplot(111, axisbg='#FFFFCC')

		thematrix = self.preparemats()
		plt.title('Click start and endpoints of the linear ROIs you whish.',fontsize=14)		
		if logscaling:
			theimage = plt.imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = plt.imshow(thematrix)
		# Choose a color palette
		theimage.set_cmap(colormap)
		plt.colorbar()

		therois   = [] # list of lists of tuples with (x,y) coordinates
		theroisx  = [] # list of lists of x-indices
		theroisy  = [] # list of lists of y-indices

		cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

		for m in range(numrois):
			plt.clf()
			if m>0:
				if therois[m-1]:
					for index in therois[m-1]:
						thematrix[index[1],index[0]] = np.amax(thematrix)
			print ("chose two points as endpoints for linear ROI No. %s" % (m+1))

			if logscaling:
				theimage = plt.imshow(np.log(np.absolute(thematrix)))
			else:
				theimage = plt.imshow(thematrix)	
			theimage.set_cmap(colormap)
			plt.draw()	

			endpoints = []
			endpts    = plt.ginput(2)
			for point in endpts:
				point = [round(element) for element in point]
				endpoints.append(np.array(point))
			roix = np.arange(endpoints[0][0],endpoints[1][0])
			roiy = [round(num) for num in np.polyval(np.polyfit([endpoints[0][0],endpoints[1][0]],[endpoints[0][1],endpoints[1][1]],1),roix)]

			theheight = np.arange(-height,height)
			
			eachroi  = []
			eachroix = []
			eachroiy = []
			for n in range(len(roix)):
				for ind in theheight:
					eachroi.append((roix[n],roiy[n]+ind))
					eachroix.append(roix[n])
					eachroiy.append(roiy[n]+ind)
			therois.append(eachroi)
			theroisx.append(eachroix)
			theroisy.append(eachroiy)

			self.rois.height.append(len(theheight)) # save the hight of each roi (e.g. for imaging)
			self.rois.width.append(endpoints[1][0]-endpoints[0][0]) # save the width of each roi
			del endpoints
			del eachroi
			del eachroix
			del eachroiy
			
		self.rois.kind    = 'linear'
		self.rois.inds  = therois
		self.rois.xinds = theroisx
		self.rois.yinds = theroisy
		self.rois.roiNumber = numrois

	def getzoomrois(self,numrois,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on plot
		"""
		# clear all existing rois
		self.deleterois()

		plt.clf()
		thematrix = self.preparemats()
		limmax    = np.shape(thematrix)[1]
		therois   = [] # list of lists of tuples with (x,y) coordinates
		theroisx  = [] # list of lists of x-indices
		theroisy  = [] # list of lists of y-indices
		plt.title('Zoom around ROI, change to shell, and press Return')
		if logscaling:
			theimage = plt.imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = plt.imshow(thematrix)	
		theimage.set_cmap(colormap)
		plt.ion()	
		for m in range(numrois):
			plt.clf()
			if m>0:
				if therois[m-1]:
					for index in therois[m-1]:
						thematrix[index[0],index[1]] = np.amax(thematrix)
			plt.title('Zoom around ROI, change to shell, and press Return')
			if logscaling:
				theimage = plt.imshow(np.log(np.absolute(thematrix)))
			else:
				theimage = plt.imshow(thematrix)	
			theimage.set_cmap(colormap)
			plt.draw()	
						
			thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
			wait = raw_input(thestring)
								
			limits = np.floor(plt.axis())
			
			inds = limits < 1
			limits[inds] = 1
			inds = limits > limmax
			limits[inds] = limmax
			print 'selected limits are: ', limits

			T = np.zeros(np.shape(thematrix))
			T[limits[3]:limits[2],:] = 1
			t = np.zeros(np.shape(thematrix))
			t[:,limits[0]:limits[1]] = 1
			indsx,indsy = np.where(t*T == 1)
			eachroi  = []
			eachroix = []
			eachroiy = []
			for n in range(len(indsx)):
				eachroi.append((indsy[n],indsx[n]))
				eachroix.append(indsx[n])
				eachroiy.append(indsy[n])
			therois.append(eachroi)
			theroisx.append(eachroix)
			theroisy.append(eachroiy)
		plt.show()
		self.rois.kind      = 'zoom'
		self.rois.inds      = therois
		self.rois.xinds     = theroisx
		self.rois.yinds     = theroisy
		self.rois.roiNumber = numrois

	def getzoomrois_frommatrix(self,matrix,numrois,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on plot
		"""
		# clear all existing rois
		self.deleterois()

		plt.clf()
		thematrix = matrix
		limmax = np.shape(thematrix)[1]
		therois = []
		plt.title('Zoom around ROI, change to shell, and press Return')
		if logscaling:
			theimage = plt.imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = plt.imshow(thematrix)	
		theimage.set_cmap(colormap)
		plt.ion()
		for m in range(numrois):
			
			plt.title('Zoom around ROI, change to shell, and press Return')
			if logscaling:
				theimage = plt.imshow(np.log(np.absolute(thematrix)))
			else:
				theimage = plt.imshow(thematrix)	
			theimage.set_cmap(colormap)
			plt.ion()

			thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
			wait = raw_input(thestring)

			limits = np.floor(plt.axis())
			plt.axis([0.0,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM,0.0])

			inds = limits < 1
			limits[inds] = 1
			inds = limits > limmax
			limits[inds] = limmax
			print 'selected limits are: ', limits

			T = np.zeros(np.shape(thematrix))
			T[limits[3]:limits[2],:] = 1
			t = np.zeros(np.shape(thematrix))
			t[:,limits[0]:limits[1]] = 1
			indsx,indsy = np.where(t*T == 1)
			eachroi = []
			for n in range(len(indsx)):
				eachroi.append((indsy[n],indsx[n]))
			therois.append(eachroi)

		plt.show()
		self.rois.kind      = 'zoom'
		self.rois.inds      = therois
		self.rois.roiNumber = numrois

	def getlinedgerois(self,index,numrois=9,logscaling=True,height=5,colormap='jet'):
		"""
		define ROIs by clicking two points in difference picture
		"""
		# clear all existing rois
		self.deleterois()

		plt.clf()
		thematrix = self.prepareedgemats(index) # change this to use an index
		if logscaling:
			theimage = plt.imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = plt.imshow(thematrix)

		# Choose a color palette
		theimage.set_cmap(colormap)
		plt.colorbar()

		therois = []
		for m in range(numrois):
			print ("chose two points as endpoints for linear ROI No. %s" % (m+1))
			endpoints = []
			endpts    = plt.ginput(2)
			for point in endpts:
				point = [round(element) for element in point]
				endpoints.append(np.array(point))
			roix = np.arange(endpoints[0][0],endpoints[1][0])
			roiy = [round(num) for num in np.polyval(np.polyfit([endpoints[0][0],endpoints[1][0]],[endpoints[0][1],endpoints[1][1]],1),roix)]
			del endpoints
			theheight = np.arange(-height,height)
			eachroi = []
			for n in range(len(roix)):
				for ind in theheight:
					eachroi.append((roix[n],roiy[n]+ind))
			therois.append(eachroi)
			del eachroi			
		plt.show()
		self.rois.kind      = 'linear'
		self.rois.inds      = therois
		self.rois.roiNumber = numrois

	def getzoomedgerois(self,index,numrois=9,logscaling=True,colormap='jet'):
		"""
		define ROIs by zooming in on plot
		"""
		# clear all existing rois
		self.deleterois()

		plt.clf()
		thematrix = self.prepareedgemats(index)
		limmax = np.shape(thematrix)[1]
		therois = []
		title('Zoom around ROI, change to shell, and press Return')
		if logscaling:
			theimage = plt.imshow(np.log(np.absolute(thematrix)))
		else:
			theimage = plt.imshow(thematrix)
		theimage.set_cmap(colormap)
		plt.ion()
		for m in range(numrois):

			plt.title('Zoom around ROI, change to shell, and press Return')
			if logscaling:
				theimage = plt.imshow(np.log(np.absolute(thematrix)))
			else:
				theimage = plt.imshow(thematrix)
			theimage.set_cmap(colormap)
			plt.ion()
						
			thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
			wait = raw_input(thestring)
								
			limits = np.floor(axis())
			plt.axis([0.0,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM,0.0])

			inds = limits < 1
			limits[inds] = 1
			inds = limits > limmax
			limits[inds] = limmax
			print 'selected limits are: ', limits

			T = np.zeros(np.shape(thematrix))
			T[limits[3]:limits[2],:] = 1
			t = np.zeros(np.shape(thematrix))
			t[:,limits[0]:limits[1]] = 1
			indsx,indsy = np.where(t*T == 1)
			eachroi = []
			for n in range(len(indsx)):
				eachroi.append((indsy[n],indsx[n]))
			therois.append(eachroi)
		plt.show()
		self.rois.kind      = 'zoom'
		self.rois.inds      = therois
		self.rois.roiNumber = numrois

	def getautorois(self,kernel_size=5,colormap='jet',thresholdfrac=100):
		"""
		define ROIs by choosing a threshold through a slider bar on the plot window
		"""
		# clear all existing rois
		self.deleterois()

		plt.clf()
		thematrix = self.preparemats() # the starting matrix to plot
		plt.title('Crank the threshold and close the plotting window, when you are satisfied.',fontsize=14)		
		ax = plt.subplot(111) 
		plt.subplots_adjust(left=0.05, bottom=0.2)
		thres0 = 0 # initial threshold value
		theimage = plt.imshow(np.log(np.absolute(thematrix)))
		theimage.set_cmap(colormap)
		plt.colorbar()
		thresxcolor = 'lightgoldenrodyellow'
		thresxamp  = plt.axes([0.2, 0.10, 0.55, 0.03], axisbg=thresxcolor)
		sthres = plt.Slider(thresxamp, 'Threshold', 0.0, np.floor(np.amax(thematrix)/thresholdfrac), valinit=thres0)

		def update(val):
			thres = sthres.val
			newmatrix = signal.medfilt2d(thematrix, kernel_size=kernel_size)
			belowthres_indices = newmatrix < thres
			newmatrix[belowthres_indices] = 0			
			labelarray,numfoundrois = measurements.label(newmatrix)
			print str(numfoundrois) + ' ROIs found!'	
			# organize rois
			therois = []
			for n in range(numfoundrois):
				rawindices = np.nonzero(labelarray == n+1)
				eachroi = []				
				for m in range(len(rawindices[0])):
					eachroi.append((rawindices[1][m],rawindices[0][m]))
				therois.append(eachroi)
			self.rois.inds = therois
			self.rois.roiNumber = numfoundrois
			theimage.set_data(newmatrix)
			plt.draw()
		sthres.on_changed(update)
		plt.show()
		self.rois.kind = 'auto'

	def getautorois_eachdet(self,kernel_size=5,colormap='jet',thresholdfrac=100):
		"""
		autoroi, detector for detector for ID20 (6 detector setup) so that different 
		thresholds can be choosen for the different detectors
		define ROIs by choosing a threshold through a slider bar on the plot window
		"""
		# clear all existing rois
		self.deleterois()

		wholematrix = np.array(self.preparemats()) # full 6-detector image
		imagelims = [[0,256,0,256],[0,256,256,512],[0,256,512,768],[256,512,0,256],[256,512,256,512],[256,512,512,768]] 

		for n in range(len(imagelims)):
			plt.clf()
			# thematrix now is a single detector image (256x256 pixels)
			thematrix = wholematrix[imagelims[n][0]:imagelims[n][1],imagelims[n][2]:imagelims[n][3]]
			plt.title('Crank the threshold and close the plotting window, when you are satisfied.',fontsize=14)		
			ax = plt.subplot(111) 
			plt.subplots_adjust(left=0.05, bottom=0.2)
			thres0 = 0 # initial threshold value
			theimage = plt.imshow(np.log(np.absolute(thematrix)))
			theimage.set_cmap(colormap)
			plt.colorbar()
			thresxcolor = 'lightgoldenrodyellow'
			thresxamp  = plt.axes([0.2, 0.10, 0.55, 0.03], axisbg=thresxcolor)
			sthres = plt.Slider(thresxamp, 'Threshold', 0.0, np.floor(np.amax(thematrix)/thresholdfrac), valinit=thres0)

			# using the "container" class to pass the results of the rois from the nested function
			roi_result = container()

			def update(val):
				thres = sthres.val
				newmatrix = signal.medfilt2d(thematrix, kernel_size=kernel_size)
				belowthres_indices = newmatrix < thres
				newmatrix[belowthres_indices] = 0			
				labelarray,numfoundrois = measurements.label(newmatrix)
				print str(numfoundrois) + ' ROIs found!'	
				# organize rois
				therois = []
				for n in range(numfoundrois):
					rawindices = np.nonzero(labelarray == n+1)
					eachroi = []
					for m in range(len(rawindices[0])):
						eachroi.append((rawindices[1][m],rawindices[0][m]))
					therois.append(eachroi)
				theimage.set_data(newmatrix)
				plt.draw()
				roi_result.therois =  therois
				roi_result.numfoundrois = numfoundrois

			sthres.on_changed(update)
			plt.show()
			thestring = 'Press enter to continue.'
			wait = raw_input(thestring)
			self.rois.inds.extend(roi_result.therois)
			self.rois.roiNumber += roi_result.numfoundrois

		self.rois.kind = 'auto'

	def saverois_old(self,filename):
		"""
		save the ROIs in file with name filename
		"""	
		f = open(filename, 'wb')
		theobject = [self.rois, self.roinums, self.roikind]
		pickle.dump(theobject, f,protocol=-1)
		f.close()

	def loadrois_old(self,filename):
		"""
		load ROIs from file with name filename
		"""
		f = open(filename,'rb')
		theobject = pickle.load(f)
		f.close()
		self.rois    = theobject[0]
		self.roinums = theobject[1]
		self.roikind = theobject[2]

	def saverois(self,filename):
		"""
		save the ROIs in file with name filename
		"""	
		f = open(filename, 'wb')
		pickle.dump(self.rois, f,protocol=-1)
		f.close()

	def loadrois(self,filename):
		"""
		load ROIs from file with name filename
		"""
		f = open(filename,'rb')
		self.rois = pickle.load(f)
		f.close()
		
	def deleterois_old(self):
		"""
		delete the existing ROIs
		"""
		self.rois    = []
		self.roinums = []
		self.roikind = []

	def deleterois(self):
		"""
		delete existing ROIs (same as at initialization)
		"""
		self.rois       = container() # single variable that holds all relevant info about the defined rois
		self.rois.inds  = [] # list of lists of touples (for each ROI)
		self.rois.xinds = [] # list of lists of x-indices (for each ROI)
		self.rois.yinds = [] # list of lists of y-indices (for each ROI)
		self.rois.roiNumber = 0 # number of ROIs defined
		self.rois.kind    = [] # what kind of ROI (auto, zoom, linear)
		self.rois.height  = [] # list of hights (in pixels, for each ROI)
		self.rois.width   = [] # list of widths (in pixels, for each ROI)

	def plotrois(self):
		"""
		returns a plot of the ROI shapes
		"""
		plt.figure()
		thematrix = np.zeros((self.DET_PIXEL_NUM,self.DET_PIXEL_NUM))
		for roi in self.rois:
			for index in roi:
				thematrix[index[1],index[0]] = 1
		theimage = plt.imshow(thematrix)
		plt.show()










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
