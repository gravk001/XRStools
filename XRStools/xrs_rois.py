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
import copy

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
		self.number_of_rois = 0  # number of ROIs defined
		self.kind           = [] # keyword (e.g. 'zoom', 'line', 'auto', etc.), certain features (esp. in imaging) are only available for certain kinds of ROIs
		self.x_indices      = [] # list of numpy arrays of x-indices (for each ROI)
		self.y_indices      = [] # list of numpy arrays of y-indices (for each ROI)
		self.masks          = [] # 3D numpy array with slices of zeros and ones (same size as detector image) for each roi
		self.input_image	= [] # 2D imput image that was used to define the ROIs
		
	def load_rois_fromMasksDict(self, masksDict, newshape=None, kind="zoom"):
		self.kind=kind
		self.red_rois = masksDict
		if newshape is not None:
			self.roi_matrix = np.zeros(newshape)
		self.roi_matrix = convert_redmatrix_to_matrix( masksDict,self.roi_matrix , offsetX=0, offsetY=0)

		self.masks          = convert_roi_matrix_to_masks(self.roi_matrix)
		self.indices        = convert_matrix_rois_to_inds(self.roi_matrix)
		self.number_of_rois = np.amax(self.roi_matrix)
		self.x_indices      = convert_inds_to_xinds(self.indices)
		self.y_indices      = convert_inds_to_yinds(self.indices)

	def save_rois(self,filename):
		"""
		Saves the  roi_bj in a pick
		"""
		import pickle
		f = open(filename, 'wb')
		pickle.dump(self, f,protocol=-1)
		f.close()

	def loadrois(self,filename):
		"""
		loads a file written with the saverois-function using pickle
		filename = absolute path to file and filename
		"""
		import pickle
		f = open(filename,'rb')
		roi_obj = pickle.load(f)
		f.close()
		self.roi_matrix = roi_obj.roi_matrix
		self.red_rois       = roi_obj.red_rois   
		self.indices        = roi_obj.indices  
		self.number_of_rois = roi_obj.number_of_rois
		self.kind           = roi_obj.kind 
		self.x_indices      = roi_obj.x_indices
		self.y_indices      = roi_obj.y_indices 
		self.masks          = roi_obj.masks 
		self.input_image	= roi_obj.input_image

		#print len(roiob.indices), roiob.roi_matrix.shape
		
		#indices = swap_indices_old_rois(roiob.indices)
		#roi_matrix     = convert_inds_to_matrix(indices,roiob.input_image.shape)
		#red_rois       = convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
		#self.load_rois_fromMasksDict(red_rois)


	def load_rois(self,filename):
		"""
		Loads ROIs from a file written with the save_rois-function.
		filename = absolute path to file and filename
		"""
		import pickle
		f = open(filename,'rb')
		roi_obj = pickle.load(f)
                self.load_rois_fromMasksDict(roi_obj.red_rois)
		f.close()

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

	def get_copy(self):
		"""
		**get_copy**
		Returns a deep copy of self.
		"""
		return copy.deepcopy(self)

	def shift_rois(self,shiftVal,direction='horiz',whichroi=None):
		"""
		**shift_rois**
		Displaces the defined ROIs by the provided value.

		Args
		----
		shiftVal : int
			Value by which the ROIs should be shifted.
		direction : string
			Description of which direction to shit by.
		whichroi : sequence
			Sequence (iterable) for which ROIs should be shifted.
		"""
		the_indices = []

		if not whichrois:
			inds = range(len(self.indices))
		else:
			inds = whichroi

		if direction == 'vert':
			for roi in self.indices:
				oneroi = []
				for pixel in roi:
					oneroi.append( (pixel[0]+shiftVal,pixel[1]) )
				the_indices.append(oneroi)

		if direction == 'horiz':
			for roi in self.indices:
				oneroi = []
				for pixel in roi:
					oneroi.append( (pixel[0], pixel[1]+shiftVal) )
				the_indices.append(oneroi)

		self.indices = the_indices

		self.roi_matrix     = convert_inds_to_matrix(self.indices,self.input_image.shape)
		self.red_rois       = convert_matrix_to_redmatrix(self.roi_matrix)
		self.x_indices      = convert_inds_to_xinds(self.indices)
		self.y_indices      = convert_inds_to_yinds(self.indices)
		self.masks          = convert_roi_matrix_to_masks(self.roi_matrix)

def convert_redmatrix_to_matrix( masksDict,mask, offsetX=0, offsetY=0):
    for key, (pos,M)  in masksDict.iteritems():
        num=int("".join([c for c in key if c.isdigit()]))
        S = M.shape
        inset =    (slice(offsetY+pos[0]  , offsetY+pos[0]+S[0]   ), slice(  offsetX+pos[1]  , offsetX+pos[1]+S[1] ) )
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

def  load_rois_fromh5(h5group_tot,md):
    
    h5group = h5group_tot["rois_definition/rois_dict"]
    for key in h5group.keys():
        md[key]=[]
        md[key].append(h5group[key]["origin"][:])
        md[key].append(h5group[key]["mask"][:])

    h5data = h5group_tot["rois_definition/image"]
    shape = h5data.shape

    return shape

def  write_rois_toh5(h5group,md):
    for key in md.keys():
        h5group.require_group(key)
        h5group[key]["origin"]=md[key][0]
        h5group[key]["mask"]=md[key][1]
    
