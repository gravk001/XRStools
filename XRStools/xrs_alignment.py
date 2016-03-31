#!/usr/bin/python
# Filename: xrs_alignment.py

#/*##########################################################################
#
# The XRStools software package for XRS spectroscopy
#
# Copyright (c) 2013-2014 European Synchrotron Radiation Facility
#
# This file is part of the XRStools XRS spectroscopy package developed at
# the ESRF by the DEC and Software group and contains practical functions, 
# most of which are translated from Matlab functions from the University of
# Helsinki Electronic Structure Laboratory.
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

import math_functions
import numpy as np
from scipy import interpolate, optimize

def findAnalyzerFocus(edfmats,motorScale,roi_obj,roiNumber):
	""" **findAnalyzerFocus**

	Returns motor position that optimizes the analyzer focus subject to a 2D Gaussian fit.

	Args: 
	-----
	edfmats (np.array): 3D Numpy array containing the EDF-matrices of a scan.
	motorScale (np.array): Motor positions along the scan (analyzer tx-scan).
	roi_obj (xrs_rois.roi_object): ROI object, defined from the scan, should have some margins around the spots.
	roiNumber (int): Number indicating which ROI to be optimized.

	Returns:
	--------
	optPos (float): Motor position that optimizes the focus of the analyzer in question.

	"""
	xmin = min(roi_obj.x_indices[roiNumber])
	xmax = max(roi_obj.x_indices[roiNumber])
	ymin = min(roi_obj.y_indices[roiNumber])
	ymax = max(roi_obj.y_indices[roiNumber])
	x = np.arange(xmin,xmax)
	y = np.arange(ymin,ymax)
	xx, yy = np.meshgrid(y, x)
	sigma_1 = [] # FWHM along dimension 1
	sigma_2 = [] # FWHM along dimension 2
	# go through all images and fit Gaussian
	for ii in range(edfmats.shape[0]):
		initial_guess = (np.amax(edfmats[ii,xmin:xmax,ymin:ymax]),(ymin+ymax)/2.,(xmin+xmax)/2.,1.0,1.0,0)
		popt, pcov = optimize.curve_fit(math_functions.flat2DGaussian, (xx, yy), edfmats[ii,xmin:xmax,ymin:ymax].ravel(), p0=initial_guess)
		sigma_1.append(popt[3])
		sigma_2.append(popt[4])

	return sigma_1, sigma_2










