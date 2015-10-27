#!/usr/bin/python
# Filename: id20_imaging.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import io
from itertools import groupby
from scipy.interpolate import Rbf, RectBivariateSpline
from scipy.optimize import leastsq, fmin

from helpers import *
#from xrs_read import rois

import xrs_read

class oneD_imaging(xrs_read.read_id20):
	""" **oneD_imaging**
	Class to construct images using the 1D piercing mode.
	"""
	def __init__(self,absfilename,energycolumn='sty',monitorcolumn='kapraman',edfName=None,single_image=True):
		try:
			self.path     = os.path.split(absfilename)[0] + '/'
			self.filename = os.path.split(absfilename)[1]
		except IOError:
			print('IOError! No such SPEC file, please check SPEC-filename.')
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
		self.encolumn       = energycolumn.lower()
		self.monicolumn     = monitorcolumn.lower()
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
		self.roi_obj = [] # an instance of the roi_object class from the xrs_rois module (new)
		# images
		self.twoDimages = {} # dictionary: one entry per scan, each scan has a list of one image per ROI

	def loadscan_2Dimages(self, scannumbers,scantype='sty'):
		# check if there are ROIs defined
		if not self.roi_obj:
			print 'Please define some ROIs first.'
			return
		# make sure scannumbers are iterable (list)
		if not isinstance(scannumbers,list):
			scannums = []
			scannums.append(scannumbers)
		else:
			scannums = scannumbers 

		for number in scannums:	# go through all scans
			scanname = 'Scan%03d' % number
			data, motors, counters, edfmats = self.readscan(number)
			position = counters[scantype.lower()]
			self.twoDimages[scanname] = []
			for ii in range(len(self.roi_obj.x_indices)): # go through all ROIs
				# construct the image
				roixinds = self.roi_obj.x_indices[ii]
				roiyinds = self.roi_obj.y_indices[ii]
				# go through all edf files of the scan, sum over the height of the roi and stack the resulting lines into a matrix
				axesrange = [0,roiyinds[-1],position[-1],position[0]]
				imageMat  = (np.sum(edfmats[:,np.amin(roixinds):np.amax(roixinds)+1,np.amin(roiyinds):np.amax(roiyinds)+1],axis=1))
				imageInst = LRimage(imageMat,position,roiyinds)
				self.twoDimages[scanname].append(imageInst)



class imageset:
	"""
	class to make SR-images from list of LR-images
	"""
	def __init__(self):
		self.list_of_images = []
		self.xshifts  = []
		self.yshifts  = []
		self.shifts   = []
		self.srimage  = []
		self.srxscale = []
		self.sryscale = []
		self.refimagenum = []

	def estimate_xshifts(self,whichimage=None):
		if not whichimage:
			ind = 0 # first image in list_of_images is the reference image
		else:
			ind = whichimage
		origx   = self.list_of_images[ind].xscale
		origy   = self.list_of_images[ind].yscale
		origim  = self.list_of_images[ind].matrix
		xshifts = []		
		for image in self.list_of_images:
			newx  = image.xscale
			newy  = image.yscale
			newim = image.matrix
			xshifts.append(estimate_xshift(origx,origy,origim,newx,newy,newim))
		self.refimagenum = ind
		self.xshifts = xshifts

	def estimate_yshifts(self,whichimage=None):
		if not whichimage:
			ind = 0 # first image in list_of_images is the reference image
		else:
			ind = whichimage
		origx   = self.list_of_images[ind].xscale
		origy   = self.list_of_images[ind].yscale
		origim  = self.list_of_images[ind].matrix
		yshifts = []		
		for image in self.list_of_images:
			newx  = image.xscale
			newy  = image.yscale
			newim = image.matrix
			yshifts.append(estimate_yshift(origx,origy,origim,newx,newy,newim))
		self.refimagenum = ind
		self.yshifts = yshifts

	def estimate_shifts(self,whichimage=None):
		if not whichimage:
			ind = 0 # first image in list_of_images is the reference image
		else:
			ind = whichimage
		origx   = self.list_of_images[ind].xscale
		origy   = self.list_of_images[ind].yscale
		origim  = self.list_of_images[ind].matrix
		shifts = []		
		for image in self.list_of_images:
			newx  = image.xscale
			newy  = image.yscale
			newim = image.matrix
			shifts.append(estimate_shift(origx,origy,origim,newx,newy,newim))
		self.refimagenum = ind
		self.shifts = shifts

	def interpolate_xshift_images(self,scaling,whichimages=None):
		if not whichimages:
			inds = range(len(self.list_of_images))
		elif not isinstance(whichimages,list):
			inds = []
			inds.append(whichimages)
		else:
			inds = whichimages

		newim = np.zeros((len(inds),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling,np.shape(self.list_of_images[inds[0]].matrix)[1]))
		newx  = np.linspace(self.list_of_images[inds[0]].xscale[0]-self.xshifts[inds[0]],self.list_of_images[inds[0]].xscale[-1]-self.xshifts[inds[0]],len(self.list_of_images[inds[0]].xscale)*scaling)
		
		newy  = self.list_of_images[inds[0]].yscale

		for n in range(len(inds)):
			print self.xshifts[inds[n]]
			oldim = self.list_of_images[inds[n]].matrix
			oldx  = self.list_of_images[inds[n]].xscale-self.xshifts[inds[n]]
			oldy  = self.list_of_images[inds[n]].yscale
			newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

		self.srimage = np.sum(newim,axis=0)
		self.srxscale = newx
		self.sryscale = newy

	def interpolate_yshift_images(self,scaling,whichimages=None):
		if not whichimages:
			inds = range(len(self.list_of_images))
		elif not isinstance(whichimages,list):
			inds = []
			inds.append(whichimages)
		else:
			inds = whichimages

		newim = np.zeros((len(inds),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling,np.shape(self.list_of_images[inds[0]].matrix)[1]))
		newx  = self.list_of_images[0].xscale
		newy  = np.linspace(self.list_of_images[inds[0]].yscale[0]-self.yshifts[inds[0]],self.list_of_images[inds[0]].yscale[-1]-self.yshifts[inds[0]],len(self.list_of_images[inds[0]].yscale)*scaling)
		for n in range(len(inds)):
			oldim = self.list_of_images[inds[n]].matrix
			oldx  = self.list_of_images[inds[n]].xscale
			oldy  = self.list_of_images[inds[n]].yscale-self.yshifts[inds[n]]
			newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

		self.srimage = np.sum(newim,axis=0)
		self.srxscale = newx
		self.sryscale = newy

	def interpolate_shift_images(self,scaling,whichimages=None):
		if not whichimages:
			inds = range(len(self.list_of_images))
		elif not isinstance(whichimages,list):
			inds = []
			inds.append(whichimages)
		else:
			inds = whichimages

		if len(scaling)<2:
			scaling = [scaling, scaling]
		print inds, self.list_of_images[inds[0]].xscale[0], self.shifts[inds[0]], self.list_of_images[inds[0]].xscale[-1]
		newim = np.zeros((len(self.list_of_images),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling[0],np.shape(self.list_of_images[inds[0]].matrix)[1]*scaling[1]))
		newx  = np.linspace(self.list_of_images[inds[0]].xscale[0]-self.shifts[inds[0]][0],self.list_of_images[inds[0]].xscale[-1]-self.shifts[inds[0]][0],len(self.list_of_images[inds[0]].xscale)*scaling[0])
		newy  = np.linspace(self.list_of_images[inds[0]].yscale[0]-self.shifts[inds[0]][1],self.list_of_images[inds[0]].yscale[-1]-self.shifts[inds[0]][1],len(self.list_of_images[inds[0]].yscale)*scaling[1])

		for n in range(len(inds)):
			oldim = self.list_of_images[inds[n]].matrix
			oldx  = self.list_of_images[inds[n]].xscale-self.shifts[inds[n]][0]
			oldy  = self.list_of_images[inds[n]].yscale-self.shifts[inds[n]][1]
			newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

		self.srimage = np.sum(newim,axis=0)
		self.srxscale = newx
		self.sryscale = newy

	def plotSR(self):
		X,Y = pylab.meshgrid(self.srxscale,self.sryscale)
		pylab.pcolor(X,Y,np.transpose(self.srimage))
		pylab.show(block=False)

	def plotLR(self,whichimage):
		if not isinstance(whichimage,list):
			inds = []
			inds.append(whichimage)
		else:
			inds = list(whichimage)

		for ind in inds:
			X,Y = pylab.meshgrid(self.list_of_images[ind].xscale,self.list_of_images[ind].yscale)
			pylab.figure(ind)
			pylab.pcolor(X,Y,np.transpose(self.list_of_images[ind].matrix))
			pylab.show(block=False)

	def save(self):
		pass

	def load(self):
		pass

	def loadkimberlite(self,matfilename):
		data_dict   = io.loadmat(matfilename)
		sorted_keys = sorted(data_dict.keys())
		sy    = data_dict['sy'][0]
		allsx = []
		for key in sorted_keys[3:12]:
			allsx.append(data_dict[key][0])

		allmats = []
		for key in sorted_keys[13:]:
			allmats.append(data_dict[key])

		alllengths = []
		for sx in allsx:
			alllengths.append(len(sx))
		
		ind = np.where(alllengths == np.max(alllengths))[0][0]
		# spline everything onto longest sx-scale
		for n in range(len(allmats)):
			print np.shape(allsx[n]), np.shape(sy), np.shape(allmats[n])
			ip         = RectBivariateSpline(allsx[n],sy,allmats[n])
			allmats[n] = ip(allsx[ind],sy)
			allsx[n]   = allsx[ind]

		allimages = []
		for n in range(len(allmats)):
			allimages.append(image(allmats[n],allsx[n],sy))
		self.list_of_images = allimages

	def loadhe3070(self,matfilename):
		data_dict  = io.loadmat(matfilename)
		sy    = data_dict['det'][0][0]['sy'][0][0][0]
		allsx = []
		allmats = []
		for n in range(9):
			allsx.append(np.reshape(data_dict['det'][0][n]['sx'][0][0]-data_dict['det'][0][n]['sx'][0][0][0],len(data_dict['det'][0][n]['sx'][0][0],)))
			allmats.append(data_dict['det'][0][n]['img'][0][0])

		alllengths = []
		for sx in allsx:
			alllengths.append(len(sx))

		ind = np.where(alllengths == np.max(alllengths))[0][0]
		for n in range(len(allmats)):
			print np.shape(allsx[n]), np.shape(sy), np.shape(np.transpose(allmats[n]))
			ip         = RectBivariateSpline(allsx[n],sy,np.transpose(allmats[n]))
			allmats[n] = ip(allsx[ind],sy)
			allsx[n]   = allsx[ind]

		allimages = []
		for n in range(len(allmats)):
			allimages.append(image(allmats[n],allsx[n],sy))
		self.list_of_images = allimages

class image:
	"""
	Container class to hold info of a single LR-image to be put togther in a SR-image by the imageset class
	"""
	def __init__(self,matrix,xscale,yscale):
		self.name	= []
		self.matrix	= matrix
		self.xscale	= xscale
		self.yscale	= yscale
		self.tth	= []

	def plotimage(self):
		pass

	def shiftx(self):
		pass

	def shifty(self):
		pass

	def save(self):
		pass

	def load(self):
		pass


def interpolate_image(oldx,oldy,oldIM,newx,newy):
	"""
	2d interpolation
	"""
	interp = RectBivariateSpline(oldx,oldy,oldIM)
	return interp(newx,newy)

def estimate_xshift(x1,y1,im1,x2,y2,im2):
	"""
	estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
	"""
	funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2+a,y2))**2.0)
	res = leastsq(funct,0.0)
	return res[0]

def estimate_yshift(x1,y1,im1,x2,y2,im2):
	"""
	estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
	"""
	funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2,y2+a))**2.0)
	res = leastsq(funct,0.0)
	return res[0]

def estimate_shift(x1,y1,im1,x2,y2,im2):
	"""
	estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
	"""
	funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2+a[0],y2+a[1]))**2.0)
	res = fmin(funct,[0.0,0.0],disp=0)
	return res


class LRimage:
	"""
	container class to hold info of a single LR-image to be put togther in a SR-image by the imageset class
	"""
	def __init__(self,matrix,xscale,yscale):
		self.name	= []
		self.matrix	= matrix
		self.xscale	= xscale
		self.yscale	= yscale
		self.tth	= []

	def plotimage(self):
		pass

	def shiftx(self):
		pass

	def shifty(self):
		pass

	def save(self):
		pass

	def load(self):
		pass

