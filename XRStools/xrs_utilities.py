#!/usr/bin/python
# Filename: xrs_utilities.py

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

import os
import math

import numpy as np
import array as arr
import matplotlib.pyplot as plt
import pickle
import traceback
import sys

from matplotlib.widgets import Cursor
from itertools import groupby
from scipy.integrate import trapz
from scipy import interpolate, signal, integrate, constants, optimize
from re import findall
from scipy.ndimage import measurements
from scipy.optimize import leastsq, fmin
from scipy.interpolate import Rbf, RectBivariateSpline
from scipy.integrate import odeint

# data_installation_dir = os.path.join( os.path.dirname(os.path.abspath(__file__)),"..","..","..","..","share","xrstools","data")
data_installation_dir = os.path.abspath('.')
# os.path.join(getattr(install_cmd, 'install_lib'),"xrstools"+version,"..","..","..","..","share","xrstools","data")


class maxipix_det:
	"""
	Class to store some useful values from the detectors used. To be used for arranging the ROIs.
	"""
	def __init__(self,name,spot_arrangement):
		self.name = name
		assert spot_arrangement in ['3x4','vertical'], 'unknown ROI arrangement, select \'3x4\' or \'vertical\'.'
		self.spot_arrangement = spot_arrangement

		self.roi_indices   = []
		self.roi_x_indices = []
		self.roi_y_indices = []
		self.roi_x_means   = []
		self.roi_y_means   = []
		self.pixel_size    = [256,256]
		self.PIXEL_RANGE   = {'VD': [0,256,0,256],  'VU': [0,256,256,512],  'VB': [0,256,512,768],
							'HR': [256,512,0,256],'HL': [256,512,256,512],'HB': [256,512,512,768]}
	def get_pixel_range(self):
		return self.PIXEL_RANGE[self.name]

	def get_det_name(self):
		return self.name


def energy(d,ba):
	"""
	% ENERGY  Calculates energy corrresponing to Bragg angle for given d-spacing
	%         function e=energy(dspace,bragg_angle)
	%
	%	  dspace for reflection
	%	  bragg_angle in DEG
	%
	%         KH 28.09.93
	"""
	hc = 12.3984191 # CODATA 2002 physics.nist.gov/constants
	return (2.0*d*np.sin(ba/180.0*np.pi)/hc)**(-1)

def dspace(hkl=[6,6,0],xtal='Si'):
	"""
	% DSPACE Gives d-spacing for given xtal
	%	 d=dspace(hkl,xtal)
	%	 hkl can be a matrix i.e. hkl=[1,0,0 ; 1,1,1];
	%	 xtal='Si','Ge','LiF','InSb','C','Dia','Li' (case insensitive)
	%	 if xtal is number this is user as a d0
	%
	%	 KH 28.09.93 
	%        SH 2005
	%
	"""
	# create a database of lattice constants (could be a shelf)
	xtable = {}
	xtable['SI'] = 5.43102088
	xtable['GE'] = 5.657
	xtable['SIXOP'] = 5.430919
	xtable['SIKOH'] = 5.430707
	xtable['LIF'] = 4.027
	xtable['INSB'] = 6.4784
	xtable['C'] = 6.708
	xtable['DIA'] = 3.57
	xtable['LI'] = 3.41
	xtable['TCNE'] = 9.736
	xtable['CU'] = 3.61
	xtable['PB'] = 4.95
	xtable['NA'] = 4.2906
	xtable['AL'] = 4.0495

	if isinstance(xtal,str):
		try:
			a0 = xtable[xtal.upper()]
		except KeyError:
			print 'Lattice constant is not in database'
			return
	else: 
		a0 = xtal # if number is provided, it's taken as lattice constant

	return a0/np.sqrt(np.sum(np.array(hkl)**2.0))

def bragg(hkl,e,xtal='Si'):
	"""
	% BRAGG  Calculates Bragg angle for given reflection in RAD
	%	  output=bangle(hkl,e,xtal)
	%  	  hkl can be a matrix i.e. hkl=[1,0,0 ; 1,1,1];
	%	  e=energy in keV
	%	  xtal='Si', 'Ge', etc. (check dspace.m) or d0 (Si default)
	%
	%	  KH 28.09.93
	%
	"""
	hc = 12.3984191 # CODATA 2002 recommended value, physics.nist.gov/constants
	return np.real(np.arcsin((2.0*dspace(hkl,xtal)*e/hc)**(-1.0)))

def braggd(hkl,e,xtal='Si'):
	"""
	# BRAGGD  Calculates Bragg angle for given reflection in deg
	#	  Call BRAGG.M
	#	  output=bangle(hkl,e,xtal)
	#  	  hkl can be a matrix i.e. hkl=[1,0,0 ; 1,1,1];
	#	  e=energy in keV
	#	  xtal='Si', 'Ge', etc. (check dspace.m) or d0 (Si default)
	#
	#	  KH 28.09.93
	"""
	return bragg(hkl,e,xtal)/np.pi*180.0

def addch(xold,yold,n,n0=0,errors=None):
	"""
	# ADDCH     Adds contents of given adjacent channels together
	#
	#           [x2,y2] = addch(x,y,n,n0)
	#           x  = original x-scale  (row or column vector)
	#           y  = original y-values (row or column vector)
	#           n  = number of channels to be summed up
	#	        n0 = offset for adding, default is 0
	#           x2 = new x-scale 
	#           y2 = new y-values
	#
	#           KH 17.09.1990
	#	    Modified 29.05.1995 to include offset
	"""
	n0=n0-np.fix(n0/n)*n
	if n0<0:
		 n0 = (n + n0)
	datalen = np.floor( (len(xold) - n0) / n)

	xnew = np.zeros(np.min([datalen,len(xold)]))
	ynew = np.zeros(np.min([datalen,len(xold)]))
	errnew = np.zeros(np.min([datalen,len(xold)]))
	for i in range(int(datalen)):
		xnew[i] = np.sum(xold[i*n+n0:i*n+n+n0])/n
		ynew[i] = np.sum(yold[i*n+n0:i*n+n+n0])/n
		if np.any(errors):
			errnew[i] = np.sqrt(np.sum(errors[i*n+n0:i*n+n+n0]**2.0))
			return xnew, ynew, errnew
	return xnew, ynew

def fwhm(x,y):
	"""
	finds full width at half maximum of the curve y vs. x
	returns 
	f  = FWHM
	x0 = position of the maximum
	"""
	if x[-1] < x[0]:
		x = np.flipud(x)
		y = np.flipud(y)

	y0 = np.amax(y)
	i0 = np.where(y == y0)
	x0 = x[i0]

	i1 = np.where(np.logical_and(y>y/3.0, x<x0))[0]
	i2 = np.where(np.logical_and(y>y/3.0, x>x0))[0]

	if len(y[i1])==0 or len(y[i2])==0:
		return 0,0
	#f  = interpolate.interp1d(y[i1],x[i1], bounds_error=False, fill_value=0.0)
	#x1 = f(y0/2.0)
	#f  = interpolate.interp1d(y[i2],x[i2], bounds_error=False, fill_value=0.0)
	#x2 = f(y0/2.0)
	x1 = np.interp(y0/2.0,y[i1],x[i1])
	x2 = np.interp(y0/2.0,y[i2],x[i2])
	fwhm = x2 - x1
	x0 = np.mean([x2, x1])
	return 2.0*fwhm, x0

def gauss(x,x0,fwhm):
    # area-normalized gaussian
    sigma = fwhm/(2*np.sqrt(2*np.log(2)));
    y = np.exp(-(x-x0)**2/2/sigma**2)/sigma/np.sqrt(2*np.pi)
    return y

def convg(x,y,fwhm):
	"""
	Convolution with Gaussian
	x  = x-vector
	y  = y-vector
	fwhm = fulll width at half maximum of the gaussian with which y is convoluted
	"""
	dx = np.min(np.absolute(np.diff(x)))
	x2 = np.arange(np.min(x)-1.5*fwhm, np.max(x)+1.5*fwhm, dx)
	xg = np.arange(-np.floor(2.0*fwhm/dx)*dx, np.floor(2.0*fwhm/dx)*dx, dx)
	yg = gauss(xg,[0,fwhm])
	yg = yg/np.sum(yg)
	y2 = spline2(x,y,x2)
	c  = np.convolve(y2,yg, mode='full')
	n  = int( np.floor(np.max(np.shape(xg))/2))
	c  = c[n:len(c)-n+1] # not sure about the +- 1 here
	f  = interpolate.interp1d(x2,c)
	return f(x)

def spline2(x,y,x2):
	"""
	Extrapolates the smaller and larger valuea as a constant
	"""
	xmin = np.min(x)
	xmax = np.max(x)
	imin = x == xmin
	imax = x == xmax
	f  = interpolate.interp1d(x,y, bounds_error=False, fill_value=0.0)
	y2 = f(x2)
	i     = np.where(x2<xmin)
	y2[i] = y[imin]
	i     = np.where(x2>xmax)
	y2[i] = y[imax]
	return y2

def pz2e1(w2,pz,th):
	"""Calculates the incident energy for a specific scattered photon and momentum value.

	Returns the incident energy for a given photon energy and scattering angle.
	This function is translated from Keijo Hamalainen's Matlab implementation (KH 29.05.96).

	Args:
	-----
	w2 (float): scattered photon energy in [keV]
	pz (np.array): pz scale in [a.u.]
	th (float): scattering angle two theta in [deg]

	Returns:
	--------
	w1 (np.array): incident energy in [keV]
	"""
	pz  = np.array(pz)
	w   = np.array(np.arange(w2/4.0,4.0*w2,w2/5000.0))
	p   = e2pz(w,w2,th)[0]
	tck = interpolate.UnivariateSpline(p,w)
	w1  = tck(pz)
	return w1

def e2pz(w1,w2,th):
	"""Calculates the momentum scale and the relativistic Compton cross section 
	correction according to P. Holm, PRA 37, 3706 (1988).

	This function is translated from Keijo Hamalainen's Matlab implementation (KH 29.05.96).

	Args:
	-----
	w1 (float or np.array): incident energy in [keV]
	w2 (float or np.array): scattered energy in [keV]
	th (float): scattering angle two theta in [deg]
	returns:
	pz (float or np.array): momentum scale in [a.u.]
	cf (float or np.array): cross section correction factor such that:
	    J(pz) = cf * d^2(sigma)/d(w2)*d(Omega) [barn/atom/keV/srad]
	"""
	w1  = np.array(w1)    # make sure arrays are used
	w2  = np.array(w2)           
	m   = constants.value('electron mass energy equivalent in MeV')*1e3 #511.003      # Mass in natural units
	th  = math.radians(th) # th/180.0*np.pi  # Angle to radians
	alp = constants.value('fine-structure constant') #1.0/137.036  # Fine structure constant
	r0  = constants.value('classical electron radius') #2.8179e-15   # Electron radius
	q   = np.sqrt(w1**2.0 + w2**2.0-2.0*w1*w2*np.cos(th))                        # Momentum transfer    
	pz  = q/2.0 - (w1-w2) * np.sqrt(1.0/4.0 + m**2.0/(2.0*w1*w2*(1.0-np.cos(th)))) # In natural units
	E   = np.sqrt(m**2.0+pz**2.0)
	A   = ((w1-w2)*E-w1*w2*(1.0-np.cos(th)))/q
	D   = (w1-w2*np.cos(th))*A/q
	R   = w1*(E-D)
	R2  = R-w1*w2*(1-np.cos(th))
	chi = R/R2 + R2/R + 2.0*m**2.0 * (1.0/R-1.0/R2) + m**4.0 * (1.0/R-1.0/R2)**2.0
	cf  = 2.0*w1*q*E/(m**2.0*r0**2.0*w2*chi)
	cf  = cf*(1.0e-28*(m*alp)) # Cross section now in barns/atom/keV/srad
	pz  = pz/(m*alp)           # pz to atomic units (a.u.)
	return pz, cf

def momtrans_au(e1,e2,tth):
	""" Returns the momentum transfer (in a.u.).

	Calculates the momentum transfer in atomic units for two given
	energies e1 and e1 (in keV) and the scattering angle tth (two theta).

	Args:
	----- 
	e1 (float or np.array): incident energy in [keV], can be a single value or a vector
	e2 (float or np.array): scattered energy in [keV], can be a single value or a vector
	tth (float): scattering angle two theta in [deg]

	Returns:
	--------
	q (float or np.array): momentum transfer [a.u.], single value or vector depending on input
	"""
	e1    = np.array(e1*1.0e3/13.60569172/2.0)
	e2    = np.array(e2*1.0e3/13.60569172/2.0)
	th    = math.radians(tth)#tth/180.0*np.pi
	hbarc = 137.03599976
	q     = 1/hbarc*np.sqrt(e1**2.0+e2**2.0-2.0*e1*e2*np.cos(th));
	return q

def vrot(v,vaxis,phi):
	""" **vrot**
	Rotates a vector around a given axis.

	Args:
	-----
	v (np.array): vector to be rotated
	vaxis (np.array): rotation axis
	phi (float): angle [deg] respecting the right-hand rule 

	Returns:
	--------
	v2 (np.array): new rotated vector

	Function by S. Huotari (2007) adopted to Python.
	"""
	h = vaxis[0]
	k = vaxis[1]
	l = vaxis[2]
	alpha = np.arctan2(k,h)
	if np.absolute(alpha)>np.finfo(float).eps:
		h2 = np.cos(alpha)*(h+k*np.tan(alpha))
	else:
		h2 = h
	v2 = np.array([h2, 0.0, l])
	ca = np.cos(alpha)
	sa = np.sin(alpha)
	R1 = np.array([[ca, sa, 0.0], [-sa, ca, 0.0], [0.0, 0.0, 1.0]])
	beta = np.radians(vangle(v2,np.array([0.0, 0.0, 1.0])))
	cb = np.cos(beta)
	sb = np.sin(beta)
	R2 = np.array([[cb, 0.0, -sb], [0.0, 1.0, 0.0], [sb, 0.0, cb]])
	phi = np.radians(phi)
	cp = np.cos(phi)
	sp = np.sin(phi)
	R3 = np.array([[cp, -sp, 0.0], [sp, cp, 0.0], [0.0, 0.0, 1.0]])
	v2 = np.dot(R3,np.dot(R2,np.dot(R1,v))) 
	v2 = np.dot(np.linalg.inv(R1),np.dot(np.linalg.inv(R2),v2))
	return v2	

def vangle(v1, v2):
	""" **vangle**
	Calculates the angle between two cartesian vectors v1 and v2 in degrees.

	Args:
	-----
	v1 (np.array): first vector.
	v2 (np.array): second vector.

	Returns:
	--------
	th (float): angle between first and second vector.

	Function by S. Huotari, adopted for Python.
	"""
	return np.arccos(np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2))/np.pi*180.0;

def cixs(theta,phi,G=np.array([-1.0,-1.0,-1.0]),x_vec=np.array([0.0, -1.0, 1.0]),crystal='Si',analRefl=[4,4,4]):
	""" **cixs**
	Calculates q1 for given Theta and Phi angle for coherent IXS experiments.

	Args:
	-----
	theta (float): scattering angle.
	phi (float): rotation angle in the plane perpendicular to G-vector.
	G (np.array): G-vector.
	x_vec (np.array): vector in the plane perpendicular to G-vector.
	crystal (str): wich crystal to use,
	analRefl (list): which analyzer reflection to use

	Returns:
	--------
	kin (np.array):  incident beam vector (K_0)
	kout (np.array): Bragg diffracted beam vector (K_h)
	qout2 (np.array): inelastically scattered beam vector (K')
	q1 (np.array): momentum transfer 1
	q2 (np.array): momentum transfer 2
	"""
	hc = 12.3984191 # CODATA 2002 recommended value, physics.nist.gov/constants
	zz = G
	G = 2.0*np.pi*G/dspace([1.0, 0.0, 0.0])
	xx = vrot(x_vec,G,theta)
	yy = vrot(xx,zz,90.0)
	a  = dspace([1, 0, 0],xtal=crystal)
	Eout = energy(dspace(analRefl),88.0)
	lambdaout = hc/Eout
	E = Eout+0.02
	lam = hc/E
	kin  = vrot(xx,yy,braggd(G,E))
	kin  = kin/np.linalg.norm(kin)*2.0*np.pi/lam
	kout =-vrot(kin,G,180.0)
	qout2 = vrot(kin,yy,-braggd(G,E))
	qout2 = qout2/np.linalg.norm(qout2)*2.0*np.pi/lambdaout
	qout2 = vrot(qout2,G,phi)
	q1 = kin-qout2
	q2 = kout-qout2
	return kin, kout, qout2, q1, q2, G

def cixs2(theta,phi,xi,G=np.array([-1.0,-1.0,-1.0]),x_vec=np.array([0.0, -1.0, 1.0]),crystal='Si',analRefl=[4,4,4]):
	""" **cixs2**
	Calculates q1 and q2 for given Theta, Phi, and Xi angle for coherent IXS experiments with asymmetric q1, q2.

	Args:
	-----
	theta (float): rotation angle around G.
	phi (float): rotation angle around vector in plane perpendicular to G-vector.
	xi (float): rotation angle for asymmetric q1, q2.
	G (np.array): G-vector.
	x_vec (np.array): vector in the plane perpendicular to G-vector.
	crystal (str): wich crystal to use,
	analRefl (list): which analyzer reflection to use

	Returns:
	--------
	kin (np.array):  incident beam vector (K_0)
	kout (np.array): Bragg diffracted beam vector (K_h)
	qout2 (np.array): inelastically scattered beam vector (K')
	q1 (np.array): momentum transfer 1
	q2 (np.array): momentum transfer 2
	"""
	hc = 12.3984191 # CODATA 2002 recommended value, physics.nist.gov/constants
	zz = G
	xx = vrot(x_vec,G,theta)
	yy = vrot(xx,zz,90.0)
	a  = dspace([1, 0, 0],xtal=crystal)
	Eout = energy(dspace(analRefl),89.0)
	lambdaout = hc/Eout
	E = Eout+0.02
	lam = hc/E
	kin  = vrot(xx,yy,braggd(G,E))
	kin  = kin/np.linalg.norm(kin)*2.0*np.pi/lam
	kout =-vrot(kin,G,180.0)
	qout2 = vrot(kin,yy,-braggd(G,E)+xi)
	qout2 = qout2/np.linalg.norm(qout2)*2.0*np.pi/lambdaout
	qout2 = vrot(qout2,G,phi)
	q1 = kin-qout2
	q2 = kout-qout2
	return kin, kout, qout2, q1, q2

def readbiggsdata(filename,element):
	"""
	Reads Hartree-Fock Profile of element 'element' from values tabulated 
	by Biggs et al. (Atomic Data and Nuclear Data Tables 16, 201-309 (1975))
	as provided by the DABAX library (http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/ComptonProfiles.dat).
	input:
	filename = path to the ComptonProfiles.dat file (the file should be distributed with this package)
	element  = string of element name
	returns:
	data     = the data for the according element as in the file:
		#UD  Columns: 
		#UD  col1: pz in atomic units 
		#UD  col2: Total compton profile (sum over the atomic electrons
		#UD  col3,...coln: Compton profile for the individual sub-shells
	occupation = occupation number of the according shells
	bindingen  = binding energies of the accorting shells
	colnames   = strings of column names as used in the file
	"""
	elementid = '#S'
	sizeid    = '#N'
	occid     = '#UOCCUP'
	bindingid = '#UBIND'
	colnameid = '#L'
	data = []
	f = open(filename,'r')
	istrue = True
	while istrue:
		line = f.readline()
		if line[0:2] == elementid:
			if line.split()[-1] == element:
				line = f.readline()
				while line[0:2] != elementid:
					if line[0:2] == sizeid:
						arraysize = int(line.split()[-1])
						line = f.readline()
					if line[0:7] == occid:
						occupation = line.split()[1:]
						line = f.readline()
					if line[0:6] == bindingid:
						bindingen = line.split()[1:]	
						line = f.readline()
					if line[0:2] == colnameid:
						colnames = line.split()[1:]
						line = f.readline()
					if line[0]== ' ':
						data.append([float(n) for n in line.strip().split()])
						#data = np.zeros((31,arraysize))
						line = f.readline()
				break
	length = len(data)
	data = (np.reshape(np.array(data),(length,arraysize)))
	return data, occupation, bindingen, colnames

def makepzprofile(element,filename=os.path.join(data_installation_dir,'data/ComptonProfiles.dat')):
	"""
	constructs compton profiles of element 'element' on pz-scale 
	(-100:100 a.u.) from the Biggs tables provided in 'filename'
	input:
	element   = element symbol (e.g. 'Si', 'Al', etc.)
	filename  = path and filename to tabulated profiles
	returns:
	pzprofile = numpy array of the CP:
		1. column: pz-scale
		2. ... n. columns: compton profile of nth shell
	binden     = binding energies of shells
	occupation = number of electrons in the according shells
	"""
	theory,occupation,binden,colnames = readbiggsdata(filename,element)
	# first spline onto a rough grid:
	roughpz = np.logspace(0.01,2,65)-1
	roughtheory      = np.zeros((len(roughpz),len(binden)+2))
	roughtheory[:,0] = roughpz     
	for n in range(len(binden)+1):
		intf               = interpolate.interp1d(theory[:,0],theory[:,n+1])
		roughtheory[:,n+1] = intf(roughpz)
	pzscale   = np.linspace(-100,100,num=4000)
	pzprofile      = np.zeros((len(pzscale),len(binden)+1))
	pzprofile[:,0] = pzscale     
	# mirror, spline onto fine grid
	for n in range(len(binden)):
		intf             = interpolate.splrep(roughtheory[:,0],roughtheory[:,n+2],s=0.000000001,k=2) # skip the column with the total J for now #try interp1d with bounds_error=False and fill_value=0.0
		pzprofile[:,n+1] = interpolate.splev(abs(pzscale),intf,der=0)
	# normalize to one electron, multiply by number of electrons
	for n in range(len(binden)):
		normval = integrate.trapz(pzprofile[:,n+1],pzprofile[:,0])
		pzprofile[:,n+1] = pzprofile[:,n+1]/normval*long(occupation[n])
	binden     = [float(en) for en in binden]
	occupation = [float(val) for val in occupation]
	return pzprofile, binden, occupation

def makeprofile(element,filename=os.path.join(data_installation_dir,'data/ComptonProfiles.dat'),E0=9.69,tth=35.0,correctasym=None):
	"""
	takes the profiles from 'makepzprofile()', converts them onto eloss 
	scale and normalizes them to S(q,w) [1/eV]
	input:
	element  = element symbol (e.g. 'Si', 'Al', etc.)
	filename = path and filename to tabulated profiles
	E0       = scattering energy [keV]
	tth      = scattering angle  [deg]
	returns:
	enscale = energy loss scale
	J = total CP
	C = only core contribution to CP
	V = only valence contribution to CP
	q = momentum transfer [a.u.]
	"""
	pzprofile,binden,occ = makepzprofile(element,filename)
	# convert to eloss scale
	enscale = ((np.flipud(pz2e1(E0,pzprofile[:,0],tth))-E0)*1e3)
        q = momtrans_au(enscale/1000.0+E0,E0,tth)
	# add asymmetry if needed (2p1/2 and 2p3/2 for Z > 35 (Br))
	asymmetry = np.flipud(HRcorrect(pzprofile,occ,q));  # asymmetry flipped for conversion to e-loss scale (???)
	if correctasym:
		pzprofile[:,1:4] = pzprofile[:,1:4] + asymmetry*correctasym

	# discard profiles below zero
	hfprofile = pzprofile[np.nonzero(enscale.T>=0)[0],:]
	q         = q[np.nonzero(enscale.T>=0)[0]] #q[:,np.nonzero(enscale.T>=0)[0]]
	enscale   = enscale[np.nonzero(enscale.T>=0)[0]] #enscale[:,np.nonzero(enscale.T>=0)[0]]
	hfprofile[:,0] = enscale
	# cut at edges
	for n in range(len(binden)):
		hfprofile[np.where(enscale<binden[n]),n+1] = 0 
	# convert J(pz) to S(q,w) via J(pz)=N_electrons*hartree*q*S(q,w) and
	# normalize using the f-sum rule (sum(S(q,w)*w)=f)
	# convert to a.u.
	hartree = 1.0/constants.physical_constants['electron volt-hartree relationship'][0]
	enscaleh = enscale/hartree # eloss in a.u.
	# normalize to one then multiply by N_el*q**2/2
	for n in range(len(binden)):
		hfprofile[:,n+1] = hfprofile[:,n+1]/(integrate.trapz(np.multiply(hfprofile[:,n+1],enscaleh),enscaleh))
		hfprofile[:,n+1] = np.multiply(hfprofile[:,n+1],(q**2.0)/2.0)*occ[n]
	# convert back to [1/eV] and sum up
	# total profile J and valence V (all edges )
	J = np.zeros((len(enscale)))
	V = np.zeros((len(enscale)))
	for n in range(len(binden)):
		if binden[n] < enscale[-1]:
			J += hfprofile[:,n+1]/hartree
			if binden[n] < 10:
				V += hfprofile[:,n+1]/hartree
	C = J - V
	return enscale,J,C,V,q

def makeprofile_comp(formula,filename=os.path.join(data_installation_dir,'data/ComptonProfiles.dat'),E0=9.69,tth=35,correctasym=None):
	"""
	returns the compton profile of a chemical compound with formula 'formula'
	input:
	formula = string of a chemical formula (e.g. 'SiO2', 'Ba8Si46', etc.)
	filename = path and filename to tabulated profiles
	E0       = scattering energy [keV]
	tth      = scattering angle  [deg]
	returns:
	eloss = energy loss scale
	J = total CP
	C = only core contribution to CP
	V = only valence contribution to CP
	q = momentum transfer [a.u.]
	"""
	elements,stoichiometries = parseformula(formula)
	
	if not np.any(correctasym):
		correctasym = np.zeros(len(elements))
		
	eloss,J,C,V,q = makeprofile(elements[0],filename,E0,tth,correctasym[0])

	for n in range(len(elements[1:])):
		eloss,j,c,v,q = makeprofile(elements[n+1],filename,E0,tth,correctasym[n+1])
		J += j
		C += c
		V += v
	return eloss, J,C,V,q

def makeprofile_compds(formulas,concentrations=None,filename=os.path.join(data_installation_dir,'data/ComptonProfiles.dat'),E0=9.69,tth=35.0,correctasym=None):
	"""
	returns sum of compton profiles from a lost of chemical compounds weighted by the given concentration
	"""
	# if correctasym is not given, no HR correction is applied 
	if not np.any(concentrations):
		concentrations = np.ones(len(formulas))/len(formulas)
	if not np.any(correctasym):
		correctasym = []
		for formula in formulas:
			elements,stoichiometries = parseformula(formula)
			correctasym.append(np.zeros(len(elements)))
	
	eloss,J,C,V,q = makeprofile_comp(formulas[0],filename,E0,tth,correctasym[0])
	if len(formulas)>1:
		J = J*concentrations[0]
		C = C*concentrations[0]
		V = V*concentrations[0]
		for n in range(len(formulas[1:])):
			eloss,j,c,v,q = makeprofile_comp(formulas[n+1],filename,E0,tth,correctasym[n+1])
			J += j*concentrations[n+1]
			C += c*concentrations[n+1]
			V += v*concentrations[n+1]
	return eloss,J,C,V,q

def HRcorrect(pzprofile,occupation,q):
	""" Returns the first order correction to filled 1s, 2s, and 2p Compton profiles.

	Implementation after Holm and Ribberfors (citation ...).

	Args: 
	-----
	pzprofile (np.array): Compton profile (e.g. tabulated from Biggs) to be corrected (2D matrix). 
	occupation (list): electron configuration.
	q (float or np.array): momentum transfer in [a.u.].

	Returns:
	--------
	asymmetry (np.array):  asymmetries to be added to the raw profiles (normalized to the number of electrons on pz scale)
	"""
	# prepare output matrix
	if len(occupation) == 1:
		asymmetry = np.zeros((len(pzprofile[:,0]),1))
	elif len(occupation) == 2:
		asymmetry = np.zeros((len(pzprofile[:,0]),2))
	elif len(occupation) >= 3:
		asymmetry = np.zeros((len(pzprofile[:,0]),3))

	# take care for the cases where 2p levels have spin-orbit split taken into account in the Biggs table
	if len(occupation)>3 and occupation[2]==2 and occupation[3]==4:
		pzprofile[:,3] = pzprofile[:,3] + pzprofile[:,4]
		occupation[2] = 6
    
	# 1s 
	if occupation[0] < 2:
		pass
	else:
		# find gamma1s lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2
		fitfct  = lambda a: (np.absolute(np.max(pzprofile[:,1])-np.max(occupation[0]*8.0*a**5.0/3.0/np.pi/(a**2.0+pzprofile[:,0]**2.0)**3.0)))
		res = optimize.leastsq(fitfct,np.sum(occupation))
		gamma1s = res[0][0]
		# calculate j0 and j1
		j0 = occupation[0]*8.0*gamma1s**5.0/3.0/np.pi/((gamma1s**2.0+pzprofile[:,0]**2.0)**3.0)
		j1 = 2.0*gamma1s*np.arctan2(pzprofile[:,0],gamma1s)-3.0/2.0*pzprofile[:,0] 
		j1 = j1/q*j0
		asymmetry[:,0] = j1
	# 2s
	if len(occupation)>1:
		if occupation[1] < 2:
			pass
		else:
			# find gamma2s
			fitfct  = lambda a: (np.absolute(np.max(pzprofile[:,2])-np.max(occupation[1]*((a**4.0-10.0*a**2.0*pzprofile[:,0]**2 + 40.0*pzprofile[:,0]**4.0)*128.0*a**5.0/15.0/np.pi/(a**2.0 + 4.0*pzprofile[:,0]**2.0)**5.0))))
			res = optimize.leastsq(fitfct,np.sum(occupation)*2.0/3.0)
			gamma2s = res[0][0]
			# calculate j0 and j1
			j0 = occupation[1]*(gamma2s**4.0-10.0*gamma2s**2.0*pzprofile[:,0]**2.0+40.0*pzprofile[:,0]**4.0)*128.0*gamma2s**5.0/15.0/np.pi/(gamma2s**2.0 + 4.0*pzprofile[:,0]**2.0)**5.0
			j1 = 2.0*gamma2s*np.arctan2(2.0*pzprofile[:,0],gamma2s)-5.0/4.0*(gamma2s**4.0+48.0*pzprofile[:,0]**4.0)/(gamma2s**4.0-10.0*gamma2s**2.0*pzprofile[:,0]**2.0+40.0*pzprofile[:,0]**4.0)*pzprofile[:,0] 
			j1 = j1/q*j0
			asymmetry[:,1] = j1
	# 2p
	if len(occupation)>2:
		if occupation[2] < 6:
			pass
		else:
			forgamma = 3.0*pzprofile[:,3]/np.trapz(pzprofile[:,3],pzprofile[:,0]) # 2p correction is defined for 3 electrons in the 2p shell
			# find gamma2p
			fitfct = lambda a: (np.absolute(np.max(forgamma)-np.max(((a**2.0+20.0*pzprofile[:,0]**2.0)*64.0*a**7.0/5.0/np.pi/(a**2.0+4.0*pzprofile[:,0]**2.0)**5.0))))
			res = optimize.leastsq(fitfct,np.sum(occupation)*1.0/3.0)
			gamma2p = res[0][0]
			# calculate j0 and j1
			j0 = 2.0*(gamma2p**2.0+20.0*pzprofile[:,0]**2.0)*64.0*gamma2p**7.0/5.0/np.pi/(gamma2p**2.0+4.0*pzprofile[:,0]**2.0)**5.0
			j1 = 2.0*gamma2p*np.arctan2(2.0*pzprofile[:,0],gamma2p)-2.0/3.0*pzprofile[:,0]*(10.0*gamma2p**2.0+60.0*pzprofile[:,0]**2.0)/(gamma2p**2.0+20.0*pzprofile[:,0]**2.0)
			j1 = j1/q*j0
			asymmetry[:,2] = j1
	return asymmetry

def parseformula(formula):
	"""Parses a chemical sum formula.

	Parses the constituing elements and stoichiometries from a given 
	chemical sum formula.

	Args:
	-----
	formula (string): string of a chemical formula (e.g. 'SiO2', 'Ba8Si46', etc.)

	Returns:
	--------
	elements (list): list of strings of constituting elemental symbols.
	stoichiometries (list): list of according stoichiometries in the same order as 'elements'.
	"""
	elements = []
	stoichiometries = []
	splitted = findall(r'([A-Z][a-z]*)(\d*)',formula)
	elements.extend([element[0] for element in splitted])
	stoichiometries.extend([(int(element[1]) if element[1] else 1) for element in splitted])
	return elements,stoichiometries

def element(z):
	"""Converts atomic number into string of the element symbol and vice versa.

	Returns atomic number of given element, if z is a string of the 
	element symbol or string of element symbol of given atomic number z.

	Args:
	-----
	z (string or int): string of the element symbol or atomic number. 

	Returns:
	Z (string or int): string of the element symbol or atomic number.
	"""
	zs = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
              'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
              'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo',
              'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba',
              'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
              'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
              'At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf',
              'Es','Fm','Md','No','Lr','Ku']
	if isinstance(z,str):
		try:
			Z = zs.index(z)+1
		except:
			Z = None
			print 'Given element ' + z + ' unknown.'
	elif isinstance(z,int):
		if z > 0 and z < 105:
			Z = zs[z-1]
		else:
			print 'Element Z = '+ str(z) +' unknown.'
	else:
		print 'type '+ type(z) + 'not supported.'	
	return Z

def myprho(energy,Z,logtablefile=os.path.join(data_installation_dir,'data/logtable.dat')):
	"""Calculates the photoelectric, elastic, and inelastic absorption of 
	an element Z 

	Calculates the photelectric , elastic, and inelastic absorption of an element Z. 
	Z can be atomic number or element symbol.

	Args:
	-----
	energy (np.array): energy scale in [keV].
	Z (string or int): atomic number or string of element symbol.

	Returns:
	--------
	murho (np.array): absorption coefficient normalized by the density.
	rho (float): density in UNITS?
	m (float): atomic mass in UNITS?

	"""
	en = np.array([])
	en = np.append(en,energy) 
	logtable = np.loadtxt(logtablefile)
	# find the right places in logtable
	if not isinstance(Z,int):
		Z = element(Z)
	try:
		ind = list(logtable[:,0]).index(Z)
	except:
		print 'no such element in logtable.dat'
	c     = np.array(logtable[ind:ind+5,:]) # 5 lines that corresponds to the element
	le    = np.log(en) # logarithm of the energy
	mr    = np.exp(c[1,3]+le*(c[2,3]+le*(c[3,3]+le*c[4,3])))  # extract mu from loglog table
	i     = np.where(en<=c[0,3])
	l     = le[i]
	mr[i] = np.exp(c[1,2]+l*(c[2,2]+l*(c[3,2]+l*c[4,2])))
	i     = np.where(en<c[0,2])
	l     = le[i]
	mr[i] = np.exp(c[1,1]+l*(c[2,1]+l*(c[3,1]+l*c[4,1])))
	i     = np.where(en<c[0,1])
	l     = le[i]
	# mu
	mu      = np.zeros((len(mr),3))
	mu[:,0] = mr
	mu[i,0] = np.exp(c[1,0]+l*(c[2,0]+l*(c[3,0]+l*c[4,0]))) # photoelectric absorption
	mu[:,1] = np.exp(c[1,4]+le*(c[2,4]+le*(c[3,4]+le*c[4,4]))) # elastic absorption
	mu[:,2] = np.exp(c[1,5]+le*(c[2,5]+le*(c[3,5]+le*c[4,5]))) # inelastic abssorption
	#
	m = c[0,4] # atomic mass
	murho = mu*0.602252/m # mu/rho
	rho = c[0,5]
	return murho, rho, m

def mpr(energy,compound):
	"""Calculates the photoelectric, elastic, and inelastic absorption of 
	a chemical compound.

	Calculates the photoelectric, elastic, and inelastic absorption of a
	chemical compound.

	Args:
	-----
	energy (np.array): energy scale in [keV].
	compound (string): chemical sum formula (e.g. 'SiO2')

	Returns:
	-------- 
	murho (np.array): absorption coefficient normalized by the density.
	rho (float): density in UNITS?
	m (float): atomic mass in UNITS?

	"""
	en   = np.array([])
	en   = np.append(en,energy) # turn energy into an iterable array
	z,w  = parseformula(compound)
	mr   = np.zeros((len(en),3)) # 1. photoelectric absorption, 2. elastic absorption, 3. inelastic absorption
	rhov = np.zeros((len(z),1))
	mv   = np.zeros((len(z),1))
	for i in range(len(z)):
		tmp,rho,m = myprho(en,z[i])
		m         = m*w[i] # weigh atomic masses by stoichiometry.
		mv[i]     = m
		rhov[i]   = rho
		mr       += tmp*m # sum up individual mu/rho
	mtot = sum(mv)
	mr   = mr/mtot
	mr 	 = np.sum(mr,1)
	return mr, rhov, mv

def mpr_compds(energy,formulas,concentrations,E0,rho_formu):
	"""Calculates the photoelectric, elastic, and inelastic absorption of 
	a mix of compounds.

	Returns the photoelectric absorption for a sum of different chemical 
	compounds.

	Args:
	-----
	energy (np.array): energy scale in [keV].
	formulas (list of strings): list of chemical sum formulas

	Returns:
	--------
	murho (np.array): absorption coefficient normalized by the density.
	rho (float): density in UNITS?
	m (float): atomic mass in UNITS?

	"""
	en  = np.array([]) # turn energy into an iterable array
	en  = np.append(en,energy)
	e0  = np.array([])
	e0  = np.append(e0,E0)
	mu_tot_in  = np.zeros((len(en)))
	mu_tot_out = np.zeros((len(e0))) # should also work for series of E0's 
	for n in range(len(formulas)):
		mu_tot_in += mpr(en,formulas[n])[0]*concentrations[n]*rho_formu[n]
		mu_tot_out += mpr(e0,formulas[n])[0]*concentrations[n]*rho_formu[n]
	return mu_tot_in, mu_tot_out

def abscorr2(mu1,mu2,alpha,beta,samthick):
	"""Calculates absorption correction for given mu1 and mu2.
	Multiply the measured spectrum with this correction factor.

	This is a translation of Keijo Hamalainen's Matlab function (KH 30.05.96).

	Args:
	-----
	mu1 (np.array): absorption coefficient for the incident energy in [1/cm].
	mu2 (np.array): absorption coefficient for the scattered energy in [1/cm].
	alpha (float): incident angle relative to plane normal in [deg].
	beta (float): exit angle relative to plane normal [deg]
	              (for transmission geometry use beta < 0).
	samthick (float): sample thickness in [cm].

	Returns:
	--------
	ac (np.array): absorption correction factor. Multiply this with your measured spectrum.

	"""
	cosa = math.cos(math.radians(alpha))
	cosb = math.cos(math.radians(beta))
	if beta >= 0: # reflection geometry
		ac =  cosa*(mu1/cosa + mu2/cosb)/(1.0 - np.exp(-mu1*samthick/cosa - mu2*samthick/cosb))
	elif np.absolute(mu1/cosa - mu2/cosb).any() > np.spacing(1): # transmission geometry
		ac = -cosa*(mu1/cosa - mu2/cosb)/(np.exp(-mu1*samthick/cosa) - np.exp(-mu2*samthick/cosb))
	else:
		ac = cosa/(samthick*np.exp(-mu1*samthick/cosa))
	return ac

def absCorrection(mu1,mu2,alpha,beta,samthick,geometry='transmission'):
	"""
        **absCorrection**

        Calculates absorption correction for given mu1 and mu2.
	Multiply the measured spectrum with this correction factor.
        This is a translation of Keijo Hamalainen's Matlab function (KH 30.05.96).

	Args
	----
	mu1 : np.array
            Absorption coefficient for the incident energy in [1/cm].
	mu2 : np.array
            Absorption coefficient for the scattered energy in [1/cm].
	alpha : float
            Incident angle relative to plane normal in [deg].
	beta : float
            Exit angle relative to plane normal [deg].
	samthick : float 
            Sample thickness in [cm].
        geometry : string, optional
            Key word for different sample geometries ('transmission', 'reflection', 'sphere'). 
            If *geometry* is set to 'sphere', no angular dependence is assumed.

	Returns
	-------
	ac : np.array
            Absorption correction factor. Multiply this with your measured spectrum.

	"""
	cosa = np.cos(math.radians(alpha))
	cosb = np.cos(math.radians(beta))

        # reflection geometry
	if geometry == 'reflection':
                if beta >= 90.0:
                        print('WARNING: are you sure about the beta angle?')
		ac =  cosa*(mu1/cosa + mu2/cosb)/(1.0 - np.exp(-mu1*samthick/cosa - mu2*samthick/cosb))

        # transmission geometry
        elif geometry == 'transmission' and np.absolute(mu1/cosa - mu2/cosb).any() > np.spacing(1):
		ac = -cosa*(mu1/cosa - mu2/cosb)/(np.exp(-mu1*samthick/cosa) - np.exp(-mu2*samthick/cosb))
	elif geometry == 'transmission' and np.absolute(mu1/cosa - mu2/cosb).any() <= np.spacing(1):
		ac = cosa/(samthick*np.exp(-mu1*samthick/cosa))

        # spherical sample
        elif geometry == 'sphere':
                ac = (mu1 + mu2)/(1.0 - np.exp(-mu1*samthick -mu2*samthick))

	return ac

def gettransmission(energy,formulas,concentrations,densities,thickness):
	"""
	returns the transmission through a sample composed of chemical formulas 
	with certain densities mixed to certain concentrations, and a thickness
	"""
	en  = np.array([]) # turn energy into an iterable array
	en  = np.append(en,energy)
	if not isinstance(formulas,list):
		theformulas = []
		theformulas.append(formulas)
	else:
		theformulas = formulas
	if not isinstance(concentrations,list):
		theconcentrations = []
		theconcentrations.append(concentrations)
	else:
		theconcentrations = concentrations
	if not isinstance(densities,list):
		thedensities = []
		thedensities.append(densities)
	else:
		thedensities = densities
	# get mu
	mu_tot = np.zeros((len(en)))
	for n in range(len(theformulas)):
		 mu_tot += mpr(en,theformulas[n])[0]*theconcentrations[n]*thedensities[n]
	return np.exp(-mu_tot*thickness)

def plottransmission(energy,formulas,concentrations,densities,thickness):
	"""
	opens a plot with the transmission plotted along the given energy vector
	"""
	if not isinstance(formulas,list):
		theformulas = []
		theformulas.append(formulas)
	else:
		theformulas = formulas
	if not isinstance(concentrations,list):
		theconcentrations = []
		theconcentrations.append(concentrations)
	else:
		theconcentrations = concentrations
	if not isinstance(densities,list):
		thedensities = []
		thedensities.append(densities)
	else:
		thedensities = densities
	transmission = gettransmission(energy,formulas,concentrations,densities,thickness)
	plt.plot(energy,transmission)
	titlestring = 'transmission of: ' + ' '.join(formulas)
	plt.title(titlestring)
	plt.xlabel('energy [keV]')
	plt.ylabel('transmission [%]')
	plt.grid(False)
	plt.show()

def getpenetrationdepth(energy,formulas,concentrations,densities):
	"""
	returns the penetration depth of a mixture of chemical formulas
	with certain concentrations and densities
	"""
	en  = np.array([]) # turn energy into an iterable array
	en  = np.append(en,energy)
	if not isinstance(formulas,list):
		theformulas = []
		theformulas.append(formulas)
	else:
		theformulas = formulas
	if not isinstance(concentrations,list):
		theconcentrations = []
		theconcentrations.append(concentrations)
	else:
		theconcentrations = concentrations
	if not isinstance(densities,list):
		thedensities = []
		thedensities.append(densities)
	else:
		thedensities = densities
	# get mu
	mu_tot = np.zeros((len(en)))
	for n in range(len(theformulas)):
		 mu_tot += mpr(en,theformulas[n])[0]*theconcentrations[n]*thedensities[n]
	return 1.0/mu_tot

def plotpenetrationdepth(energy,formulas,concentrations,densities):
	"""
	opens a plot window of the penetration depth of a mixture of chemical formulas
	with certain concentrations and densities plotted along the given energy vector
	"""
	if not isinstance(formulas,list):
		theformulas = []
		theformulas.append(formulas)
	else:
		theformulas = formulas
	if not isinstance(concentrations,list):
		theconcentrations = []
		theconcentrations.append(concentrations)
	else:
		theconcentrations = concentrations
	if not isinstance(densities,list):
		thedensities = []
		thedensities.append(densities)
	else:
		thedensities = densities
	pendepth = getpenetrationdepth(energy,formulas,concentrations,densities)
	plt.plot(energy,pendepth)
	titlestring = 'penetration depth of: ' + ' '.join(formulas)
	plt.title(titlestring)
	plt.xlabel('energy [keV]')
	plt.ylabel('penetration depth [cm]')
	plt.grid(False)
	plt.show()

def sumx(A):
	"""
	Short-hand command to sum over 1st dimension of a N-D matrix (N>2) and to squeeze it to N-1-D matrix.
	"""
	return np.squeeze(np.sum(A,axis=0))

def specread(filename,nscan):
	"""
	reads scan "nscan" from SPEC-file "filename"
	INPUT: 	filename = string with the SPEC-file name
		   	nscan    = number (int) of desired scan 
	OUTPUT: data     =
			motors   =
			counters = dictionary
	"""
	scannid   = '#S'
	countid   = '#L'
	motorid   = '#P'
	data      = []
	motors    = []
	counterss = []    
	f = open(filename,'r')
	while True:
		line = f.readline()
		if not line: break
		if line[0:2] == scannid:
			if int(line.split()[1]) == nscan:
				line = '##'+line
				while line and line[0:2]!='#S':
					line = f.readline() 
					if not line:
						break                    
					if line[0:2] == countid:
						cline = '  '+line[2:]
						counterss = [n.strip() for n in filter(None,cline.split('  ')[1:])]
					if line[0:2] == motorid:
						motors.append([float(n) for n in line.strip().split()[1:]])                    
					if line[0] != '#':
						data.append([float(n) for n in line.strip().split()])
	data.pop(-1) # the trailing empty line                    
	f.close()
	# put the data into a dictionary with entries from the counterss
	counters = {}
	for n in range(len(counterss)):
		counters[counterss[n].lower()] = [row[n] for row in data] # data[:,n]
	return data, motors, counters

def edfread(filename):
	"""
	reads edf-file with filename "filename"
	OUTPUT:	data = 256x256 numpy array
	"""	
	# get some info from header
	f = open(filename,'rb').readlines()
	counter = 0
	predata = []
	for entry in f:
		counter += 1
		if entry.strip().split()[0] == '}':
			break
	for entry in f[:counter]:
		if entry.strip().split()[0] == 'Dim_1':
			dim1 = int(entry.strip().split()[2])
		if entry.strip().split()[0] == 'Dim_2':
			dim2 = int(entry.strip().split()[2])
		if entry.strip().split()[0] == 'Size':
			size = int(entry.strip().split()[2])
		if entry.strip().split()[0] == 'UnsignedShort':
			type_code = 'H'
		if entry.strip().split()[0] == 'SignedInteger':
			type_code = 'i'
	length = 0
	for line in f:
		length += len(line)
	headerlength = (length-size)/2			
	# get the data
	f = open(filename,'rb')
	predata = arr.array(type_code)
	predata.fromfile(f,(headerlength+dim1*dim2)) # this prevents the header (1024 characters long) to end up in the 256x256 picture
	data = np.reshape(predata[headerlength:],(dim1,dim2)) # this prevents the header (1024 characters long) to end up in the 256x256 picture
	f.close()
	return data

def edfread_test(filename):
	"""
	reads edf-file with filename "filename"
	OUTPUT:	data = 256x256 numpy array

	here is how i opened the HH data: 
	data = np.fromfile(f,np.int32)
	image = np.reshape(data,(dim,dim))

	"""	
	# get some info from header
	f = open(filename,'rb').readlines()
	counter = 0
	predata = []
	for entry in f:
		counter += 1
		if entry.strip().split()[0] == '}':
			break
	for entry in f[:counter]:
		if entry.strip().split()[0] == 'Dim_1':
			dim1 = int(entry.strip().split()[2])
		if entry.strip().split()[0] == 'Dim_2':
			dim2 = int(entry.strip().split()[2])
		if entry.strip().split()[0] == 'Size':
			size = int(entry.strip().split()[2])
	length = 0
	for line in f:
		length += len(line)
	headerlength = (length-size)/2			
	# get the data
	f = open(filename,'rb')
	predata = arr.array('H')
	predata.fromfile(f,(headerlength+dim1*dim2)) # this prevents the header (1024 characters long) to end up in the 256x256 picture
	data = np.reshape(predata[headerlength:],(dim2,dim1)) # this prevents the header (1024 characters long) to end up in the 256x256 picture
	f.close()
	return data

def momtrans_au(e1,e2,tth):
	"""
	Calculates the momentum transfer in atomic units
	input: 
	e1  = incident energy  [keV]	
	e2  = scattered energy [keV]
	tth = scattering angle [deg]
	returns:
	q   = momentum transfer [a.u.] (corresponding to sin(th)/lambda)
	"""
	e1    = np.array(e1*1e3/13.60569172/2)
	e2    = np.array(e2*1e3/13.60569172/2)
	th    = np.radians(tth)#tth/180.0*numpy.pi
	hbarc = 137.03599976
	q     = 1/hbarc*np.sqrt(e1**2.0+e2**2.0-2.0*e1*e2*np.cos(th));
	return q

def momtrans_inva(e1,e2,tth):
	"""
	Calculates the momentum transfer in inverse angstrom
	input: 
	e1  = incident energy  [keV]	
	e2  = scattered energy [keV]
	tth = scattering angle [deg]
	returns:
	q   = momentum transfer [a.u.] (corresponding to sin(th)/lambda)
	"""
	e = 1.602e-19
	c = 2.9979e8 
	hbar = 6.626e-34/2/np.pi

	e1    = np.array(e1*1e3*e/c/hbar)
	e2    = np.array(e2*1e3*e/c/hbar)
	th    = np.radians(tth)
	q     = np.sqrt(e1**2+e2**2-2*e1*e2*np.cos(th))/1e10
	return q

def energy_monoangle(angle,d=5.4307/np.sqrt(11)):
	"""
	% ENERGY  Calculates energy corrresponing to Bragg angle for given d-spacing
	%         function e=energy(dspace,bragg_angle)
	%
	%         dspace for reflection (defaulf for Si(311) reflection)
	%         bragg_angle in DEG
	%
	%         KH 28.09.93
	%
	"""
	hc = 12.3984191 # CODATA 2002 physics.nist.gov/constants
	e  = (2.0*d*sin(angle/180.0*np.pi)/hc)**(-1.0)
	return e

def find_center_of_mass(x,y):
	"""
	Returns the center of mass (first moment) for the given curve y(x)
	"""
	deno = np.trapz(y,x)
	if deno==0.0:
		return 0.0
		# print "*** print_tb:"
		# traceback.print_stack()
		# print " DENO==0!"
		# return 0.0
	return np.trapz(y*x,x)/deno


def odefctn(y,t,abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1):
	"""
	#%    [T,Y] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2,...) passes the additional
	#%    parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
	#%    all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
	#%    no options are set.   
	"""
	#print 'shape of y is ' , np.shape(y), np.shape(t)
	fcomp = 1.0/(complex(0,-lex)) * (-2.0*((abb0*(abb8 + abb7*sgbeta*t) + abb1) + complex(0,y0))*(y[0] + complex(0,y[1])) + c1*(1.0 + (y[0] + complex(0,y[1]))**2.0))
	return fcomp.real,fcomp.imag

def taupgen(e, hkl = [6,6,0], crystals = 'Si', R = 1.0, dev = np.arange(-50.0,150.0,1.0), alpha = 0.0):
	"""
	% TAUPGEN          Calculates the reflectivity curves of bent crystals
	%
	% function [refl,e,dev]=taupgen_new(e,hkl,crystals,R,dev,alpha);
	%
	%              e = fixed nominal energy in keV
	%            hkl = reflection order vector, e.g. [1 1 1]
	%       crystals = crystal string, e.g. 'si' or 'ge'
	%              R = bending radius in meters
	%            dev = deviation parameter for which the 
	%                  curve will be calculated (vector) (optional)
	%          alpha = asymmetry angle 
	% based on a FORTRAN program of Michael Krisch
	% Translitterated to Matlab by Simo Huotari 2006, 2007
	% Is far away from being good matlab writing - mostly copy&paste from
	% the fortran routines. Frankly, my dear, I don't give a damn. 
	% Complaints -> /dev/null
	"""
	prefix = data_installation_dir+'/'
	path = prefix + 'data/chitables/chitable_' # path to chitables
	# load the according chitable (tabulated)
	hkl_string = str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2]))
	filestring = path + crystals.lower() + hkl_string + '.dat'
	chi = np.loadtxt(filestring)

	# good for 1 m bent crystals in backscattering
	ystart = -50.0 # start value of angular range in arcsecs
	yend   = 150.0 # end value of angular range in arcsecs
	ystep  = 1.0   # step width in arcsecs

	if len(chi[:,0]) == 1:
		print ' I will only  calculate for the following energy: ' + '%.4f' % chi[0,0] + ' keV!!!'
	else:
		if e < np.min(chi[:,0]) or e > np.max(chi[:,0]):
			print 'Energy outside of the range in ' + filestring
			return

		chi0r = np.interp(e,chi[:,0],chi[:,1])
		chi0i = np.interp(e,chi[:,0],chi[:,2])
		chihr = np.interp(e,chi[:,0],chi[:,3])
		chihi = np.interp(e,chi[:,0],chi[:,4])

	th = braggd(hkl,e,crystals)
	lam = 12.3984191/e/10.0 # wavelength in nm

	reflcorr = 0.0
	chi0 = complex(chi0r,chi0i)
	chih = complex(chihr,chihi)

	if crystals.upper() == 'SI':
		s13 = -0.278
	elif crystals.upper() == 'GE':
		s13 = -0.273
	else:
		print 'Poisson ratio for this crystal not defined'
		return

	s15 = -0.0 # s15/s11
	dsp = dspace(hkl,crystals)/10.0 # dspace

	dwf    = 1.0 # dwf = 0.899577 # debye-waller factor
	radius = R # meridional bending radius
	rsag   = R*np.sin(np.radians(th))**2.0 # sagittal bending radius
	thick  = 500.0 # thickness in micrometers #rsag = R

	lam      = lam*1e-9
	dsp      = dsp*1e-9
	alpha    = np.radians(alpha) # alpha in rad
	thick    = thick*1e-6
	ystart   = ystart/3600.0/180.0*np.pi
	yend     = yend/3600.0/180.0*np.pi
	ystep    = ystep/3600.0/180*np.pi
	dev      = dev/3600.0/180.0*np.pi
	reflcorr = reflcorr/3600.0/180.0*np.pi

	thetab = np.arcsin(lam/(2.0*dsp))
	cpol   = 1.0 # cpol=0.5*(1+cos(2*thetab).^2) # cpol=cos(2*thetab).^2

	# gamma0 = sin(thetab+alpha) # normal convention
	# gammah = -sin(thetab-alpha) # normal convention
	gammah = -np.sin(thetab + alpha) # Krisch et al. convention (really!)
	gamma0 = np.sin(thetab - alpha) # Krisch et al. convention (I'm not kidding!!)

	beta  = gamma0/np.abs(gammah)
	gamma = gammah/gamma0

	a0 = np.sqrt(1-gamma0**2.0)
	ah = np.sqrt(1-gammah**2.0)

	mu = -2.0*np.pi/lam*chi0i

	tdepth = 1.0/mu/(1.0/np.abs(gamma0)+1.0/np.abs(gammah))

	lex = lam*np.sqrt(gamma0*np.abs(gammah))/(np.pi*chihr)

	y0 = chi0i*(1.0+beta)/(2.0*np.sqrt(beta)*chihr)

	pfried = -chihi/chihr

	c1 = cpol*dwf* complex(1.0,pfried)

	#abbreviation concerning the deviation parameter y
	abb0 = -np.sqrt(beta)/2.0/chihr
	abb1 = chi0r*(1.0+beta)/(2.0*np.sqrt(beta)*chihr)

	#abbreviations concerning the deformation field

	abb2 = gamma0*gammah*(gamma0-gammah)
	abb3 = 1.0 + 1.0/(gamma0*gammah)
	abb4 = s13*(1.0 + radius/rsag)
	abb5 = (ah - a0)/(gamma0 - gammah)*s15 
	abb6 = 1.0/(np.abs(cpol)*chihr*np.cos(thetab)*radius)
	abb7 = 2.0*np.abs(cpol)*chihr*np.cos(thetab)/gamma0

	#   a spectrometer based on a spherical diced analyzer crystal with a 1-m bending radius in nearly backscattering conditions utilizing a strain gradient beta
	sgbeta = abb6*(abb2*(abb3 - abb4 + abb5))

	nstep=len(dev)
	eta  = np.zeros_like(dev)
	abb8z = np.zeros_like(dev)
	refl  = np.zeros_like(dev)
	refl1 = np.zeros_like(dev)
	refl2 = np.zeros_like(dev)
	for l in range(nstep):
		# actual value of the deviation angle
		# dev[l] = ystart + (l - 1)*ystep

		# deviation parameter
		abb8   = -2.0*np.sin(2.0*thetab)*dev[l]
		eta[l] = (dev[l]*np.sin(2.0*thetab)+np.abs(chi0.real)/2.0*(1.0-gamma))/(np.abs(cpol)*np.sqrt(np.abs(gamma))*np.sqrt(chih*chih))
		eta[l] = eta[l].real

		ndiff = 2
		xend = 0
		x = np.max([-10.0*tdepth, -thick])
		y = np.array([0.0, 0.0])
		h = xend
		abb8z[l] = abb8

		# in this point call the subroutine
		#     [T,Y] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2,...) passes the additional
		#    parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
		#    all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
		#    no options are set.   
		#print 'the fucking shape of y is ', np.shape(y)
		T = np.arange(x,xend,1e-8)
		Y = odeint(odefctn,y,T,args=(abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1)) 

		# normalized reflectivity at this point
		refl[l] = np.sum(Y[-1,:]**2.0)
		refl1[l] = Y[-1,0]
		refl2[l] = Y[-1,1]

	de = dev * e * 1.0e6 /np.tan(thetab)

	lam    = lam *1.0e+09        
	dsp    = dsp*1.0e+09        
	alpha  = alpha/np.pi*180.0        
	ystart = ystart/4.848136811e-06
	yend   = yend/4.848136811e-06   
	ystep  = ystep/4.848136811e-06
	dev    = dev/4.848136811e-06 # dev in arcsecs
	
	dev = dev/3600.0 # in degrees
	thb = th
	th  = thb + dev
	e0  = e
	e   = energy(dspace(hkl,crystals),th)-e0
	e = e*1e6

	dev = dev*3600.0 # back to arcsecs

	return refl,e,dev,e0

def readfio(prefix, scannumber, repnumber=0):
	"""
	if repnumber = 0:
	reads a spectra-file (name: prefix_scannumber.fio)
	if repnumber > 1:
	reads a spectra-file (name: prefix_scannumber_rrepnumber.fio)
	"""
	suffix = '.fio'
	filename = prefix + '%05d' % scannumber + suffix
	if repnumber > 0:
		filename = prefix + '%05d' % scannumber + 'r' + '%d' % repnumber + suffix

	# analyze structure of file	
	fid     = open(filename,'r')

	colnameflag  = ' Col'
	colstartflag = '%d'	
	colnames = []

	linenum  = 0
	for line in fid:
		linenum +=1
		if colnameflag in line: 
			colnames.append(line.strip())
		if colnameflag in line: 
			startline = linenum
	fid.close()
	thefile = open(filename,'r').readlines()
	data = []
	for line in thefile[(len(colnames)+startline):]:
		data.append([float(x) for x in line.strip().split()])

	return np.array(data), colnames

def energy_monoangle(angle,d=5.4307/np.sqrt(11)):
	"""
	% ENERGY  Calculates energy corrresponing to Bragg angle for given d-spacing
	%         function e=energy(dspace,bragg_angle)
	%
	%         dspace for reflection (defaulf for Si(311) reflection)
	%         bragg_angle in DEG
	%
	%         KH 28.09.93
	%
	"""
	hc = 12.3984191 # CODATA 2002 physics.nist.gov/constants
	e  = (2.0*d*sin(angle/180.0*np.pi)/hc)**(-1.0)
	return e

def readp01image(filename):
	"""
	reads a detector file from PetraIII beamline P01
	"""
	dim = 256
	f = open(filename,'rb')
	data = np.fromfile(f,np.int32)
	#	predata = arr.array('H')
	#	predata.fromfile(f,(dim*dim))
	image = np.reshape(data,(dim,dim))
	f.close()
	return image

def readp01scan(prefix,scannumber):
	"""
	reads a whole scan from PetraIII beamline P01 (experimental)
	"""
	print ("parsing files of scan No. %s" % scannumber)
	#fioname = prefix + 'online/hasylab_' + "%05d" % scannumber + '.fio'
	fioprefix = prefix + 'online/ixs_scan_'
	fiodata = readfio(fioprefix,scannumber)[0]

	mats1 = np.zeros((np.shape(fiodata)[0],256,256))
	mats2 = np.zeros((np.shape(fiodata)[0],256,256))
	mats  = np.zeros((np.shape(fiodata)[0],256,256*2))

	for n in range(np.shape(fiodata)[0]):
		matname1 = prefix + 'ixs_scan_' + "%05d" % scannumber + '/mdpxa/ixs_scan_' + "%05d" % scannumber + '_a_' + "%05d" % (n+1)
		matname2 = prefix + 'ixs_scan_' + "%05d" % scannumber + '/mdpxa/ixs_scan_' + "%05d" % scannumber + '_b_' + "%05d" % (n+1)

		mats1[n,:,:] = readp01image(matname1)
		mats2[n,:,:] = readp01image(matname2)
		mats[n,:,0:256] = mats1[n,:,:]
		mats[n,:,256:]  = mats2[n,:,:]
	return fiodata, mats, mats1, mats2

def readp01scan_rep(prefix,scannumber,repetition):
	"""
	reads a whole scan with repititions from PetraIII beamline P01 (experimental)
	"""
	print ("parsing files of scan No. %s" % scannumber)
	#fioname = prefix + 'online/hasylab_' + "%05d" % scannumber + 'r' + "%1d" % repetition + '.fio'
	fioprefix = prefix + 'online/ixs_scan_'
	fiodata = readfio(fioprefix,scannumber,repetition)[0]

	mats1 = np.zeros((np.shape(fiodata)[0],256,256))
	mats2 = np.zeros((np.shape(fiodata)[0],256,256))
	mats  = np.zeros((np.shape(fiodata)[0],256,256*2))

	for n in range(np.shape(fiodata)[0]):
		matname1 = prefix + 'ixs_scan_' + "%05d" % scannumber + 'r' + "%1d" % repetition  + '/mdpxa/ixs_scan_' + "%05d" % scannumber + 'r' + "%1d" % repetition + '_a_' + "%05d" % (n+1)
		matname2 = prefix + 'ixs_scan_' + "%05d" % scannumber + 'r' + "%1d" % repetition  + '/mdpxa/ixs_scan_' + "%05d" % scannumber + 'r' + "%1d" % repetition + '_b_' + "%05d" % (n+1)

		mats1[n,:,:] = readp01image(matname1)
		mats2[n,:,:] = readp01image(matname2)
		mats[n,:,0:256] = mats1[n,:,:]
		mats[n,:,256:]  = mats2[n,:,:]
	return fiodata, mats, mats1, mats2


