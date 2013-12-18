#!/usr/bin/python
# Filename: extraction.py

import os
import math

import numpy as np
import array as arr
import matplotlib.pyplot as plt

from itertools import groupby
from scipy.integrate import trapz
from scipy import interpolate, signal, integrate, constants, optimize
from re import findall
from scipy.ndimage import measurements
from scipy.optimize import leastsq, fmin
from scipy.interpolate import Rbf, RectBivariateSpline

def fwhm(x,y):
	"""
	finds full width at half maximum of the curve y vs. x
	returns 
	f  = FWHM
	x0 = position of the maximum
	"""
	if x[-1] < x[0]:
		x = np.flipud(x)
		y = flipud(y)

	y0 = np.amax(y)
	i0 = np.where(y == y0)[0][0] # this is because it returns several values sometimes (always take the first), needs improvement
	x0 = x[i0]

	i1 = np.where(np.logical_and(y>y/3.0, x<x0))[0]
	i2 = np.where(np.logical_and(y>y/3.0, x>x0))[0]

	f  = interpolate.interp1d(y[i1],x[i1], bounds_error=False, fill_value=0.0)
	x1 = f(y0/2.0)
	f  = interpolate.interp1d(y[i2],x[i2], bounds_error=False, fill_value=0.0)
	x2 = f(y0/2.0)
	fwhm = x2 - x1
	x0 = np.mean([x2, x1])
	return 2.0*fwhm, x0

def gauss(x,x0,fwhm):
	"""
	area-normalized gaussian
	x    = x-axis vector
	x0   = position 
	fwhm = full width at half maximum
	"""
	sigma = fwhm/(2*np.sqrt(2*np.log(2)));
	y = np.exp(-(x-x0)**2/2/sigma**2)/sigma/np.sqrt(2*np.pi)
	return y

def gauss_areanorm(x,x0,fwhm):
	"""
	area-normalized gaussian
	x    = x-axis vector
	x0   = position 
	fwhm = full width at half maximum
	"""
	sigma = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
	y = np.exp(-(x-x0)**2.0/2.0/sigma**2)/sigma/np.sqrt(2.0*np.pi)
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
	Extrapolates the smaller and larger values as a constant
	x  = old x-axis
	y  = old y-axis
	x2 = new x-axis
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

def readxas(filename):
	"""
	function output = readxas(filename)%[e,p,s,px,py,pz] = readxas(filename)

	% READSTF   Load StoBe fort.11 (XAS output) data
	%
	%   [E,P,S,PX,PY,PZ] = READXAS(FILENAME)
	%
	%      E        energy transfer [eV]
	%      P        dipole transition intensity
	%      S        r^2 transition intensity 
	%      PX       dipole transition intensity along x
	%      PY       dipole transition intensity along y
	%      PZ       dipole transition intensity along z
	%
	%                   as line diagrams.
	%
	%                             T Pylkkanen @ 2011-10-17 
	"""
	# Open file
	f = open(filename,'r')
	lines = f.readlines()
	f.close()
	data = []
	for line in lines[1:]:
		data.append([float(x) for x in line.replace('D', 'e').strip().split()])

	data      = np.array(data)
	data[:,0] = data[:,0]*27.211384565719481 # convert from a.u. to eV
	

	data[:,3] = 2.0/3.0*data[:,3]**2.0*e/27.211384565719481 # osc(x)
	data[:,4] = 2.0/3.0*data[:,4]**2.0*e/27.211384565719481 # osc(y)
	data[:,5] = 2.0/3.0*data[:,5]**2.0*e/27.211384565719481 # osc(z)

	return data

def broaden_diagram(e,s,params=[1.0, 1.0, 537.5, 540.0],npoints=1000):
	"""
	function [e2,s2] = broaden_diagram2(e,s,params,npoints)

	% BROADEN_DIAGRAM2   Broaden a StoBe line diagram
	%
	%  [ENE2,SQW2] = BROADEN_DIAGRAM2(ENE,SQW,PARAMS,NPOINTS)
	%
	%   gives the broadened spectrum SQW2(ENE2) of the line-spectrum
	%   SWQ(ENE). Each line is substituted with a Gaussian peak,
	%   the FWHM of which is determined by PARAMS. ENE2 is a linear
	%   scale of length NPOINTS (default 1000).
	%
	%    PARAMS = [f_min f_max emin max]
	%
	%     For ENE <= e_min, FWHM = f_min.
	%     For ENE >= e_max, FWHM = f_min.
	%     FWHM increases linearly from [f_min f_max] between [e_min e_max].
	%
	%                               T Pylkkanen @ 2008-04-18 [17:37]
	"""
	f_min = params[0]
	f_max = params[1]
	e_min = params[2]
	e_max = params[3]

	e2   = np.linspace(np.min(e)-10.0,np.max(e)+10.0,npoints);
	s2   = np.zeros_like(e2) 
	fwhm = np.zeros_like(e)

	# FWHM: Constant  -- Linear -- Constant
	A    = (f_max-f_min)/(e_max-e_min)
	B    = f_min - A*e_min
	fwhm = A*e + B
	inds = e <= e_min
	fwhm[inds] = f_min
	inds = e >= e_max	
	fwhm[inds] = f_max

	for i in range(len(s)):
		s2 += s[i]*gauss(e2,e[i],fwhm[i])

	return e2, s2

def broaden_linear(spec,params=[0.8, 8, 537.5, 550],npoints=1000):
	"""
	broadens a spectrum with a Gaussian of width params[0] below 
	params[2] and width params[1] above params[3], width increases 
	linear in between.
	returns two-column numpy array of length npoints with energy and the broadened spectrum
	"""
	evals = spec[:,0]
	sticks= spec[:,1]
	f_min = params[0]
	f_max = params[1]
	e_min = params[2]
	e_max = params[3]
	e2    = np.linspace(np.min(evals)-10.0,np.max(evals)+10.0,npoints)
	s2    = np.zeros(len(e2))
	fwhm  = np.zeros(len(evals))
	# FWHM: Constant  -- Linear -- Constant
	A    = (f_max-f_min)/(e_max-e_min)
	B    = f_min - A*e_min
	fwhm = A*evals + B
	fwhm[evals <= e_min] = f_min
	fwhm[evals >= e_max] = f_max
	for n in range(len(sticks)):
		s2 = s2 + sticks[n]*gauss(e2,evals[n],fwhm[n])
	spectrum = np.zeros((len(e2),2))
	spectrum[:,0] = e2
	spectrum[:,1] = s2
	return spectrum

def load_stobe_specs(prefix,postfix,fromnumber,tonumber,step,stepformat=2):
	"""
	load a bunch of StoBe calculations, which filenames are made up of the 
	prefix, postfix, and the counter in the between the prefix and postfix 
	runs from 'fromnumber' to 'tonumber' in steps of 'step' (number of digits 
	is 'stepformat')
	"""
	numbers = np.linspace(fromnumber,tonumber,(tonumber-fromnumber + step)/step)
	filenames = []
	precision = '%0'+str(stepformat)+'d'
	for number in numbers:
		thenumber = precision % number
		thefilename = prefix+thenumber+postfix
		filenames.append(thefilename)
	specs = []
	for filename in filenames:
		try:
			specs.append(readxas(filename))
		except:
			print 'found no file: ' + filename
	return specs

def load_erkale_spec(filename):
	"""
	returns an erkale spectrum
	"""
	spec = np.loadtxt(filename)
	return spec

def load_erkale_specs(prefix,postfix,fromnumber,tonumber,step,stepformat=2):
	"""
	returns a list of erkale spectra
	"""
	numbers = np.linspace(fromnumber,tonumber,(tonumber-fromnumber + step)/step)
	filenames = []
	precision = '%0'+str(stepformat)+'d'
	for number in numbers:
		thenumber = precision % number
		thefilename = prefix+thenumber+postfix
		filenames.append(thefilename)
	specs = []
	for filename in filenames:
		try:
			specs.append(load_erkale_spec(filename))
		except:
			print 'found no file: ' + filename
	return specs

def cut_spec(spec,emin=None,emax=None):
	"""
	deletes lines of matrix with first column smaller than emin and larger than emax 
	"""
	if not emin:
		emin = spec[0,0]
	if not emax:
		emax = spec[-1,0]
	spec = spec[spec[:,0]>emin]
	spec = spec[spec[:,0]<emax]
	return spec

# extraction
def constant(x,a):
	"""
	returns a constant
	x = x-axis
	a = value of the constant 
	returns vector with same length as x with entries a
	"""
	x = np.array(x)
	y = np.zeros(len(x)) + a
	return y

def linear(x,a):
	"""
	returns a linear function y = ax + b
	x = x-axis
	a = list of parameters, a[0] = slope, a[1] = y-offset
	it is better to use numpy's polyval
	"""
	x = np.array(x)
	y = a[0]*x + a[1]
	return y

def lorentz(x,a):
	"""
	% LORENTZ  The Lorentzian function
	%
	%          y  = lorentz(x,a)
	%          x      = x-vector
	%          a[0]   = peak position
	%          a[1]   = FWHM
	%          a[2]   = Imax (if not defined, Imax = 1)
	"""
	y = a[2]*(((x-a[0])/(a[1]/2.0))**2.0+1.0)**(-1.0)
	return y

def pearson7(x,a):
	"""
	returns a pearson function
	a[0] = Peak position
	a[1] = FWHM
	a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
	a[3] = Peak intensity
	a[4] = Background
	"""
	x = np.array(x)
	y = a[3] * (1.0+(2.0**(1.0/a[2])-1.0) * (2.0*(x-a[0])/a[1])**2.0)**(-a[2])+a[4]
	return y

def pearson7_forcurvefit(x,a0,a1,a2,a3,a4):
	"""
	special version of person7 for use in constrained/bound minimizations
	returns a pearson function
	a[0] = Peak position
	a[1] = FWHM
	a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
	a[3] = Peak intensity
	a[4] = Background
	"""
	x = np.array(x)
	y = a3 * (1.0+(2.0**(1.0/a2)-1.0) * (2.0*(x-a0)/a1)**2.0)**(-a2)+a4
	return y

def pearson7_zeroback(x,a):
	"""
	returns a pearson function but without y-offset
	a[0] = Peak position
	a[1] = FWHM
	a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
	a[3] = Peak intensity
	"""
	x = np.array(x)
	y = a[3] * (1.0+(2.0**(1.0/a[2])-1.0) * (2.0*(x-a[0])/a[1])**2.0)**(-a[2])
	return y

def gauss(x,a):
	"""
	returns a gaussian with peak value normalized to unity
	a[0] = peak position
	a[1] = Full Width at Half Maximum
	"""
	y = np.exp(-np.log(2.0)*((x-a[0])/a[1]*2.0)**2.0)
	return y

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

def e2pz(w1,w2,th):
	"""
	Calculates the momentum scale and the relativistic Compton cross section 
	correction according to P. Holm, PRA 37, 3706 (1988).
	from KH 29.05.96
	input:
	w1 = incident energy  [keV]
	w2 = scattered energy [keV]
	th = scattering angle [deg]
	returns:
	pz = momentum scale   [a.u.]
	cf = cross section correction factor such that:
	J(pz) = cf * d^2(sigma)/d(w2)*d(Omega) [barn/atom/keV/srad]
	"""
	w1  = np.array(w1)    # make sure arrays are used
	w2  = np.array(w2)           
	m   = constants.value('electron mass energy equivalent in MeV')*1e3 #511.003      # Mass in natural units
	th  = math.radians(th) # th/180.0*np.pi  # Angle to radians
	alp = constants.value('fine-structure constant') #1.0/137.036  # Fine structure constant
	r0  = constants.value('classical electron radius') #2.8179e-15   # Electron radius
	q   = np.sqrt(w1**2 + w2**2-2*w1*w2*np.cos(th))                        # Momentum transfer    
	pz  = q/2.0 - (w1-w2) * np.sqrt(1/4.0 + m**2/(2*w1*w2*(1-np.cos(th)))) # In natural units
	E   = np.sqrt(m**2+pz**2)
	A   = ((w1-w2)*E-w1*w2*(1-np.cos(th)))/q
	D   = (w1-w2*np.cos(th))*A/q
	R   = w1*(E-D)
	R2  = R-w1*w2*(1-np.cos(th))
	chi = R/R2 + R2/R + 2.0*m**2 * (1/R-1/R2) + m**4 * (1/R-1/R2)**2
	cf  = 2.0*w1*q*E/(m**2*r0**2*w2*chi)
	cf  = cf*(1e-28*(m*alp))			                             # Cross section now in barns/atom/keV/srad
	pz  = pz/(m*alp)                                                 # pz to atomic units (a.u.)
	return pz, cf

def pz2e1(w2,pz,th):
	"""
	Calculates the incident energy value corresponding
	specific scattered photon and momentum value.
	from KH 29.05.96
	input: 
	w2 = scattered energy [keV]
	pz = momentum value   [a.u.]
	th = scattering angle [deg]
	returns:
	w1 = incident energy  [keV]
	"""
	pz  = np.array(pz)
	w   = np.array(np.arange(w2/4.0,4.0*w2,w2/5000.0))
	p   = e2pz(w,w2,th)[0]
	tck = interpolate.UnivariateSpline(p,w)
	w1  = tck(pz)
	return w1

# theory
def e2pz(w1,w2,th):
	"""
	Calculates the momentum scale and the relativistic Compton cross section 
	correction according to P. Holm, PRA 37, 3706 (1988).
	from KH 29.05.96
	input:
	w1 = incident energy  [keV]
	w2 = scattered energy [keV]
	th = scattering angle [deg]
	returns:
	pz = momentum scale   [a.u.]
	cf = cross section correction factor such that:
	J(pz) = cf * d^2(sigma)/d(w2)*d(Omega) [barn/atom/keV/srad]
	"""
	w1  = np.array(w1)    # make sure arrays are used
	w2  = np.array(w2)           
	m   = constants.value('electron mass energy equivalent in MeV')*1e3 #511.003      # Mass in natural units
	th  = math.radians(th) # th/180.0*np.pi  # Angle to radians
	alp = constants.value('fine-structure constant') #1.0/137.036  # Fine structure constant
	r0  = constants.value('classical electron radius') #2.8179e-15   # Electron radius
	q   = np.sqrt(w1**2 + w2**2-2*w1*w2*np.cos(th))                        # Momentum transfer    
	pz  = q/2.0 - (w1-w2) * np.sqrt(1/4.0 + m**2/(2*w1*w2*(1-np.cos(th)))) # In natural units
	E   = np.sqrt(m**2+pz**2)
	A   = ((w1-w2)*E-w1*w2*(1-np.cos(th)))/q
	D   = (w1-w2*np.cos(th))*A/q
	R   = w1*(E-D)
	R2  = R-w1*w2*(1-np.cos(th))
	chi = R/R2 + R2/R + 2.0*m**2 * (1/R-1/R2) + m**4 * (1/R-1/R2)**2
	cf  = 2.0*w1*q*E/(m**2*r0**2*w2*chi)
	cf  = cf*(1e-28*(m*alp))			                             # Cross section now in barns/atom/keV/srad
	pz  = pz/(m*alp)                                                 # pz to atomic units (a.u.)
	return pz, cf

def pz2e1(w2,pz,th):
	"""
	Calculates the incident energy value corresponding
	specific scattered photon and momentum value.
	from KH 29.05.96
	input: 
	w2 = scattered energy [keV]
	pz = momentum value   [a.u.]
	th = scattering angle [deg]
	returns:
	w1 = incident energy  [keV]
	"""
	pz  = np.array(pz)
	w   = np.array(np.arange(w2/4.0,4.0*w2,w2/5000.0))
	p   = e2pz(w,w2,th)[0]
	tck = interpolate.UnivariateSpline(p,w)
	w1  = tck(pz)
	return w1

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
	th    = math.radians(tth)#tth/180.0*np.pi
	hbarc = 137.03599976
	q     = 1/hbarc*np.sqrt(e1**2.0+e2**2.0-2.0*e1*e2*np.cos(th));
	return q

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

def makepzprofile(element,filename='./xrstools/things/ComptonProfiles.dat'):
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
		intf             = interpolate.splrep(roughtheory[:,0],roughtheory[:,n+2],s=0.000000001,k=3) # skip the column with the total J for now #try interp1d with bounds_error=False and fill_value=0.0
		pzprofile[:,n+1] = interpolate.splev(abs(pzscale),intf,der=0)
	# normalize to one electron, multiply by number of electrons
	for n in range(len(binden)):
		normval = integrate.trapz(pzprofile[:,n+1],pzprofile[:,0])
		pzprofile[:,n+1] = pzprofile[:,n+1]/normval*long(occupation[n])
	binden     = [float(en) for en in binden]
	occupation = [float(val) for val in occupation]
	return pzprofile, binden, occupation

def makeprofile(element,filename='./xrstools/things/ComptonProfiles.dat',E0=9.69,tth=35.0,correctasym=None):
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
	q         = q[:,np.nonzero(enscale.T>=0)[0]]
	enscale   = enscale[:,np.nonzero(enscale.T>=0)[0]]
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
	
def parseformula(formula):
	"""
	parses the constituing elements and stoichiometries from a given formula
	input:
	formula = string of a chemical formula (e.g. 'SiO2', 'Ba8Si46', etc.)
	returns:
	elements        = list of constituting elemental symbols
	stoichiometries = list of according stoichiometries in same order as 'elements'
	"""
	elements = []
	stoichiometries = []
	splitted = findall(r'([A-Z][a-z]*)(\d*)',formula)
	elements.extend([element[0] for element in splitted])
	stoichiometries.extend([(int(element[1]) if element[1] else 1) for element in splitted])
	return elements,stoichiometries

def makeprofile_comp(formula,filename='./xrstools/things/ComptonProfiles.dat',E0=9.69,tth=35,correctasym=None):
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

def makeprofile_compds(formulas,concentrations=None,filename='./xrstools/things/ComptonProfiles.dat',E0=9.69,tth=35.0,correctasym=None):
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
	"""
	returns first correction to a 1s, 2s, and 2p compton profiles after Holm and
	Ribberfors. 
	INPUT: pzprofile (e.g. tabulated from Biggs), 
	       occupation (electron configuration), 
	       q (momentum transfer in a.u. )
	OUTPUT: asymmetries to be added to the raw profiles (normalized to the number of electrons on pz scale)
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
		print 'WARNING: cannot calculate HR correction for unfilled 1s shell!'
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
			print 'WARNING: cannot calculate HR correction for unfilled 2s shell!'
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
			print 'WARNING: cannot calculate HR correction for unfilled 2p shell!'
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

def element(z):
	"""
	returns atomic number of given element, if z is a string of the 
	element symbol or string of element symbol of given atomic number z
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
			print 'Given element ' + z + ' unknown'
	elif isinstance(z,int):
		if z > 0 and z < 105:
			Z = zs[z-1]
		else:
			print 'Element Z = '+ str(z) +' unknown'
	else:
		print 'type '+ type(z) + 'not supported'	
	return Z

def myprho(energy,Z,logtablefile='./xrstools/things/logtable.dat'):
	"""
	calculates the photoelectric, elastic, and inelastic absorption of 
	an element Z (can be atomic number or symbol) on the energyscal e1
	input:
	energy = energy scale [keV] (i am guessing keV only right now)
	Z      = atomic number or string of element symbol
	returns: 
	murho = absorption coefficient normalized by the density
	rho   = density
	m     = atomic mass
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
	return murho,rho,m

def mpr(energy,compound):
	"""
	calculates the photoelectric, elastic, and inelastic absorption of 
	an element Z (can be atomic number or symbol) on the energyscal e1
	input:
	energy = energy scale [keV]
	Z      = atomic number or string of element symbol
	returns: 
	murho = absorption coefficient normalized by the density
	rho   = density
	m     = atomic mass
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
	return mr,rhov,mv

#print mpr(1,'SiO2')[0]

def mpr_compds(energy,formulas,concentrations,E0,rho_formu):
	"""
	calculates the photoelectric, elastic, and inelastic absorption of 
	a mix of compounds (list of strings with chemical formulas and list 
	of concentrations) on the energyscal e1
	input:
	energy = energy scale [keV]
	Z      = atomic number or string of element symbol
	returns: 
	murho = absorption coefficient normalized by the density
	rho   = density
	m     = atomic mass
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

#print mpr_compds([1,2,3],['SiO2','CO2'],[0.5,0.5],9.86,[1 1])

def abscorr2(mu1,mu2,alpha,beta,samthick):
	"""
	Calculates absorption correction for given mu1 and mu2
	Multiply the measured spectrum with this correction factor
	corrfac= abscorr2(mu1,mu2,alpha,beta,l)
	mu1    = absorption coefficient for the incident energy (1/cm)
	mu2    = absorption coefficient for the scattered energy (1/cm)
	alpha  = incident angle relative to plane normal (deg)
	beta   = exit angle relative to plane normal (deg)
	           (for transmission geometry use beta < 0)
	l      = sample thickness (cm)
	KH 30.05.96
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

# xrs_read
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
	length = 0
	for line in f:
		length += len(line)
	headerlength = (length-size)/2			
	# get the data
	f = open(filename,'rb')
	predata = arr.array('H')
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
	allscans = sorted(allscans,key=lambda x:x.gettype())
	rawgroups = []
	results = groupby(allscans,lambda x:x.gettype())
	print 'the following scangroups were found:'
	for key,thegroup in results:
		print key
		rawgroups.append(list(thegroup))
	return rawgroups

def makegroup(groupofscans,grouptype=None):
	"""
	takes a group of scans, sums up the signals and monitors, estimates poisson errors, and returns an instance of the scangroup class (turns several instances of the "scan" class into an instance of the "scangroup" class)
	"""
	if not grouptype:
		grouptype = groupofscans[0].gettype() # the scantype of the sum of scans will be the same as the first from the list	
	theenergy   = groupofscans[0].energy # all scans are splined onto energy grid of the first scan
	thesignals  = np.zeros(groupofscans[0].getshape())
	theerrors   = np.zeros(groupofscans[0].getshape())
	themonitors = np.zeros(np.shape(groupofscans[0].monitor))
	# add up signals from each scan for each ROI
	for scan in groupofscans:
		f  = interpolate.interp1d(scan.energy,scan.monitor, bounds_error=False, fill_value=0.0)
		moni = f(theenergy)
		themonitors += moni
		for n in range(thesignals.shape[-1]):
			f = interpolate.interp1d(scan.energy,scan.signals[:,n], bounds_error=False, fill_value=0.0)
			signal = f(theenergy)
			thesignals[:,n] += signal
	# Poisson errors
	for n in range(thesignals.shape[-1]):
		theerrors[:,n]  = np.sqrt(thesignals[:,n])
	# and normalize	to the monitor
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
		grouptype = groupofscans[0].gettype() # the scantype of the sum of scans will be the same as the first from the list	
	theenergy   = groupofscans[0].energy
	thesignals  = np.zeros(groupofscans[0].getshape())
	theerrors   = np.zeros(groupofscans[0].getshape())
	themonitors = np.zeros(np.shape(groupofscans[0].monitor))
	# add up signals from each scan for each ROI
	for scan in groupofscans:
		themonitors += scan.monitor
		for n in range(thesignals.shape[-1]):
			thesignals[:,n] += scan.signals[:,n]
	# Poisson errors
	for n in range(thesignals.shape[-1]):
		theerrors[:,n]  = np.sqrt(thesignals[:,n])
	# and normalize to the monitor
	for n in range(thesignals.shape[-1]):
		thesignals[:,n] = thesignals[:,n]/themonitors
		theerrors[:,n]  = theerrors[:,n]/themonitors
	group = scangroup(theenergy,thesignals,theerrors,grouptype)
	return group

def makegroup_test(groupofscans,grouptype=None,elasticnum=None,cenom=None):
	"""
	takes a group of scans, sums up the signals and monitors, estimates poisson errors, and returns an instance of the scangroup class (turns several instances of the "scan" class into an instance of the "scangroup" class), uses the center of mass of the elastic scan closest in scannumber (but smaller) to correct each scan in energy loss to it's according elastic 
	"""
	if not grouptype:
		grouptype = groupofscans[0].gettype() # the scantype of the sum of scans will be the same as the first from the list

	theenergy   = groupofscans[0].energy
	theeloss    = np.zeros_like(groupofscans[0].energy)
	thesignals  = np.zeros(groupofscans[0].getshape())
	theerrors   = np.zeros(groupofscans[0].getshape())
	themonitors = np.zeros(np.shape(groupofscans[0].monitor))
	theresolution = []
	# go through all scans
	for scan in groupofscans:
		# do scan by scan energy correction if elastic numbers are provided
		if elasticnum:
			if grouptype == 'elastic':
				cenom_all = scan.cenom # elastic scans have their own E0
				theresolution.append(scan.resolution) # elastic scans have resolution
			# find the elastic that is closest but smaller in scannumber
			else:
				ind = np.where(np.array(elasticnum)<scan.number)[0]
				if len(ind)>0:
					elastic_ind = min(range(len(np.array(elasticnum)[ind])), key=lambda i: abs(elasticnum[i]-scan.number))
					cenom_all   = cenom[elastic_ind]
		# if no elasticnumber are provided
		else:
			cenom_all   = 0 # needs changing
		# sum up the monitor
		f  = interpolate.interp1d(scan.energy,scan.monitor, bounds_error=False, fill_value=0.0)
		moni = f(theenergy)
		themonitors += moni
		# master eloss scale
		theeloss = (groupofscans[0].energy-cenom_all[0])*1.0e3
		# go through all ROIs and spline onto master eloss scale
		for n in range(thesignals.shape[-1]):
			x = (scan.energy-cenom_all[n])*1.0e3 # insert trivial values far below and far above energy to avoid interpolation artefacts
			x = np.insert(x,0,-1e10,)
			x = np.insert(x,-1,1e10)
			y = scan.signals[:,n]
			y = np.insert(y,0,0)
			y = np.insert(y,-1,0)
			f = interpolate.interp1d(x,y,bounds_error=False, fill_value=0.0)
			signal = f(theeloss)
			thesignals[:,n] += signal
	# Poisson errors
	for n in range(thesignals.shape[-1]):
		theerrors[:,n]  = np.sqrt(thesignals[:,n])
	# and normalize	to the monitor
	for n in range(thesignals.shape[-1]):
		thesignals[:,n] = thesignals[:,n]/themonitors
		theerrors[:,n]  = theerrors[:,n]/themonitors
	#avarage over all resolutions if there were elastic scans involved
	if theresolution:
		resolution = np.mean(np.array(theresolution),axis=0)
	else:
		resolution = []
	group = scangroup(theenergy,thesignals,theerrors,grouptype,eloss=theeloss,resolution=resolution)
	return group

def appendscans(groups):
	"""
	append groups of scans ordered by their first energy value. long scans are inserted into gaps that at greater than four times the 		grid of the finer scans
	"""
	# find all types of groups	
	grouptypes = [key for key in groups.keys()]
	# deal with cases where there is a long scan
	if 'long' in grouptypes:
		allgroups = [] # all groups BUT the long one
		for grouptype in grouptypes:
			if grouptype != 'long':
				allgroups.append(groups[grouptype])
		allgroups.sort(key = lambda x:x.getestart()) # sort all groups but long by their start-energy
		
		# go through all groups and see if there is a gap between them, if yes, see if the long scan can be inserted into the gap
		# first things before the first group
		theenergy = allgroups[0].energy
		thesignals = allgroups[0].signals
		theerrors = allgroups[0].errors
		inds = np.where(groups['long'].energy<theenergy[0])
		theenergy  = np.append(groups['long'].energy[inds],theenergy)
		thesignals = np.append(np.squeeze(groups['long'].signals[inds,:]),thesignals,0)
		theerrors  = np.append(np.squeeze(groups['long'].errors[inds,:]),theerrors,0)
		# now go through the other groups and see if there is a gap and part of the long scan to insert
		for n in range(len(allgroups[1:])):
			if allgroups[n+1].getestart()-allgroups[n].geteend() > 2.0*allgroups[n+1].getmeanegridspacing(): # worry only if scan groups are far apart
				inds = np.where(np.logical_and(groups['long'].energy>allgroups[n].geteend(),groups['long'].energy<allgroups[n+1].getestart()))
				theenergy  = np.append(theenergy,groups['long'].energy[inds])
				thesignals = np.append(thesignals,np.squeeze(groups['long'].signals[inds,:]),0)
				theerrors  = np.append(theerrors,np.squeeze(groups['long'].errors[inds,:]),0)
			theenergy  = np.append(theenergy,allgroups[n+1].energy)
			thesignals = np.append(thesignals,allgroups[n+1].signals,0)
			theerrors  = np.append(theerrors,allgroups[n+1].errors,0)
		# see if there is more long scan after the last group
		inds = np.where(groups['long'].energy>theenergy[-1])
		theenergy  = np.append(theenergy,groups['long'].energy[inds])
		thesignals = np.append(thesignals,np.squeeze(groups['long'].signals[inds,:]),0)
		theerrors  = np.append(theerrors,np.squeeze(groups['long'].errors[inds,:]),0)

	# deal with the cases where there is no long scan	
	if not 'long' in grouptypes:
		allgroups = []
		for group in groups:
			allgroups.append(groups[group])
		allgroups.sort(key = lambda x:x.getestart()) # sort the groups by their start-energy
		# energy scale
		theenergy = allgroups[0].energy
		for n in range(len(allgroups[1:])):
			theenergy = np.append(theenergy,allgroups[n+1].energy,0)
		# signals
		thesignals = allgroups[0].signals
		for n in range(len(allgroups[1:])):
			thesignals = np.append(thesignals,allgroups[n+1].signals,0)
		# errors
		theerrors = allgroups[0].errors
		for n in range(len(allgroups[1:])):
			theerrors = np.append(theerrors,allgroups[n+1].errors,0)
	return theenergy, thesignals, theerrors

def appendscans_test(groups):
	"""
	append groups of scans ordered by their first energy value. long scans are inserted into gaps that at greater than four times the 		grid of the finer scans
	"""
	# find all types of groups	
	grouptypes = [key for key in groups.keys()]
	# deal with cases where there is a long scan
	if 'long' in grouptypes:
		allgroups = [] # all groups BUT the long one
		for grouptype in grouptypes:
			if grouptype != 'long':
				allgroups.append(groups[grouptype])
		allgroups.sort(key = lambda x:x.getestart()) # sort all groups but long by their start-energy
		
		# go through all groups and see if there is a gap between them, if yes, see if the long scan can be inserted into the gap
		# first things before the first group
		theenergy  = allgroups[0].energy
		theeloss   = allgroups[0].eloss
		thesignals = allgroups[0].signals
		theerrors  = allgroups[0].errors
		inds1       = np.where(groups['long'].energy<theenergy[0])
		theenergy  = np.append(groups['long'].energy[inds1],theenergy)
		inds2       = np.where(groups['long'].eloss<theenergy[0])
		print inds1 == inds2
		theeloss   = np.append(groups['long'].eloss[inds2],theeloss)
		thesignals = np.append(np.squeeze(groups['long'].signals[inds,:]),thesignals,0)
		theerrors  = np.append(np.squeeze(groups['long'].errors[inds,:]),theerrors,0)
		# now go through the other groups and see if there is a gap and part of the long scan to insert
		for n in range(len(allgroups[1:])):
			if allgroups[n+1].getestart()-allgroups[n].geteend() > 2.0*allgroups[n+1].getmeanegridspacing(): # worry only if scan groups are far apart
				inds = np.where(np.logical_and(groups['long'].energy>allgroups[n].geteend(),groups['long'].energy<allgroups[n+1].getestart()))
				theenergy  = np.append(theenergy,groups['long'].energy[inds])
				thesignals = np.append(thesignals,np.squeeze(groups['long'].signals[inds,:]),0)
				theerrors  = np.append(theerrors,np.squeeze(groups['long'].errors[inds,:]),0)
			theenergy  = np.append(theenergy,allgroups[n+1].energy)
			thesignals = np.append(thesignals,allgroups[n+1].signals,0)
			theerrors  = np.append(theerrors,allgroups[n+1].errors,0)
		# see if there is more long scan after the last group
		inds = np.where(groups['long'].energy>theenergy[-1])
		theenergy  = np.append(theenergy,groups['long'].energy[inds])
		thesignals = np.append(thesignals,np.squeeze(groups['long'].signals[inds,:]),0)
		theerrors  = np.append(theerrors,np.squeeze(groups['long'].errors[inds,:]),0)

	# deal with the cases where there is no long scan	
	if not 'long' in grouptypes:
		allgroups = []
		for group in groups:
			allgroups.append(groups[group])
		allgroups.sort(key = lambda x:x.getestart()) # sort the groups by their start-energy
		# energy scale
		theenergy = allgroups[0].energy
		theeloss  = allgroups[0].eloss
		for n in range(len(allgroups[1:])):
			theenergy = np.append(theenergy,allgroups[n+1].energy,0)
			theeloss  = np.append(theeloss,allgroups[n+1].eloss,0)
		# signals
		thesignals = allgroups[0].signals
		for n in range(len(allgroups[1:])):
			thesignals = np.append(thesignals,allgroups[n+1].signals,0)
		# errors
		theerrors = allgroups[0].errors
		for n in range(len(allgroups[1:])):
			theerrors = np.append(theerrors,allgroups[n+1].errors,0)
		# resolution (stored in the elastic line group)
		if allgroups[0].grouptype == 'elastic':
			theresolution = allgroups[0].resolution
		else:
			theresolution = []
	return theenergy, theeloss, thesignals, theerrors, theresolution

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
	hbar = 6.626e-34/2/pi

	e1    = np.array(e1*1e3*e/c/hbar)
	e2    = np.array(e2*1e3*e/c/hbar)
	th    = np.radians(tth)
	q     = np.sqrt(e1**2+e2**2-2*e1*e2*np.cos(phi))/1e10
	return q

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

class scan:
	"""
	this is a container class for inelastic x-ray scattering scans with 2D detectors	
	each scan is an instance of this class, the scans of an experiment can then be 
	grouped and ordered by this class's attributes scantype, starting energy, etc. 
	"""
	def __init__(self,mats,num,en,monoa,moni,counts,mots,data,scantype='generic'):
		# rawdata
		self.edfmats  = np.array(mats)
		self.number   = num
		self.scantype = scantype
		self.energy   = np.array(en)
		self.monoangle= np.array(monoa)
		self.monitor  = np.array(moni)
		# some things maybe not imediately necessary
		self.counters = np.array(counts)
		self.motors   = mots
		self.specdata = np.array(data)
		# data (to be filled after defining rois)
		self.eloss    = []
		self.signals  = []
		self.errors   = []
		self.cenom    = [] # center of mass (only filled if scantype is 'elastic')
		self.resolution = [] # resolution (only filled if scantype is 'elastic')
		
	def applyrois(self,rois):
		"""
		sums up each 2D matrix of a scan over the indices given in rois,
		i.e. turns the 3D matrix of a scan (stack of 2D detector images)
		into a matrix of size (len(energy),number of rois)
		rois are a list of tuples
		"""
		data       = np.zeros((len(self.edfmats),len(rois)))
		cenom      = []
		resolution = []
		edfmats    = self.edfmats
		for n in range(len(rois)): # each roi (default is 9)
			for m in range(len(edfmats)): # each energy point along the scan
				for l in range(len(rois[n])): # each pixel on the detector
					data[m,n] += edfmats[m,rois[n][l][1],rois[n][l][0]]
		
		# determine the center of mass adn resolution, if it's an elastic scan
		if self.scantype == 'elastic':
			for n in range(len(rois)): # each roi
				cenom.append(trapz(data[:,n]*self.energy,x=self.energy)/trapz(data[:,n],x=self.energy))
				try:
					FWHM,x0 = fwhm((self.energy - cenom[n])*1e3,data[:,n])
					resolution.append(FWHM)
				except:
					resolution.append(0)
		# pass variables back to the class attributes
		self.signals    = np.array(data)
		self.cenom      = cenom
		self.resolution = resolution

	def applyrois_old(self,rois):
		"""
		sums up each 2D matrix of a scan over the indices given in rois,
		i.e. turns the 3D matrix of a scan (stack of 2D detector images)
		into a matrix of size (len(energy),number of rois)
		rois are a list of tuples
		this seems a bit faster than the old version
		"""
		data     = np.zeros((len(self.edfmats),len(rois)))
		roixinds = []
		roiyinds = []
		cenom   = []
		resolution = []
		edfmats = self.edfmats
		for r in range(len(rois)):
			for n in range(len(rois[r])):
				roixinds.append(rois[r][n][0])
				roiyinds.append(rois[r][n][1])
			data[:,r] = np.sum(np.sum(edfmats[:,np.amin(roiyinds):np.amax(roiyinds)+1,np.amin(roixinds):np.amax(roixinds)+1],axis=1),axis=1)
		# determine the center of mass, if it's an elastic scan
		if self.scantype == 'elastic':
			for n in range(len(rois)): # each roi
				cenom.append(trapz(data[:,n]*self.energy,x=self.energy)/trapz(data[:,n],x=self.energy))
				try:
					FWHM,x0 = fwhm((self.energy - cenom[n])*1e3,data[:,n])
					resolution.append(FWHM)
				except:
					resolution.append(0.0) # need something more sofisticated to fit FWHM
		self.signals = np.array(data)
		self.cenom   = cenom
		self.resolution = resolution

	def gettype(self):
		return self.scantype

	def getnumber(self):
		return self.number

	def getshape(self):
		if not self.signals.any():
			print 'please apply the ROIs first.'
			return
		else:
			return np.shape(self.signals)

	def getnumofrois(self):
		if not self.signals.any():
			print 'please apply the ROIs first.'
			return
		else:
			return np.shape(self.signals)[1]

class scangroup:
	"""
	container class for scans with the same 'scantype' attribute
	each group of scans will summed into an instance of this class
	different goups of scans will then be stitched together based on 
	their type, starting energy, etc. 
	"""
	def __init__(self,energy,signals,errors,grouptype='generic',eloss=[],resolution=[]):
		self.energy    = energy
		self.eloss     = eloss
		self.signals   = signals
		self.errors    = errors
		self.grouptype = grouptype
		self.resolution = resolution

	def gettype(self):
		return self.grouptype

	def getmeanenergy(self):
		return np.mean(self.energy)

	def getestart(self):
		return self.energy[0]

	def geteend(self):
		return self.energy[-1]

	def getmeanegridspacing(self):
		return np.mean(np.diff(self.energy))

	def getmaxediff(self):
		return (self.energy[-1]-self.energy[0])

# superresolution/imaging

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

