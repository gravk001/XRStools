#!/usr/bin/python
# Filename: stobe_analyze.py

import os
import warnings
from copy import deepcopy

import numpy as np
import array as arr

from itertools import groupby
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from scipy import constants

from pylab import *
from scipy import signal
from scipy.ndimage import measurements
import matplotlib.pyplot as plt

__metaclass__ = type # new style classes

def gauss1(x,x0,fwhm):
	"""
	returns a gaussian with peak value normalized to unity
	a[0] = peak position
	a[1] = Full Width at Half Maximum
	"""
	y = np.exp(-np.log(2.0)*((x-x0)/fwhm*2.0)**2.0)
	return y

def gauss(x,x0,fwhm):
    # area-normalized gaussian
    sigma = fwhm/(2*np.sqrt(2*np.log(2)));
    y = np.exp(-(x-x0)**2/2/sigma**2)/sigma/np.sqrt(2*np.pi)
    return y

def gauss_areanorm(x,x0,fwhm):
	"""
	area-normalized gaussian
	"""
	sigma = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
	y = np.exp(-(x-x0)**2.0/2.0/sigma**2)/sigma/np.sqrt(2.0*np.pi)
	return y

def convg(x,y,fwhm):
	"""
	Convolution with Gaussian	
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
        s2 = s2 + sticks[n]*gauss1(e2,evals[n],fwhm[n])
    spectrum = np.zeros((len(e2),2))
    spectrum[:,0] = e2
    spectrum[:,1] = s2
    return spectrum

def load_stobe_specs(prefix,postfix,fromnumber,tonumber,step,stepformat=2):
	"""
	load a bunch of StoBe calculations, which filenames are made up of the prefix, postfix, and the counter in the between the prefix and postfix runs from 'fromnumber' to 'tonumber' in steps of 'step' (number of digits is 'stepformat')
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
    spec = np.loadtxt(filename)
    return spec

def load_erkale_specs(prefix,postfix,fromnumber,tonumber,step,stepformat=2):
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
    if not emin:
        emin = spec[0,0]
    if not emax:
        emax = spec[-1,0]
    spec = spec[spec[:,0]>emin]
    spec = spec[spec[:,0]<emax]
    return spec



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
		s2 += s[i]*gauss1(e2,e[i],fwhm[i])

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


class stobe:
	"""
	class to analyze StoBe results
	"""
	def __init__(self,prefix,postfix,fromnumber,tonumber,step,stepformat=2):
		self.energy   = [] # array of final energy scale for all snapshots of this run
		self.signal   = [] # array of averaged and smoothed results     
		self.rawspecs = load_stobe_specs(prefix,postfix,fromnumber,tonumber,step,stepformat) # list of all raw stick spectra:rawspecs[n][energy,sticks]
		self.broadened= [] # list of broadened stick spectra 

	def cut_rawspecs(self,emin=None,emax=None):
		cutspecs = []
		for spec in self.rawspecs:
			cutspecs.append(cut_spec(spec,emin,emax))
		self.rawspecs = cutspecs

	def broaden_lin(self,params=[0.8, 8, 537.5, 550],npoints=1000):
		for spec in self.rawspecs:
			self.broadened.append(broaden_linear(spec,params,npoints))

	def sum_specs(self):
		self.energy   = self.broadened[0][:,0] # first spectrum defines energy scale
		self.signal   = np.zeros(np.shape(self.energy))
		for spec in self.broadened:
			f = interp1d(spec[:,0], spec[:,1],bounds_error=False, kind='cubic', fill_value=0.0)
			self.signal += f(self.energy)

	def norm_area(self,emin,emax):
		inds = np.where(np.logical_and(self.energy>=emin,self.energy<=emax))
		norm = np.trapz(self.signal[inds],self.energy[inds])
		self.signal = self.signal/norm
			

class erkale:
	"""
	class to analyze ERKALE results
	"""
	def __init__(self,prefix,postfix,fromnumber,tonumber,step,stepformat=2):
		self.energy   = [] # array of final energy scale for all snapshots of this run
		self.sqw      = [] # array of averaged and smoothed results        
		self.rawspecs = load_erkale_specs(prefix,postfix,fromnumber,tonumber,step,stepformat) # list of all raw stick spectra
		self.broadened= []
		self.norm     = [] # results of normalization

	def cut_rawspecs(self,emin=None,emax=None):
		cutspecs = []
		for spec in self.rawspecs:
			cutspecs.append(cut_spec(spec,emin,emax))
		self.rawspecs = cutspecs

	def cut_broadspecs(self,emin=None,emax=None):
		cutspecs = []
		for spec in self.broadened:
			cutspecs.append(cut_spec(spec,emin,emax))
		self.broadened = cutspecs

	def broaden_lin(self,params=[0.8, 8, 537.5, 550],npoints=1000):
		for spec in self.rawspecs:
			self.broadened.append(broaden_linear(spec,params,npoints))

	def sum_specs(self):
		self.energy   = self.broadened[0][:,0] # first spectrum defines energy scale
		self.sqw      = np.zeros(np.shape(self.energy))
		for spec in self.broadened:
			f = interp1d(spec[:,0], spec[:,1],bounds_error=False, kind='cubic', fill_value=0.0)
			self.sqw += f(self.energy)

	def norm_area(self,emin=None,emax=None):
		if not emin:
			emin = self.energy[0]
		if not emax:
			emax = self.energy[-1]

		inds = np.where(np.logical_and(self.energy>=emin,self.energy<=emax))[0]
		self.sqw = self.sqw/np.trapz(self.sqw[inds],self.energy[inds])
    
	def norm_max(self):
		pass

	def plot_spec(self):
		plt.plot(self.energy,self.sqw)
		plt.show()


#### calculations
# cut clusters from xyz files
# cut clusters from gro files

################################### 
# reading function for cowan's code output

class xyzAtom:
        def __init__(self,name,coordinates,number):
                self.name        = name
                self.coordinates = np.array(coordinates)
                self.x_coord     = self.coordinates[0]
                self.y_coord     = self.coordinates[1]
                self.z_coord     = self.coordinates[2]
                self.number      = number

        def getNorm(self):
                return np.linalg.norm(self.coordinates)

	def getCoordinates(self):
		return self.coordinates

	def translateSelf(self, vector):
		try:
			self.coordinates += vector
        	        self.x_coord     += vector[0]
        	        self.y_coord     += vector[1]
        	        self.z_coord     += vector[2]
		except ValueError:
			print('Vector must be 3D np.array!') 

class xyzMolecule:
        def __init__(self,xyzAtoms):
                self.xyzAtoms = xyzAtoms

        def getCoordinates(self):
                return [atom.getCoordinates() for atom in self.xyzAtoms]

        def getCoordinates_name(self,name):
                atoms = []
                for atom in self.xyzAtoms:
                        if atom.name == name:
                                atoms.append(atom.coordinates)
                return atoms

        def get_atoms_by_name(self,name):
                atoms = []
                for atom in self.xyzAtoms:
                        if atom.name == name:
                                atoms.append(atom)
                        else:
                                pass
                if len(atoms) == 0:
                        print('Found no atoms with given name in box.')
                        return
                return atoms

        def getGeometricCenter(self):
                for_average = np.zeros((len(self.xyzAtoms),3))
                for ii in range(len(self.xyzAtoms)):
                        for_average[ii, :] = self.xyzAtoms[ii].coordinates
                return np.mean(for_average,axis = 0)

	def translateSelf(self,vector):
		for atom in self.xyzAtoms:
			atom.translateSelf(vector)

	def scatterPlot(self):
		from mpl_toolkits.mplot3d import Axes3D
		fig = figure()
		ax  = Axes3D(fig)
		x_vals = [coord[0] for coord in self.getCoordinates()]
		y_vals = [coord[1] for coord in self.getCoordinates()]
		z_vals = [coord[2] for coord in self.getCoordinates()]
		ax.scatter(x_vals, y_vals, z_vals)
		show()

class xyzBox:
        def __init__(self,xyzAtoms,boxLength=None):
                self.xyzMolecules = []
                self.xyzAtoms     = xyzAtoms
                self.n_atoms      = len(self.xyzAtoms)
                self.title        = None
		self.boxLength    = boxLength

	def setBoxLength(self,boxLength,angstrom=True):
		if angstrom:
			self.boxLength = boxLength
		else:
			self.boxLength = boxLength*constants.physical_constants['atomic unit of length'][0]*10**10

        def writeBox(self, filename):
                writeXYZfile(filename, self.n_atoms, self.title, self.xyzAtoms)

        def getCoordinates(self):		
                return [atom.getCoordinates() for atom in self.xyzAtoms]

        def get_atoms_by_name(self,name):
                # find oxygen atoms
                atoms = []
                for atom in self.xyzAtoms:
                        if atom.name == name:
                                atoms.append(atom)
                        else:
                                pass
                if len(atoms) == 0:
                        print('Found no atoms with given name in box.')
                        return
                return atoms

	def scatterPlot(self):
		from mpl_toolkits.mplot3d import Axes3D
		fig = figure()
		ax  = Axes3D(fig)
		#coordinates = self.getCoordinates()
		#print type(coordinates), coordinates
		x_vals = [coord[0] for coord in self.getCoordinates()]
		y_vals = [coord[1] for coord in self.getCoordinates()]
		z_vals = [coord[2] for coord in self.getCoordinates()]
		cla()
		ax.scatter(x_vals, y_vals, z_vals)
		draw()

        def get_OO_neighbors(self,Roocut=3.6):
                """
                Returns list of numbers of nearest oxygen neighbors in readius 'Roocut'
                """
                o_atoms = self.get_atoms_by_name('O')
		if not self.boxLength:
        	        return count_OO_neighbors(o_atoms,Roocut)
		else:
	       	        return count_OO_neighbors(o_atoms,Roocut,boxLength=self.boxLength)

	def get_OO_neighbors_pbc(self,Roocut=3.6):
		o_atoms = self.get_atoms_by_name('O')
		return count_OO_neighbors_pbc(o_atoms,Roocut,boxLength=self.boxLength)

	def get_h2o_molecules(self,o_name='O',h_name='H'):
		o_atoms  = self.get_atoms_by_name(o_name)
		h_atoms  = self.get_atoms_by_name(h_name)
		if self.boxLength:
			self.xyzMolecules = find_H2O_molecules(o_atoms,h_atoms,boxLength=self.boxLength)
		else:
			self.xyzMolecules = find_H2O_molecules(o_atoms,h_atoms,boxLength=None)

	def get_atoms_from_molecules(self):
		if not self.xyzAtoms:
			self.xyzAtoms = []
			for molecule in self.xyzMolecules:
				for atom in molecule.xyzAtoms:
					self.xyzAtoms.append(atom)
		self.n_atoms      = len(self.xyzAtoms)

	def get_hbonds(self, Roocut=3.6, Rohcut=2.4, Aoooh=30.0):
		hb_total_sum  = 0 # total number of H-bonds in box
		hb_donor_sum  = 0 # total number of donor H-bonds in box
                o_atoms  = self.get_atoms_by_name('O')
                h_atoms  = self.get_atoms_by_name('H')
		if self.boxLength:
			h2o_mols = find_H2O_molecules(o_atoms,h_atoms,boxLength=self.boxLength)
			self.xyzMolecules = h2o_mols
		else:
			h2o_mols = find_H2O_molecules(o_atoms,h_atoms)
			test_h2o_mols = h2o_mols
		for mol1 in h2o_mols:
			for mol2 in h2o_mols:
				donor, acceptor = countHbonds_pbc(mol1,mol2,self.boxLength,Roocut=Roocut, Rohcut=Rohcut, Aoooh=Aoooh)
				hb_donor_sum  += donor
				hb_accept_sum += acceptor
		return hb_donor_sum, hb_accept_sum #hbonds, dbonds, abonds, hbondspermol

	def changeOHBondlength(self,fraction, oName='O', hName='H'):
		o_atoms  = self.get_atoms_by_name('O')
		h_atoms  = self.get_atoms_by_name('H')
		# find all H2O molecules
		if self.boxLength:
			h2o_mols = find_H2O_molecules(o_atoms,h_atoms,boxLength=self.boxLength)
		else:
			h2o_mols = find_H2O_molecules(o_atoms,h_atoms)
		# change the bond length
		new_h2o_mols = []
		for mol in h2o_mols:
			new_h2o_mols.append(changeOHBondLength(mol, fraction, boxLength=self.boxLength, oName=oName, hName=hName))
		# redefine all molecules and atoms in box
		self.xyzMolecules = new_h2o_mols
		self.xyzAtoms = []
		for mol in self.xyzMolecules:
			for atom in mol.xyzAtoms:
				self.xyzAtoms.append(atom)

def getPeriodicTestBox_molecules(Molecules,boxLength,numbershells=1):
	vectors = []
	for l in range(-numbershells,numbershells+1):
		for m in range(-numbershells,numbershells+1):
			for n in range(-numbershells,numbershells+1):
				vectors.append(np.array([l, m, n]))
	pbc_molecules = []
	for vector in vectors:
		for molecule in Molecules:
			cpAtoms = []
			for atom in molecule.xyzAtoms:
				cpAtom = deepcopy(atom)
				cpAtom.translateSelf(vector*boxLength)
				cpAtoms.append(cpAtom)
			pbc_molecules.append(xyzMolecule(cpAtoms))
	return pbc_molecules

def getPeriodicTestBox(xyzAtoms,boxLength,numbershells=1):
	vectors = []
	for l in range(-numbershells,numbershells+1):
		for m in range(-numbershells,numbershells+1):
			for n in range(-numbershells,numbershells+1):
				vectors.append(np.array([l, m, n]))
	pbc_atoms = []
	for vector in vectors:
		for atom in xyzAtoms:
			cpAtom = deepcopy(atom)
			cpAtom.translateSelf(vector*boxLength)
			pbc_atoms.append(cpAtom)
	return pbc_atoms

def getTranslVec(atom1,atom2,boxLength):
	""" **getTranslVec**
	Returns the translation vector that brings atom2 closer to atom1 in case
	atom2 is further than boxLength away.
	"""
	xcoord1 = atom1.coordinates[0]
	xcoord2 = atom2.coordinates[0]
	ycoord1 = atom1.coordinates[1]
	ycoord2 = atom2.coordinates[1]
	zcoord1 = atom1.coordinates[2]
	zcoord2 = atom2.coordinates[2]
	translVec = np.zeros(atom2.coordinates.shape)
	if xcoord1-xcoord2 > boxLength/2.0:
		translVec[0] = -boxLength/2.0
	if ycoord1-ycoord2 > boxLength/2.0:
		translVec[1] = -boxLength/2.0
	if zcoord1-zcoord2 > boxLength/2.0:
		translVec[2] = -boxLength/2.0
	return translVec

def getTranslVec_geocen(mol1COM,mol2COM,boxLength):
	""" **getTranslVec_geocen**
	"""
	translVec = np.zeros((len(mol1COM),))
	if mol1COM[0]-mol2COM[0] > boxLength/2.0:
		translVec[0] = -boxLength/2.0
	if mol1COM[1]-mol2COM[0] > boxLength/2.0:
		translVec[1] = -boxLength/2.0
	if mol1COM[2]-mol2COM[0] > boxLength/2.0:
		translVec[2] = -boxLength/2.0
	return translVec

def getDistancePbc(atom1,atom2,boxLength):
	xdist  = atom1.coordinates[0] - atom2.coordinates[0]
	xdist -= boxLength*round(xdist/boxLength)
	ydist  = atom1.coordinates[1] - atom2.coordinates[1] 
	ydist -= boxLength*round(ydist/boxLength)
	zdist  = atom1.coordinates[2] - atom2.coordinates[2] 
	zdist -= boxLength*round(zdist/boxLength)
	return np.sqrt(xdist**2.0 + ydist**2.0 + zdist**2.0)

def getDistance(atom1, atom2):
	return np.linalg.norm(atom2.getCoordinates()-atom2.getCoordinates())

def getDistVectorPbc(atom1,atom2,boxLength):
	xdist  = atom1.coordinates[0] - atom2.coordinates[0]
	xdist -= boxLength*round(xdist/boxLength)
	ydist  = atom1.coordinates[1] - atom2.coordinates[1] 
	ydist -= boxLength*round(ydist/boxLength)
	zdist  = atom1.coordinates[2] - atom2.coordinates[2] 
	zdist -= boxLength*round(zdist/boxLength)
	return np.array([xdist, ydist, zdist])

def getDistVector(atom1,atom2):
	xdist  = atom1.coordinates[0] - atom2.coordinates[0]
	ydist  = atom1.coordinates[1] - atom2.coordinates[1] 
	zdist  = atom1.coordinates[2] - atom2.coordinates[2] 
	return np.array([xdist, ydist, zdist])

def countHbonds_pbc(mol1,mol2,boxLength,Roocut=3.6, Rohcut=2.4, Aoooh=30.0):
	hb_donor  = 0
	hb_accept = 0
	# get atoms
	mol1_o = mol1.get_atoms_by_name('O')
	mol1_h = mol1.get_atoms_by_name('H')
	mol2_o = mol2.get_atoms_by_name('O')
	mol2_h = mol2.get_atoms_by_name('H')
	if getDistancePbc(mol1_o[0],mol2_o[0],boxLength) <= Roocut:
		# donor bond through first hydrogen atom of mol1
		if getDistancePbc(mol1_h[0],mol2_o[0],boxLength) <= Rohcut:
			vec1 = getDistVectorPbc(mol1_h[0],mol1_o[0],boxLength)
			vec2 = getDistVectorPbc(mol2_o[0],mol1_o[0],boxLength)
			if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
				hb_donor += 1
		# donor bond through second hydrogen atom of mol1
		if getDistancePbc(mol1_h[1],mol2_o[0],boxLength) <= Rohcut:
			vec1 = getDistVectorPbc(mol1_h[1],mol1_o[0],boxLength)
			vec2 = getDistVectorPbc(mol2_o[0],mol1_o[0],boxLength)
			if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
				hb_donor += 1
		# acceptor bond through first hydrogen atom of mol2
		if getDistancePbc(mol1_o[0],mol2_h[0],boxLength) <= Rohcut:
			vec1 = getDistVectorPbc(mol2_h[0],mol1_o[0],boxLength)
			vec2 = getDistVectorPbc(mol2_o[0],mol1_o[0],boxLength)
			if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
				hb_accept += 1
		# acceptor bond through second hydrogen atom of mol2
		if getDistancePbc(mol1_o[0],mol2_h[1],boxLength) <= Rohcut:
			vec1 = getDistVectorPbc(mol2_h[1],mol1_o[0],boxLength)
			vec2 = getDistVectorPbc(mol2_o[0],mol1_o[0],boxLength)
			if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
				hb_accept += 1
	return hb_donor, hb_accept

def countHbonds(mol1,mol2, Roocut=3.6, Rohcut=2.4, Aoooh=30.0):
	hb_donor  = 0
	hb_accept = 0
	# get the coordinates
	Aoooh = np.radians(Aoooh)
	mol1_o = mol1.getCoordinates_name('O')
	mol1_h = mol1.getCoordinates_name('H')
	mol2_o = mol2.getCoordinates_name('O')
	mol2_h = mol2.getCoordinates_name('H')
	# check O-O distance
	if np.linalg.norm(mol2_o[0] - mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0] - mol1_o[0]) <= Roocut: # check Roocut is met
		# check for donor H-bonds
		if np.linalg.norm(mol2_o[0] - mol1_h[0]) > 0.0 and np.linalg.norm(mol2_o[0] - mol1_h[0]) <= Rohcut: # check Rohcut, first H
			vec1 = mol1_h[0]-mol1_o[0]
			vec2 = mol2_o[0]-mol1_o[0]
			if np.arccos((np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
				hb_donor += 1
		if np.linalg.norm(mol2_o[0] - mol1_h[1]) > 0.0 and np.linalg.norm(mol2_o[0] - mol1_h[1]) <= Rohcut: # check Rohcut, second H
			vec1 = mol1_h[1]-mol1_o[0]
			vec2 = mol2_o[0]-mol1_o[0]
			if np.arccos((np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
				hb_donor += 1
		# check for acceptor H-bonds
		if np.linalg.norm(mol2_h[0] - mol1_o[0]) > 0.0 and np.linalg.norm(mol2_h[0] - mol1_o[0]) <= Rohcut: # check Rohcut, first H
			vec1 = mol2_o[0]-mol1_o[0]
			vec2 = mol2_h[0]-mol1_o[0]
			if np.arccos((np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
				hb_accept += 1
		if np.linalg.norm(mol2_h[1] - mol1_o[0]) > 0.0 and np.linalg.norm(mol2_h[1] - mol1_o[0]) <= Rohcut: # check Rohcut, socond H
			vec1 = mol2_o[0]-mol1_o[0]
			vec2 = mol2_h[1]-mol1_o[0]
			if np.arccos((np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
				hb_accept += 1
	return hb_donor, hb_accept

def countHbonds_orig(mol1,mol2, Roocut=3.6, Rohcut=2.4, Aoooh=30.0):
	mol1_o = mol1.getCoordinates_name('O')
	mol1_h = mol1.getCoordinates_name('H')
	mol2_o = mol2.getCoordinates_name('O')
	mol2_h = mol2.getCoordinates_name('H')
	hbnoD1 = 0
	if np.linalg.norm(mol2_o[0]-mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0]-mol1_o[0]) <= Roocut:
		if np.linalg.norm(mol2_h[0] - mol1_o[0]) <= Rohcut:
			vec1 = mol2_h[0]-mol2_o[0]
			vec2 = mol1_o[0]-mol2_o[0]
			if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
				hbnoD1 += 1
	hbnoD2 = 0
	if np.linalg.norm(mol2_o[0]-mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0]-mol1_o[0]) <= Roocut:
		if np.linalg.norm(mol2_h[1]-mol1_o[0]) <= Rohcut:
			vec1 = mol2_h[1]-mol2_o[0]
			vec2 = mol1_o[0]-mol2_o[0]
			if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
				hbnoD2 += 1
	hbnoA1 = 0
	if np.linalg.norm(mol2_o[0]-mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0]-mol1_o[0]) <= Roocut:
		if np.linalg.norm(mol2_o[0]-mol1_h[0]) <= Rohcut:
			vec1 = mol2_o[0]-mol1_o[0]
			vec2 = mol1_h[0]-mol1_o[0]
			if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
				hbnoA1 += 1
	hbnoA2 = 0
	if np.linalg.norm(mol2_o[0]-mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0]-mol1_o[0]) <= Roocut:
		if np.linalg.norm(mol2_o[0]-mol1_h[1]) <= Rohcut:
			vec1 = mol2_o[0]-mol1_o[0]
			vec2 = mol1_h[1]-mol1_o[0]
			if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
				hbnoA2 += 1
	dbonds = hbnoD1 + hbnoD2
	abonds = hbnoA1 + hbnoA2
	return dbonds, abonds

def repair_h2o_molecules_pbc(h2o_mols,boxLength):
	new_mols = []
	for mol in h2o_mols:
		o_atom  = mol.get_atoms_by_name('O')[0]
		h_atoms = mol.get_atoms_by_name('H')
		new_mol = [o_atom]
		for h_atom in h_atoms:
			cpAtom = deepcopy(h_atom)
			xdist  = o_atom.coordinates[0] - cpAtom.coordinates[0]
			ydist  = o_atom.coordinates[1] - cpAtom.coordinates[1]
			zdist  = o_atom.coordinates[2] - cpAtom.coordinates[2]
			cpAtom.translateSelf([boxLength*round(xdist/boxLength),boxLength*round(ydist/boxLength),boxLength*round(zdist/boxLength)])
			new_mol.append(cpAtom)
		new_mols.append(xyzMolecule(new_mol))
	return new_mols

def find_H2O_molecules(o_atoms,h_atoms,boxLength=None):
	h2o_molecules = []
	if not boxLength:
		warnings.warn('No box length provided, will not take PBC into account!')
		for o_atom in o_atoms:
			ho_dists = []
			for h_atom in h_atoms:
				ho_dists.append(np.linalg.norm(o_atom.coordinates - h_atom.coordinates))
			order = np.argsort(ho_dists)
			if np.linalg.norm(o_atom.coordinates - h_atoms[order[0]].getCoordinates()) <= 1.5 and np.linalg.norm(o_atom.coordinates - h_atoms[order[1]].getCoordinates()) <= 1.5:
				h2o_molecules.append(xyzMolecule([o_atom,h_atoms[order[0]],h_atoms[order[1]]]))
		return h2o_molecules
	else:
		for o_atom in o_atoms:
			ho_dists = []
			for h_atom in h_atoms:
				ho_dists.append(getDistancePbc(o_atom,h_atom,boxLength))
			order = np.argsort(ho_dists)
			h2o_molecules.append(xyzMolecule([o_atom,h_atoms[order[0]],h_atoms[order[1]]]))
		return h2o_molecules

def writeXYZfile(filename,numberOfAtoms, title, list_of_xyzAtoms):
        # create file
        xyz = open(filename,'w+')
        # write number of atoms
        xyz.write(str(numberOfAtoms) + ' \n')
        # write title
	if not title:
		title = 'None'
        xyz.write(title + ' \n')
        # write coordinates
        for atom in list_of_xyzAtoms:
                xyz.write('%4s %8f %8f %8f \n' % (atom.name, atom.x_coord, atom.y_coord, atom.z_coord))
        xyz.close()

def count_OO_neighbors(list_of_o_atoms,Roocut,boxLength=None):
	noo   = []
	if not boxLength:
		warnings.warn('No box length provided, will not take PBC into account!')
	        for atom1 in list_of_o_atoms:
	                dists = []
	                for atom2 in list_of_o_atoms:
	                        dists.append(np.linalg.norm(atom1.coordinates - atom2.coordinates))
	                noo.append(len(np.where(np.logical_and(np.sort(np.array(dists))>0.0, np.sort(np.array(dists))<=Roocut))[0]))
	        return noo
	else:
		for atom1 in list_of_o_atoms:
			dists = []
			for atom2 in list_of_o_atoms:
				dists.append(getDistancePbc(atom1,atom2,boxLength))
	                noo.append(len(np.where(np.logical_and(np.sort(np.array(dists))>0.0, np.sort(np.array(dists))<=Roocut))[0]))
		return noo

def count_OO_neighbors_pbc(list_of_o_atoms,Roocut,boxLength,numbershells=1):
	noo = []
	vectors = []
	for l in range(-numbershells,numbershells+1):
		for m in range(-numbershells,numbershells+1):
			for n in range(-numbershells,numbershells+1):
				vectors.append(np.array([l, m, n]))
	pbc_atoms = []
	for vector in vectors:
		for atom in list_of_o_atoms:
			cpAtom = deepcopy(atom)
			cpAtom.translateSelf(vector*boxLength)
			pbc_atoms.append(cpAtom)
	for atom1 in list_of_o_atoms:
		dists = []
		for atom2 in pbc_atoms:
			dists.append(np.linalg.norm(atom1.coordinates - atom2.coordinates))
		sdists = np.sort(np.array(dists))
		noo.append(len(np.where(np.logical_and(sdists>0.0, sdists<=Roocut))[0]))
	return noo

def boxParser(filename):
        """**parseXYZfile**
        Reads an xyz-style file.
        """
        atoms = []
        coordinates = []
        xyz = open(filename)
        n_atoms = int(xyz.readline())
        title = xyz.readline()
        for line in xyz:
                atom,x,y,z = line.split()
                atoms.append(atom)
                coordinates.append([float(x), float(y), float(z)])
        xyz.close()
        xyzAtoms = []
        for ii in range(n_atoms):
                xyzAtoms.append(xyzAtom(atoms[ii],coordinates[ii],ii))

        return xyzBox(xyzAtoms)

def changeOHBondLength(h2oMol, fraction, boxLength=None, oName='O', hName='H'):
	o_atom = h2oMol.get_atoms_by_name(oName)
	h_atom = h2oMol.get_atoms_by_name(hName)
	# get OH bond-vectors
	ohVectors = []
	if not boxLength:
		warnings.warn('No box length provided, will not take PBC into account!')
		for atom in h_atom:
			ohVectors.append(getDistVector(atom,o_atom[0]))
	else:
		for atom in h_atom:
			ohVectors.append(getDistVectorPbc(atom,o_atom[0],boxLength))
	# get new OH bond vectors
	new_ohVectors = []
	for vector in ohVectors:
		new_ohVectors.append( vector*fraction)
	# make new molecule
	newmol =[deepcopy(o_atom[0])]
	for vector in new_ohVectors:
		newmol.append(xyzAtom(hName,o_atom[0].getCoordinates()+vector,1))
	return xyzMolecule(newmol)


class xyzTrajectory:
        def __init__(self,xyzBoxes):
                self.xyzBoxes = xyzBoxes

        def writeRandBox(self,filename):
                ind = np.random.randint(len(self.xyzBoxes))
                self.xyzBoxes[ind].writeBox(filename)


      

def parseXYZfile(filename):
        """**parseXYZfile**
        Reads an xyz-style file.
        """
        atoms = []
        coordinates = []
        xyz = open(filename)
        n_atoms = int(xyz.readline())
        title = xyz.readline()
        for line in xyz:
                atom,x,y,z = line.split()
                atoms.append(atom)
                coordinates.append([float(x), float(y), float(z)])
        xyz.close()
        return n_atoms, title, atoms, coordinates



