#!/usr/bin/python
# Filename: extraction.py

from helpers import *

import numpy as np
import pylab
import math
from scipy import interpolate, signal, integrate, constants, optimize, ndimage
import matplotlib.pyplot as plt

class extraction:
	"""
	class for extraction of S(q,w) from an instance of the id20 class "data" and the predictthings class "theory"
	"""
	def __init__(self,data,theory,prenormrange=[5,np.inf]):
		# the data	
		self.eloss   = data.eloss
		self.signals = data.signals
		self.errors  = data.errors
		self.E0      = data.E0
		self.tth     = data.tth
		self.prenormrange = prenormrange
		# the theory 
		self.J        = theory.J
		self.C        = theory.C
		self.V        = theory.V
		self.qvals    = theory.q
		self.formulas = theory.formulas
		self.concentrations = theory.concentrations
		# output
		self.background = np.zeros(np.shape(data.signals))
		self.sqw        = np.zeros(np.shape(data.signals))
		self.pzscale      = np.flipud(np.arange(-10,10,0.05)) # definition according to Huotari et al, JPCS 62 (2001) 2205
		self.valencepz      = np.zeros((len(self.pzscale),len(self.signals[0,:])))
		self.valasymmetrypz = np.zeros((len(self.pzscale),len(self.signals[0,:])))
		self.valence      = np.zeros((len(self.eloss),len(self.signals[0,:])))
		self.valasymmetry = np.zeros((len(self.eloss),len(self.signals[0,:])))
		self.sqwav      = np.zeros(np.shape(data.eloss))
		self.sqwaverr   = np.zeros(np.shape(data.eloss))
		# rough normalization over range given by prenormrange
		for n in range(len(self.signals[0,:])):
			HFnorm = np.trapz(self.J[:,n],self.eloss)
			inds   = np.where(np.logical_and(self.eloss>=prenormrange[0],self.eloss<=prenormrange[1]))[0]
			EXPnorm = np.trapz(self.signals[inds,n],self.eloss[inds])
			self.signals[:,n] = self.signals[:,n]/EXPnorm*HFnorm
			self.errors[:,n]  = self.errors[:,n]/EXPnorm*HFnorm

	def areanorm(self,whichq,emin=None,emax=None):
		"""
		normalizes self.signals to area in between emin and emax, 
		default values cover the whole self.eloss axis
		"""		
		cols = []		
		if not isinstance(whichq,list):
			cols.append(whichq)
		else:
			cols = whichq

		if not emin:
			emin = self.eloss[0]
		if not emax:
			emax = self.eloss[-1]

		for col in cols:
			inds = np.where(np.logical_and(self.eloss>=emin,self.eloss<=emax))
			self.signals[:,col] = self.signals[:,col]/np.trapz(self.signals[inds,col],self.eloss[inds])

	def energycorrect(self,whichq,alpha,densities,samthickness):
		"""
		apply energy dependent corrections to the measured data based on scattering angles, sample material, ... 
		whichq       = single value or list of columns to apply the corrections to (index starts at zero)
		alpha        = incident beam angle (relative to sample surface normal), need to be negative for transmission geometry
		density      = single value (just one compound) or list of values (mixture of several compounds) 
		samthickness = sample thickness in [cm]		
		"""
		# make things iterable
		cols = []		
		if not isinstance(whichq,list):
			cols.append(whichq)
		else:
			cols = whichq
		denses = []		
		if not isinstance(densities,list):
			denses.append(densities)
		else:
			denses = densities

		# calculate beta (exit angle) from alpha and tth for all columns in whichq
		beta = []
		
		if alpha >=0: # reflection geometry (check this again)
			for n in range(len(cols)):
				beta.append(np.absolute(180.0 - alpha - self.tth[cols[n]]))
		else: # transmission geometry (check this again)
			for n in range(len(cols)):
				beta.append(-1.0*(np.absolute(np.absolute(alpha) - self.tth[cols[n]])))

		# absorption and self-absorption
		# mu_in and mu_out from log-log table
		mu_in, mu_out = mpr_compds(self.eloss/1e3+self.E0,self.formulas,self.concentrations,self.E0,denses)
		ac = np.zeros((len(self.eloss),len(cols)))
		for n in range(len(cols)):
			ac[:,n] = abscorr2(mu_in,mu_out,alpha,beta[n],samthickness)

		# cross section correction (use cf output from e2pz routine)
		pz = np.zeros((len(self.eloss),len(cols)))
		cf = np.zeros((len(self.eloss),len(cols)))
		for n in range(len(cols)):
			pz[:,n], cf[:,n] = e2pz(self.E0+self.eloss/1.0e3,self.E0,self.tth[n])

		for col in cols:
			self.signals[:,col] = self.signals[:,col]*ac[:,n]*cf[:,n]
			self.errors[:,col]  = self.errors[:,col]*ac[:,n]*cf[:,n]
			# normalize back to HF profiles (because we know those are on eV / [1/eV] scale (is this correct)
			HFnorm = np.trapz(self.J[:,col],self.eloss)
			inds   = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]
			EXPnorm = np.trapz(self.signals[inds,col],self.eloss[inds])
			self.signals[:,col] = self.signals[:,col]/EXPnorm*HFnorm
			self.errors[:,col]  = self.errors[:,col]/EXPnorm*HFnorm

	def removeelastic(self,whichq,range1,range2,guess=None,stoploop=True,overwrite=False):
		"""
		subtract a pearson function before starting extraction procedure, e.g. for subtracting the elastic peak tail
		guess values: 
		a[0] = Peak position
		a[1] = FWHM
		a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
		a[3] = Peak intensity
		a[4] = background
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		region1 = np.where(np.logical_and(self.eloss>=range1[0],self.eloss<=range1[1]))[0]
		region2 = np.where(np.logical_and(self.eloss>=range2[0],self.eloss<=range2[1]))[0]
		region  = np.append(region1,region2)

		plt.ion()
		for col in columns:
			if not guess: 
				guess = []
				ind = self.signals[:,col].argmax(axis=0) 
				guess.append(self.eloss[ind]) # max of signal (in range of prenorm from __init__)
				guess.append(1.0) # twice the position of the peason maximum
				guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
				guess.append(1e2) # Peak intensity
				guess.append(0.0) # background
			
			fitfct  = lambda a: self.signals[region,col] - pearson7(self.eloss[region],a)
			res     = optimize.leastsq(fitfct,guess)
			yres    = pearson7(self.eloss,res[0])

			plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
			plt.legend(('data','pearson fit','data - pearson'))
			plt.draw()

			if stoploop:
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
			if overwrite: 
				self.signals[:,col] = self.signals[:,col] - yres
		plt.ioff()

	def removeconst(self,whichq,emin,emax,ewindow=100.0):
		"""
		fits a constant as background in the range emin-emax and 
		saves the constant in self.back and the background subtracted self.signals in self.sqw
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		plt.ion()
		for col in columns:
			inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
			res  = np.polyfit(self.eloss[inds],np.transpose(self.signals[inds,col]), 0)
			yres = np.polyval(res, self.eloss)

			plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
			plt.legend(('signal','constant fit','signal - constant'))
			plt.title('Hit [enter] in the python shell to continue')			
			plt.xlabel('energy loss [eV]')
			plt.ylabel('signal [a.u.]')
			plt.xlim(emin-ewindow,emax+ewindow)
			plt.autoscale(enable=True, axis='y', tight=False)
			
			plt.draw()

			self.background[:,col] = yres
			self.sqw[:,col]        = self.signals[:,col] - yres
			_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def removelinear(self,whichq,emin,emax,ewindow=100.0,stoploop=True):
		"""
		fits a linear function as background in the range emin-emax and 
		saves the linear in self.back and the background subtracted self.signals in self.sqw
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		plt.ion()
		for col in columns:
			inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
			res  = np.polyfit(self.eloss[inds],np.transpose(self.signals[inds,col]), 1)
			yres = np.polyval(res, self.eloss)
			
			plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
			plt.legend(('signal','linear fit','signal - linear'))
			plt.title('Hit [enter] in the python shell to continue')			
			plt.xlabel('energy loss [eV]')
			plt.ylabel('signal [a.u.]')
			plt.grid(False)
			plt.xlim(emin-ewindow,emax+ewindow) 
			plt.autoscale(enable=True, axis='y')
			plt.draw()

			self.background[:,col] = yres
			self.sqw[:,col]        = self.signals[:,col] - yres
			if stoploop:
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def removepoly(self,whichq,emin,emax,polyorder=2.0,ewindow=100.0):
		"""
		fits a polynomial of order "polyorder" (default is quadratic) as background 
		in the range emin-emax and saves the polynomial in self.back and the background 
		subtracted self.signals in self.sqw
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		plt.ion()
		for col in columns:
			inds = np.where(np.logical_and(self.eloss >= emin,self.eloss <= emax))
			res  = np.polyfit(self.eloss[inds],np.transpose(self.signals[inds,col]), polyorder)
			yres = np.polyval(res, self.eloss)
			
			plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
			plt.legend(('signal','poly fit','signal - poly'))
			plt.title('Hit [enter] in the python shell to continue')			
			plt.xlabel('energy loss [eV]')
			plt.ylabel('signal [a.u.]')
			plt.xlim(emin-ewindow,emax+ewindow) 
			plt.autoscale(enable=True, axis='y')
			plt.draw()

			self.background[:,col] = yres
			self.sqw[:,col]        = self.signals[:,col] - yres
			_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def removeconstpcore(self,whichq,constregion,coreregion,weights=[5,1],guess=[0.0, 1.0],ewindow=100.0,stoploop=True):
		"""
		fit a const to the preedge and scale data to postedge
		matches the theory profiles
		fminconv: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq		

		region1 = np.where(np.logical_and(self.eloss >= constregion[0], self.eloss <= constregion[1]))
		region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
		region  = np.append(region1*weights[0],region2*weights[1])

		plt.ion()
		for col in columns:
			# first scale data to same area as core profile in region2
			corenorm = np.trapz(self.C[region2,col],self.eloss[region2])
			self.signals[:,col] = self.signals[:,col]/np.trapz(self.signals[region2,col],self.eloss[region2])*corenorm

			#fitfct  = lambda a: a[1]*self.signals[region,col] - (np.polyval([a[0]],self.eloss[region])+self.C[region,col])
			#res = optimize.leastsq(fitfct,guess)[0]
			#yres = np.polyval([res[0]],self.eloss)

			c1 = lambda a: -np.sum((a[1]*self.signals[region2,col] - (np.polyval([a[0]],self.eloss[region2])+self.C[region2,col]))**2.0) # post edge should oscillate around HF core profile
			
			fitfct  = lambda a: np.sum( (a[1]*self.signals[region1,col] - (np.polyval([a[0]],self.eloss[region1])+self.C[region1,col])) )
			cons = [c1] #, c2, c3, c4, c5, c6
			res     = optimize.fmin_cobyla(fitfct,guess,cons,maxfun=10000)
			yres = np.polyval([res[0]],self.eloss)

			plt.plot(self.eloss,self.signals[:,col]*res[1],self.eloss,yres+self.C[:,col],self.eloss,self.signals[:,col]*res[1]-yres,self.eloss,self.C[:,col])
			plt.legend(('data','fit','data - constant','core profile'))
			plt.title('Hit [enter] in the python shell to continue')			
			plt.xlabel('energy loss [eV]')
			plt.ylabel('signal [a.u.]')
			plt.xlim(constregion[0]-ewindow,coreregion[1]+ewindow) 
			plt.autoscale(enable=True, axis='y')
			plt.draw()
			# save some data for paper-plot
			#thedata = np.zeros((len(self.eloss),5))
			#thedata[:,0] = self.eloss
			#thedata[:,1] = self.signals[:,col]*res[1]
			#thedata[:,2] = yres+self.C[:,col]
			#thedata[:,3] = self.signals[:,col]*res[1]-yres
			#thedata[:,4] = self.C[:,col]
			#filename = '/home/csahle/Dropbox/tool_paper/figures/analysis/licl_linpcore_background_det' + '%s' % col + '.dat'
			#np.savetxt(filename,thedata)
			self.background[:,col] = yres
			self.sqw[:,col]        = res[1]*self.signals[:,col] - yres
			if stoploop:
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def removelinpcore(self,whichq,linregion,coreregion,weights=[5,1],guess=[0.0, 0.0, 1.0],ewindow=100.0,stoploop=True):
		"""
		fit a linear to the preedge and scale data to postedge
		matches the theory profiles
		fminconv: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq		

		region1 = np.where(np.logical_and(self.eloss >= linregion[0], self.eloss <= linregion[1]))
		region2 = np.where(np.logical_and(self.eloss >= coreregion[0], self.eloss <= coreregion[1]))
		region  = np.append(region1*weights[0],region2*weights[1])

		plt.ion()
		for col in columns:
			# first scale data to same area as core profile in region2
			corenorm = np.trapz(self.C[region2,col],self.eloss[region2])
			self.signals[:,col] = self.signals[:,col]/np.trapz(self.signals[region2,col],self.eloss[region2])*corenorm

			# then try a minimization
			#fitfct  = lambda a: (a[2]*self.signals[region,col] - (np.polyval(a[0:2],self.eloss[region])+self.C[region,col]))
			#res = optimize.leastsq(fitfct,guess)[0]
			#yres = np.polyval(res[0:2],self.eloss)

			c1 = lambda a: -np.sum((a[2]*self.signals[region2,col] - (np.polyval(a[0:2],self.eloss[region2])+self.C[region2,col]))**2.0) # post edge should oscillate around HF core profile
			c2 = lambda a: a[2] # scaling should not be negative
			fitfct  = lambda a: np.sum( (a[2]*self.signals[region1,col] - (np.polyval(a[0:2],self.eloss[region1])+self.C[region1,col])) )
			cons = [c1,c2] #, c2, c3, c4, c5, c6
			res     = optimize.fmin_cobyla(fitfct,guess,cons)
			yres = np.polyval(res[0:2],self.eloss)

			plt.plot(self.eloss,res[2]*self.signals[:,col],self.eloss,yres+self.C[:,col],self.eloss,res[2]*self.signals[:,col]-yres,self.eloss,self.C[:,col])
			plt.legend(('data','fit','data - linear','core profile'))
			plt.title('Hit [enter] in the python shell to continue')			
			plt.xlabel('energy loss [eV]')
			plt.ylabel('signal [a.u.]')
			plt.xlim(linregion[0]-ewindow,coreregion[1]+ewindow) 
			plt.autoscale(enable=True, axis='y')
			plt.draw()

			self.background[:,col] = yres
			self.sqw[:,col]        = res[2]*self.signals[:,col] - yres
			if stoploop:			
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def removepearson(self,whichq,emin,emax,guess=None,stoploop=True):
		"""
		guess values: 
		a[0] = Peak position
		a[1] = FWHM
		a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
		a[3] = Peak intensity
		a[4] = background
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		region = np.where(np.logical_and(self.eloss>=emin,self.eloss<=emax))[0]
		guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

		plt.ion()
		for col in columns:
			if not guess: 
				guess = []
				ind = np.where( self.signals[guessregion,col] == np.max(self.signals[guessregion,col]) )[0][0]
				guess.append(self.eloss[ind]) # max of signal (in range of prenorm from __init__)
				guess.append(guess[0]*2.0) # twice the position of the peason maximum
				guess.append(1000.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
				guess.append(1.0) # Peak intensity
				guess.append(0.0) # background
			
			fitfct  = lambda a: self.signals[region,col] - pearson7(self.eloss[region],a)
			res     = optimize.leastsq(fitfct,guess)
			yres    = pearson7(self.eloss,res[0])

			plt.plot(self.eloss,self.signals[:,col],self.eloss,yres,self.eloss,self.signals[:,col]-yres)
			plt.legend(('data','pearson fit','data - pearson'))
			plt.draw()

			self.background[:,col] = yres
			self.sqw[:,col]        = self.signals[:,col] - yres
			if stoploop:
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def removecoreppearson(self,whichq,pearsonrange,postrange,weights=[2,1],guess=None,stoploop=True):
		"""
		weights must be integers!
		
		guess values: 
		pearson (always zero background):
		a[0] = Peak position
		a[1] = FWHM
		a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
		a[3] = Peak intensity
		linear:
		a[4] = linear slope
		a[5] = linear background/offset
		data: 
		a[6] = scaling factor
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq		

		region1 = np.where(np.logical_and(self.eloss >= pearsonrange[0], self.eloss <= pearsonrange[1]))[0]
		region2 = np.where(np.logical_and(self.eloss >= postrange[0], self.eloss <= postrange[1]))[0]
		region  = np.append(region1*weights[0],region2*weights[1])
		guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

		plt.ion()
		for col in columns:
			if not guess: 
				guess = []
				ind = self.signals[guessregion,col].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
				guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
				guess.append(guess[0]*1.0) # once the position of the peason maximum
				guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
				guess.append(self.signals[guessregion,col][ind]) # Peak intensity
				guess.append(0.0) # linear slope
				guess.append(0.0) # linear background
				guess.append(1.0) # scaling factor for exp. data

			# some sensible boundary conditions for the fit:
			c1 = lambda a: a[1]*np.absolute(2e2  - a[1]) # FWHM should not be bigger than 200 eV and positive
			c2 = lambda a: a[2] # shape should not be negative
			c3 = lambda a: a[3] # peak intensity should not be negative
			c4 = lambda a: np.absolute(5e-1 - a[4]) # slope for linear background should be small
			c5 = lambda a: a[3] - a[5] # offset for linear should be smaller than maximum of pearson
			c6 = lambda a: a[6]*np.absolute(1e10 - a[6]) # scaling factor for the data should not be negative
			
			fitfct  = lambda a: np.sum( (a[6]*self.signals[region,col] - pearson7_zeroback(self.eloss[region],a[0:4]) - np.polyval(a[4:6],self.eloss[region]) - self.C[region,col])**8.0 )
			cons = [c1] #, c2, c3, c4, c5, c6
			res     = optimize.fmin_cobyla(fitfct,guess,cons)
			print res
			yres    = pearson7_zeroback(self.eloss,res[0:4]) + np.polyval(res[4:6],self.eloss) + self.C[:,col]

			plt.plot(self.eloss,self.signals[:,col]*res[6],self.eloss,yres,self.eloss,self.signals[:,col]*res[6]-yres,self.eloss,self.C[:,col])
			plt.legend(('scaled data','pearson + linear + core','data - (pearson + linear +core)','core'))
			plt.draw()

			self.background[:,col] = yres
			self.sqw[:,col]        = self.signals[:,col] - yres
			if stoploop:
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def removecoreppearson2(self,whichq,pearsonrange,postrange,guess=None,stoploop=True):
		"""
		guess values: 
		pearson (always zero background):
		a[0] = Peak position
		a[1] = FWHM
		a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
		a[3] = Peak intensity
		a[4] = pearson offset
		linear:
		a[5] = linear slope
		a[6] = linear background/offset
		data: 
		a[7] = scaling factor
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq		

		region1 = np.where(np.logical_and(self.eloss >= pearsonrange[0], self.eloss <= pearsonrange[1]))[0]
		region2 = np.where(np.logical_and(self.eloss >= postrange[0], self.eloss <= postrange[1]))[0]
		
		guessregion = np.where(np.logical_and(self.eloss>=self.prenormrange[0],self.eloss<=self.prenormrange[1]))[0]

		plt.ion()
		for col in columns:
			if not guess: 
				guess = []
				ind = self.signals[guessregion,col].argmax(axis=0) # find index of maximum of signal in "prenormrange" (defalt [5,inf])
				guess.append(self.eloss[guessregion][ind]) # max of signal (in range of prenorm from __init__)
				guess.append(guess[0]*1.0) # once the position of the peason maximum
				guess.append(1.0) # pearson shape, 1 = Lorentzian, infinite = Gaussian
				guess.append(self.signals[guessregion,col][ind]) # Peak intensity
				guess.append(0.0) # pearson offset
				guess.append(1.0) # scaling factor for exp. data
			#  - np.polyval(a[5:7],self.eloss[region1])
			# boundary conditions for the fit: let the post-edge region oscilate around the HF core profile
			c1 = lambda a: np.sum( (a[5]*self.signals[region2,col] - pearson7_zeroback(self.eloss[region2],a[0:5]) - self.C[region2,col])**2.0 )
			
			fitfct  = lambda a: np.sum( (a[5]*self.signals[region1,col] - pearson7_zeroback(self.eloss[region1],a[0:5]) - self.C[region1,col])**2.0 )
			cons = [c1] #, c2, c3, c4, c5, c6
			res     = optimize.fmin_cobyla(fitfct,guess,cons)

			yres    = pearson7_zeroback(self.eloss,res[0:5])

			plt.plot(self.eloss,self.signals[:,col]*res[5],self.eloss,yres,self.eloss,self.signals[:,col]*res[5]-yres ,self.eloss,self.C[:,col])

			plt.legend(('scaled data','pearson','data - pearson','core'))
			plt.draw()

			self.background[:,col] = yres
			self.sqw[:,col]        = self.signals[:,col] - yres
			if stoploop:
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def remquickval(self,whichq,corefitrange,interpolrange,convwidth,stoploop=True):
		"""
		quick and dirty way of valence profile extraction from a single spectrum. 
		works if the edge rides on the tail of the valence profile at high q. the HF core
		profile is fitted to the spectrum in the 'corefitrange' and the resulting valence 
		profile is cut out and interpolated over in the 'interpolrange'. finally the valence 
		profile is smoothed by convolution with a gaussian of FWHM 'convwidth' and subtracted
		from the original data such that the resulting S(q,w) oscillates around the HF
		core profile.
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		fitrange = np.where(np.logical_and(self.eloss >= corefitrange[0], self.eloss <= corefitrange[1]))[0]

		interprange1 = np.where(self.eloss<=interpolrange[0])[0]
		interprange2 = np.where(self.eloss>=interpolrange[1])[0]
		interprange  = np.append(interprange1,interprange2)

		plt.ion()
		for col in columns:
			fitfct   = lambda a: np.sum((self.signals[fitrange,col] - a*self.C[fitrange,col])**2.0)
			constr   = lambda a: a # scaling factor for the HF core profile should not be negative

			res = optimize.fmin_cobyla(fitfct,[1.0],cons=constr)[0]

			# subtract the HF core compton profile and interpolate through the edge
			f       = interpolate.interp1d(self.eloss[interprange],self.signals[interprange,col]-res*self.C[interprange,col], bounds_error=False, fill_value=0.0)
			valdata = f(self.eloss)
			valdata = convg(self.eloss,valdata,convwidth)

			subdata = self.signals[:,col] - valdata

			plt.plot(self.eloss,self.signals[:,col],self.eloss,self.C[:,col]*res,self.eloss,valdata,self.eloss,subdata)
			plt.legend(('data','scaled core compton','estimated valence','extracted data'))			
			plt.draw()

			self.valence[:,col] = valdata
			self.sqw[:,col]     = subdata/res # scale the extracted data back to fit the HF core profile

			if stoploop:
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()    # close the figure to show the next one
		plt.ioff()

	def extractval(self,whichq,mirror=False,linrange1=None,linrange2=None):
		"""
		extracts a valence profile from q-value(s) given in whichq by first fitting 
		the core HF profile to the data at places linrange1 and linrange2 (one, two, 
		or no ranges can be given), then subtracting the HF profile from the data. the
		resulting valence profile in the near-edge region can be replaced by a pearson 
		function (default) or by mirroring the negative side of the valence profile 
		(mirror=True). if mirror is set to False, also the asymmetry is fitted.
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		# set pz scale of first given q-val as 'master' grid
		absolutepz = e2pz(self.eloss/1e3+self.E0,self.E0,self.tth[columns[0]])[0]
		
		plt.cla()
		plt.ion()
		for col in columns: 
			# set the pz scale for each q
			pz = e2pz(self.eloss/1e3+self.E0,self.E0,self.tth[col])[0]

			if linrange1 and linrange2:
				range1   = np.where(np.logical_and(self.eloss>=linrange1[0],self.eloss<=linrange1[1]))[0]
				range2   = np.where(np.logical_and(self.eloss>=linrange2[0],self.eloss<=linrange2[1]))[0]
				linrange = np.append(range1,range2)
			elif linrange1:
				linrange = np.where(np.logical_and(self.eloss>=linrange1[0],self.eloss<=linrange1[1]))[0]
			else: 
				linrange = np.where(0.1*self.C[:,col] > self.V[:,col])[0]

			# simple minimization to subtract a linear from the data to fit ontop of the Compton profile:
			fitfct = lambda a: (self.signals[linrange,col] - np.polyval([a[0],a[1]],self.eloss[linrange]) ) - self.J[linrange,col]
			res    = optimize.leastsq(fitfct,[0.0,0.0])[0]

			# raw valence (when later extracted from the data, also a linear shoud be accounted for)
			val = self.signals[:,col] - np.polyval(res,self.eloss) - self.C[:,col]

			if mirror: # just replace the edgepart of the valence profile by the other half of the profile
				mirrorval = np.append(val[pz<=0.0],np.flipud(val[pz<=0]))
				mirrorpz  = np.append(pz[pz<=0],np.flipud(pz[pz<=0]*-1))
				order = np.argsort(mirrorpz)

				f = interpolate.interp1d(mirrorpz[order],mirrorval[order],bounds_error=False, fill_value=0.0)
				extractedval = f(absolutepz)
				
				plt.plot(absolutepz,val,absolutepz,extractedval)
				plt.legend(['exp. S(q,w) - HF core profile','mirrored extracted valence profile'])
				plt.xlabel('pz [a.u.]')
				plt.ylabel('S(q,w) [1/eV]')
				plt.draw()
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
				plt.close()

				self.valence[:,col] = extractedval
				self.valasymmetry = np.zeros_like(self.valence)

			else: # fit pearson to replace near edge part
				print ('select a point above which the valence profile should be replaced by a Pearson function')
				plt.plot(absolutepz,val)
				plt.ylim(np.amin(val)-np.amin(val)*1.1,val[absolutepz.flat[np.abs(absolutepz - 0.0).argmin()]]*1.5) 
				xyval       = pylab.ginput(1)[0]
				edgeregion  = np.where(pz < xyval[0])[0]
				start_param = [pz[val==np.amax(val[edgeregion])][0], 4.0, 1.0, np.amax(val[edgeregion]), 0.0 ]
				fitfct      = lambda a: val[edgeregion] - pearson7(pz[edgeregion],a)
				param = optimize.leastsq(fitfct,start_param)[0]
				# param   = optimize.curve_fit(pearson7_forcurvefit, pz[theregion], val[theregion],p0=start_param)[0]
				pearson = pearson7(pz,param)
				val[pz>xyval[0]] = pearson[pz>xyval[0]]
				extractedval = val;
				plt.close()

				# try fitting the valence asymmetry
				print 'trying to extract valence asymmetry!'
				pzp = absolutepz[absolutepz >= 0.0]
				pzm = absolutepz[absolutepz < 0.0]
				jvalp  = extractedval[absolutepz >=0.0 ]
				f = interpolate.UnivariateSpline(-pzm,extractedval[absolutepz<0.0])
				jvalm  = f(pzp)
				fitfct = lambda a: jvalp-jvalm - a[0]*(np.tanh(pzp/a[1])*np.exp(-(pzp/np.absolute(a[2]))**4.0))
				res    = optimize.leastsq(fitfct,[0.0,1.0,1.0])[0]
				asym   = (res[0]*(np.tanh(pz/res[1])*np.exp(-(pz/np.absolute(res[2]))**4.0)))/2.0
				plt.plot(absolutepz,extractedval,absolutepz,asym,absolutepz,extractedval+asym)
				plt.legend(['extracted valence profile','fitted valence asymmetry','asymmetry corrected valence profile'])
				plt.draw()
				self.valencepz[:,col]      = extractedval-asym
				self.valasymmetrypz[:,col] = asym

		self.pzscale = absolutepz
		plt.ioff()

	def extractval_test(self,whichq,mirror=False,linrange1=None,linrange2=None):
		"""
		extracts a valence profile from q-value(s) given in whichq by first fitting 
		the core HF profile to the data at places linrange1 and linrange2 (one, two, 
		or no ranges can be given), then subtracting the HF profile from the data. the
		resulting valence profile in the near-edge region can be replaced by a pearson 
		function (default) or by mirroring the negative side of the valence profile 
		(mirror=True). if mirror is set to False, also the asymmetry is fitted.
		"""
		
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		plt.cla()
		plt.ion()

		for col in columns: 
			# set the pz scale for each q
			pz_dmy = e2pz(self.eloss/1e3+self.E0,self.E0,self.tth[col])[0]

			# find the regions in which to fit a linear function
			if linrange1 and linrange2:
				range1   = np.where(np.logical_and(self.eloss>=linrange1[0],self.eloss<=linrange1[1]))[0]
				range2   = np.where(np.logical_and(self.eloss>=linrange2[0],self.eloss<=linrange2[1]))[0]
				linrange = np.append(range1,range2)
			elif linrange1:
				linrange = np.where(np.logical_and(self.eloss>=linrange1[0],self.eloss<=linrange1[1]))[0]
			else: 
				linrange = np.where(0.1*self.C[:,col] > self.V[:,col])[0]

			# simple minimization to subtract a linear from the data to fit ontop of the HF Compton profile:
			fitfct = lambda a: (self.signals[linrange,col] - np.polyval([a[0],a[1]],self.eloss[linrange]) ) - a[2]*self.J[linrange,col]
			res    = optimize.leastsq(fitfct,[0.0,0.0,1.0])[0]

			# raw valence (when later extracted from the data, also a linear shoud be accounted for)
			val = self.signals[:,col] - np.polyval(res[0:2],self.eloss) - res[2]*self.C[:,col]
			
			#plt.plot(self.eloss,self.signals[:,col],self.eloss,self.C[:,col]*res[2]+np.polyval(res[0:2],self.eloss))
			if mirror: # just replace the edgepart of the valence profile by the other half of the profile
				mirrorval = np.append(val[pz_dmy<=0.0],np.flipud(val[pz_dmy<=0]))
				mirrorpz  = np.append(pz_dmy[pz_dmy<=0],np.flipud(pz_dmy[pz_dmy<=0]*-1))
				order = np.argsort(mirrorpz)

				f = interpolate.interp1d(mirrorpz[order],mirrorval[order],bounds_error=False, fill_value=0.0)
				extractedval = f(pz_dmy)
				
				plt.plot(pz_dmy,val,pz_dmy,extractedval)
				plt.legend(['exp. S(q,w) - HF core profile','mirrored extracted valence profile'])
				plt.xlabel('pz [a.u.]')
				plt.ylabel('S(q,w) [1/eV]')
				plt.draw()
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
				plt.close()

				self.valence[:,col] = extractedval
				self.valasymmetry = np.zeros_like(self.valence)

			else: # fit pearson to replace near edge part
				print ('select a point above which the valence profile should be replaced by a Pearson function')
				plt.plot(pz_dmy,val)
				plt.ylim((np.amin(val[pz_dmy<2.0])*0.9,np.amax(val[pz_dmy<2.0])*1.1)) # np.amin(val)-np.amin(val)*1.1,val[pz_dmy.flat[np.abs(pz_dmy - 0.0).argmin()]]*1.5) 
				xyval       = pylab.ginput(1)[0]
				edgeregion  = np.where(pz_dmy < xyval[0])[0]
				start_param = [pz_dmy[val==np.amax(val[edgeregion])][0], 4.0, 1.0, np.amax(val[edgeregion]), 0.0 ]
				fitfct      = lambda a: val[edgeregion] - pearson7(pz_dmy[edgeregion],a)
				param = optimize.leastsq(fitfct,start_param)[0]
				# param   = optimize.curve_fit(pearson7_forcurvefit, pz_dmy[theregion], val[theregion],p0=start_param)[0]
				pearson = pearson7(pz_dmy,param)
				val[pz_dmy>xyval[0]] = pearson[pz_dmy>xyval[0]]
				extractedval = val;
				plt.close()

				# try fitting the valence asymmetry
				print 'trying to extract valence asymmetry!'
				pzp = pz_dmy[pz_dmy >= 0.0]
				pzm = pz_dmy[pz_dmy < 0.0]
				jvalp  = extractedval[pz_dmy >=0.0 ]
				f = interpolate.interp1d(-pzm,extractedval[pz_dmy<0.0],bounds_error=False, fill_value=0.0)
				jvalm  = f(pzp)
				fitfct = lambda a: jvalp-jvalm - a[0]*(np.tanh(pzp/a[1])*np.exp(-(pzp/np.absolute(a[2]))**4.0))
				res    = optimize.leastsq(fitfct,[0.0,1.0,1.0])[0]
				print res
				asym   = (res[0]*(np.tanh(pz_dmy/res[1])*np.exp(-(pz_dmy/np.absolute(res[2]))**4.0)))/2.0
				plt.plot(pz_dmy,extractedval,pz_dmy,asym,pz_dmy,extractedval+asym)
				#plt.plot(pz_dmy[pz_dmy<0],extractedval[pz_dmy<0]+asym[pz_dmy<0],-pz_dmy[pz_dmy>=0],extractedval[pz_dmy>=0]+asym[pz_dmy>=0])
				plt.legend(['extracted valence profile','fitted valence asymmetry','asymmetry corrected valence profile'])
				plt.draw()
				#self.valence[:,col]      = extractedval-asym
				#self.valasymmetry[:,col] = asym

				#print pz_dmy[0], pz_dmy[-1], absolutepz[0],absolutepz[-1]
				f = interpolate.interp1d(np.flipud(pz_dmy),np.flipud(extractedval-asym), kind='cubic',bounds_error=False, fill_value=0.0)
				#f = interpolate.UnivariateSpline(pz_dmy,extractedval-asym,k=3,s=10)
				absval = f(np.flipud(self.pzscale))
				self.valencepz[:,col] = np.flipud(absval)
				f = interpolate.interp1d(np.flipud(pz_dmy),np.flipud(asym), kind='cubic',bounds_error=False, fill_value=0.0)
				absasym = f(np.flipud(self.pzscale))
				self.valasymmetrypz[:,col] = np.flipud(absasym)
				#print absval
				#plt.cla()				
				#plt.plot(absolutepz,absval)
				#plt.cla()
				#plt.plot(absolutepz,np.interp(absolutepz,pz_dmy,extractedval-asym),pz_dmy,extractedval-asym)
				#break

			#self.pzscale[:,col] = absolutepz #pz_dmy
		plt.ioff()

	def getallvalprof(self,whichq,smoothgval=0.0,stoploop=True):
		"""
		takes the extracted valence profile extracted from q-value whichq
		and transforms them onto the other q-values 
		whichq     = column from which the valence profile was extracted
		smoothgval = FWHM used for gaussian smoothing (default is 0.0, i.e. no smoothing)
		stoploop   = boolean, plots each result if set to True
		"""
		newenergy  = np.zeros((len(self.pzscale),len(self.tth)))
		newvalence = np.zeros((len(self.eloss),len(self.tth)))
		newasym    = np.zeros((len(self.eloss),len(self.tth)))
		plt.ion()
		for n in range(len(self.tth)):
			newenergy[:,n]  = (pz2e1(self.E0,self.pzscale,self.tth[n]) - self.E0)*1e3 # each energy scale in [eV]
			f               = interpolate.interp1d(newenergy[:,n],self.valencepz[:,whichq],bounds_error=False, fill_value=0.0)
			newvalence[:,n] = f(self.eloss)/self.qvals[:,n]
			f               = interpolate.interp1d(newenergy[:,n],self.valasymmetrypz[:,whichq],bounds_error=False, fill_value=0.0)
			newasym[:,n]    = f(self.eloss)/self.qvals[:,n]

			plt.plot(self.eloss,newvalence[:,n],self.eloss,newasym[:,n])
			plt.draw()
			if stoploop:
				_ = raw_input("Press [enter] to continue.") # wait for input from the user
			plt.close()

		if smoothgval > 0.0:
			for n in range(len(self.tth)):
				self.valence[:,n] = convg(self.eloss,newvalence[:,n],smoothgval) + newasym[:,n]
			self.valasymmetry = newasym
		else:
			self.valence = newvalence + newasym
			self.valasymmetry = newasym
		plt.ioff()

	def remvalenceprof(self,whichq,eoffset=0.0):
		"""
		removes extracted valence profile from q-values given in list whichq by fitting the scaled
		valence profile plus a linear function plus the HF core compton profile to the data
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		inds = np.where(self.eloss>=self.prenormrange[0])[0]

		plt.ion()
		for col in columns:
			f = interpolate.interp1d(self.eloss+eoffset,self.valence[:,col],bounds_error=False, fill_value=0.0)
			#self.valence[:,col]  = f(self.eloss)
			thevalence = f(self.eloss)
			fitfct = lambda a: self.signals[inds,col]-self.C[inds,col]-a[0]*thevalence[inds]-np.polyval([a[1],a[2]],self.eloss[inds])
			res    = optimize.leastsq(fitfct,[6.0,0.0,0.0])[0]
			plt.plot(self.eloss,self.signals[:,col])
			plt.plot(self.eloss,self.C[:,col]+res[0]*thevalence)
			plt.plot(self.eloss,np.polyval(res[1:3],self.eloss))
			plt.plot(self.eloss,self.signals[:,col]-res[0]*thevalence-np.polyval([res[1],res[2]],self.eloss),self.eloss,self.C[:,col])
			plt.draw()
			self.sqw = self.signals[:,col]-res[0]*thevalence-np.polyval([res[1],res[2]],self.eloss)
			self.valence[:,col] = thevalence
		plt.ioff()

	def remvalenceprof_test(self,whichq,eoffset=0.0):
		"""
		removes extracted valence profile from q-values given in list whichq by fitting the scaled
		valence profile plus a linear function plus the HF core compton profile to the data
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else:
			columns = whichq

		inds = np.where(self.eloss>=self.prenormrange[0])[0]

		plt.ion()
		for col in columns:
			f = interpolate.interp1d(self.eloss+eoffset,self.valence[:,col],bounds_error=False, fill_value=0.0)
			#self.valence[:,col]  = f(self.eloss)
			thevalence = f(self.eloss)

			fitfct = lambda a: self.signals[inds,col]-self.C[inds,col]-a[0]*np.interp(self.eloss,self.eloss+a[1],thevalence)[inds]-np.polyval([a[2],a[3]],self.eloss[inds])
			res    = optimize.leastsq(fitfct,[1.0,0.0,0.0,0.0])[0]
			plt.plot(self.eloss,self.signals[:,col])
			plt.plot(self.eloss,self.C[:,col]+res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)+np.polyval(res[2:4],self.eloss))
			plt.plot(self.eloss,np.polyval(res[2:4],self.eloss))
			plt.plot(self.eloss,self.signals[:,col]-res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)-np.polyval([res[2],res[3]],self.eloss),self.eloss,self.C[:,col])
			plt.draw()
			self.sqw[:,col] = self.signals[:,col]-res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)-np.polyval([res[2],res[3]],self.eloss)
			# save some ascii files for paper
			thedata = np.zeros((len(self.eloss),6))
			thedata[:,0] = self.eloss
			thedata[:,1] = self.signals[:,col]
			thedata[:,2] = res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)+np.polyval(res[2:4],self.eloss)
			thedata[:,3] = np.polyval(res[2:4],self.eloss)
			thedata[:,4] = self.signals[:,col]-res[0]*np.interp(self.eloss,self.eloss+res[1],thevalence)-np.polyval([res[2],res[3]],self.eloss)
			thedata[:,5] = self.C[:,col]
			filename = '/home/csahle/Dropbox/tool_paper/figures/analysis/val_extraction_det' + '%s' % col + '.dat'
			np.savetxt(filename,thedata) 
		plt.ioff()

	def averageqs(self,whichq,errorweighing=True):
		"""
		averages S(q,w) over the q-values given
		whichq        = list of q-values over which to average (index starts at zero)
		errorweighing = boolean, weighs sum by errors if set to True
		"""
		if not isinstance(whichq,list):
			columns = []
			columns.append(whichq)
		else: 
			columns = whichq

		# build the matricies
		av    = np.zeros((len(self.eloss),len(whichq)))
		averr = np.zeros((len(self.eloss),len(whichq)))
		for n in range(len(columns)):
			# find data points with error = 0.0 and replace by 1.0
			inds = np.where(self.errors[:,columns[n]] == 0.0)[0]
			for ind in inds:
				self.errors[ind,columns[n]] = 1.0
			# arrange the desired columns into a matrix
			av[:,n]    = self.sqw[:,columns[n]]
			averr[:,n] = self.errors[:,columns[n]]
		# sum things up
		if errorweighing:
			self.sqwav    = np.sum( av/averr**2.0 ,axis=1)/( np.sum(1.0/averr**2.0,axis=1))
			self.sqwaverr = np.sqrt( 1.0/np.sum(1.0/(averr)**2.0,axis=1) )

		else: 
			self.sqwav    = np.sum(av,axis=1)
			self.sqwaverr = np.sqrt(np.absolute(self.sqwav)) # check this again




