import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as scp
import scipy.optimize
import copy

import numexpr
numexpr.set_num_threads(1)

import time
from scipy import interpolate

import numpy.lib.recfunctions

from scripts.Injector import Injector
from scripts.PDF import PDF
from scripts.RandomTools import RandomTools

import resource

class LLh(PDF, Injector):
	"""A class for the Log Likelihhod.
	"""
	
	def __init__(self, **kwargs):
		#Bins for energy (Tev?)
		self.EnergyBins = np.linspace(1., 10., 40+1)
		#Bins for sin declination (not evenly spaced)
		self.sinDecBins = np.unique(np.concatenate([
								np.linspace(-1., -0.9, 2+1),
								np.linspace(-0.9, -0.2, 8+1),
								np.linspace(-0.2, 0.2, 15+1),
								np.linspace(0.2, 0.9, 12+1),
								np.linspace(0.9, 1., 2+1),
								]))
		#Provides grid values for spectral index (gamma)
		self.GammaGrid = np.linspace(0, 5., 20+1)
		#Sets precision
		self.precision = .1
		
		#Sets weights as nan for some crazy reason
		self.SeasonWeight = np.nan
		self.AcceptanceWeight = np.nan
		
		#Sets "_..." to False
		self._ReturnInjectorNExp = False
		
		self._g1 = np.nan
		self._w_cache = np.nan
		self.w_cache = dict()
		self.spline_cache = dict()
		self.w_cache_BG = dict()
		self.w_cache_Sig = dict()
		
		#Produces a set (i.e no duplicates) of datapoints for gamma
		#This is best on 33 points between 0.9 and 4.1, with each modified by _around(i)
		#Useful for different recisions, where rounding errors might lead to duplicates in set
		self.GammaSupportPoints = set([self._around(i) for i in np.linspace(0.9, 4.1, 30+3)])
		
		self.N_all = np.nan
		self.n_select = np.nan
		self._exp = np.nan
		
		#Sets start and end date for data taking (typo irrelevant)
		self.DataStart = kwargs['StartDataTakingMJD']
		self.DataEnd = kwargs['EndDataTankingMJD']
		
		#Sets the livetime, and the total Season length
		self.Livetime = kwargs['Livetime']
		self.SeasonTimeSpan = self.DataEnd - self.DataStart
		
		#Loads "Experimental" and "Monte Carlo" data sets (stored in .npy format)
		self.exp = np.load(kwargs['ExpPath'])
		self.mc = np.load(kwargs['MCPath'])
		
		#Does something...
		self._sources2 = np.load(kwargs['SourcePath'])
		self._sources = np.lib.recfunctions.append_fields(self._sources2, 'TimeNorm', data=np.ones_like(self._sources2['ra']))
		del(self._sources2)
		
		try:
			DecBins = np.load(kwargs['AcceptanceWeightPath2']+'_bins_dec.npy')
			GammaBins = np.load(kwargs['AcceptanceWeightPath2']+'_bins_gamma.npy')
			values = np.load(kwargs['AcceptanceWeightPath2']+'_values.npy')
			f = interpolate.interp2d(DecBins, GammaBins, values, kind='linear')
			self.AcceptanceFitFunc = f
			del(f)
		except:
			print('No Acceptance Files')
		
		self.UseEnergy = kwargs['UseEnergy']
		self.FitGamma = kwargs['FitGamma']
		self.FixedGamma = kwargs['FixedGamma']
		
		self.UseBox = kwargs['UseBox']
	  
		self.FitWeights = kwargs['FitWeights']
		
		self.UseTime = kwargs['UseTime']
		self.TimeModel = kwargs['TimeModel']
		if self.UseTime:
			self.DiscDelay = 10.
			self.sources['discoverydate_mjd'] = self.sources['discoverydate_mjd'] - self.DiscDelay
			
			if self.TimeModel=='Box':
				self.TimeBoxLenght = kwargs['TimeBoxLenght']
				
			if self.TimeModel=='BoxPre':
				self.TimeBoxLenght = kwargs['TimeBoxLenght']

			if self.TimeModel=='Decay':
				self.Model_tpp = kwargs['Model_tpp']
				self.DecayModelLenght = kwargs['DecayModelLenght']
				self.Model_Length = self.DecayModelLenght
		
		self.SmearInjection = kwargs['SmearInjection']
		self.MissTiming = kwargs['MissTiming']
			   
		self.InjectionBandMask = dict()
		
		if 'InjectionGamma' in kwargs.keys():
			self.InjectionGamma = kwargs['InjectionGamma']
		else:
			self.InjectionGamma = 2.
			self.WeightsInject = self.GetWeights(self.mc, self.InjectionGamma, )
		print 'Injection Spectrum', self.InjectionGamma
	
		self.WeightCache = np.nan
		
		self.Unblind = False
		
		if self.Unblind==False:
			np.random.shuffle(self.exp['timeMJD'])  
		else:
			print('WARNING: Running in Unblinding Mode')
			
		
		self.BootStrap = False
		try:
			self.BootStrap = kwargs['BootStrap']
		except:
			pass
		
#=====================================================================================
#Running Part
#=====================================================================================    
	@property
	def sources(self, ):
		return self._sources
	
	
	def PrepareFakeDataSetAndEvalautePDF(self, k, ):
		self.w_cache = dict()
		self.w_cache_Sig = dict()
		self._ev = np.nan
		self._ev_B = np.nan
		self._ev_S = np.nan
		self.EnergyWeightCache = np.nan
		self.SoB = np.nan
			
		if self.UseBox == True:
			self.w_cache_BG = dict()
			self._ev = self.SelectEventsInBand(self.sources[0], self.exp)
			self._ev = self.scramble_exp_data(self._ev)
			self._ev = self.SelectEventsInBox(self.sources[0], self._ev, )
		if self.UseBox == False:
			self._ev = self.scramble_exp_data(self.exp)   
		if self.UseEnergy == True:
			self.GenerateBGWeightDictForAllGamma(self._ev)
							
		self.sig_events = self.generate_sig_events(self.sources, self.mc, k, )
		if self.UseEnergy == True:
			self.GenerateSigWeightDictForAllGamma(self.sig_events)
				
		self._ev = np.concatenate( (self._ev, self.sig_events) )
		
		self.N_all = len(self.exp)+len(self.sig_events)
		self.N = len(self._ev)
		
		self.EvaluateB()
		self.EvaluateS()
		del(self._ev)
		self.SoB = self.SoB / self._ev_B
		del(self._ev_B)
	
		if self.UseEnergy == True:
			if self.FitGamma == False:
				gamma = self.FixedGamma
				self.EnergyWeightCache = np.append(self.weightFAST(gamma, self.w_cache_BG),
												   self.weightFAST(gamma,self.w_cache_Sig))
	
	def InitEverythingForMultipleTrials(self, ):
		"""Initialises 
		"""
		#Produces scrambled experimental data set
		self._ev = self.exp
		self._ev = self.scramble_exp_data(self._ev)
		#Finds a spatial PDF for the background, based on the experimental Sin Declination distribution		
		bckg_spline_space = self.create_space_BG_pdf(self._ev)
		self.bckg_spline_space = bckg_spline_space
		
		#Assigns a weight to each source, equal to 1/(r^2) for distance r
		self.sources['weight_distance'] = self.sources['distance']**(-2.)
		
		#If accounting for energy, produces Energy PDFs
		if self.UseEnergy == True:
			print('Inizilizing Energy PDFs')
			self.GenerateSplineDictForAllGamma(self.exp, self.mc)
			self.GenerateBGWeightDictForAllGamma(self._ev)    
	
		if self.UseTime == True:
			self.ComputeSourceWeightsTime()
			self.InitRandomGeneratorPDF()   
	
	def ProduceLLhFunction(self, ):
		f = lambda x: - self.TestStatNewLLhFitWeights(x, )

		return f

#=====================================================================================
#End Running Part
#=====================================================================================      
			  
	def TestStatNewLLh(self, params, ):
		N_all = self.N_all
		N = self.N       
		SoB = self.SoB
		
		if self.UseEnergy==True:
			n = params[:-1]
			gamma = params[-1]
			w = np.append(self.weightFAST(gamma, self.w_cache_BG),
						 (self.weightFAST(gamma, self.w_cache_Sig)))
		if self.UseEnergy==False:
			n = params[:]
			w = 1.
		
		n = n * self.SeasonWeight
		y = w * np.sum(self.sources['weight'][:,None] * SoB, axis=0)
		y = (y-1.) / N_all
		alpha = n * y
		funval = np.log1p(alpha)
		funval = funval.sum() + (N_all - N) * np.log1p(-n / N_all)       
		return 2.* funval
	
	
	def TestStatNewLLhFast(self, params, ):
		N_all = self.N_all
		N = self.N       
		SoB = self.SoB
		
		if self.UseEnergy==True:
			n = params[:-1]
			gamma = params[-1]
			w = np.append(self.weightFAST(gamma, self.w_cache_BG),
						 (self.weightFAST(gamma, self.w_cache_Sig)))
		if self.UseEnergy==False:
			n = params[:]
			w = 1.
		
		n = n * self.SeasonWeight
		b = self.sources['weight'][:,None]
		y = np.sum(numexpr.evaluate('(b * SoB)'), axis=0)        
		funval = numexpr.evaluate('log1p(n * (((w * y)-1.) / N_all))')
		funval = funval.sum() + (N_all - N) * np.log1p(-n / N_all)
		return 2.* funval
	  
	
	def TestStatNewLLhFitWeights(self, params, ):
		N_all = self.N_all
		N = self.N       
		SoB = self.SoB
		
		if self.UseEnergy==True:
			n = params[:-1]
			gamma = params[-1]
			w = np.append(self.weightFAST(gamma, self.w_cache_BG),
						 (self.weightFAST(gamma, self.w_cache_Sig)))
		if self.UseEnergy==False:
			n = params[:]
			w = 1.
		
		n = n * self.SeasonWeight
		b = self.sources['weight'][:,None]

		if self.FitWeights==False:   
			y = np.sum(numexpr.evaluate('(b * SoB)'), axis=0)       
			a = numexpr.evaluate('n * (w*y-1.)')
			del(w)

		if self.FitWeights==True:
			y = numexpr.evaluate('SoB')
			c = np.matrix(numexpr.evaluate('w*SoB-1.'))
			del(w)
			a = np.squeeze(np.asarray(n*c))
			del(c)
		funval = numexpr.evaluate('log1p(a / N_all)')
		
		#Make sure to delete it from memory		
		del(SoB)
		del(y)
		del(a)
		if self.UseBox==True:
			funval = funval.sum() + (N_all - N) * np.log1p(-n / N_all)       
		else:
			funval = funval.sum()
		return 2.* funval