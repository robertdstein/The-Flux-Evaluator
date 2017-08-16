import numpy as np
from scipy.stats import norm
from icecube import astro

import scipy
from scipy.interpolate import UnivariateSpline
from scipy.integrate import romberg
import scipy as scp

import sys

import numexpr
numexpr.set_num_threads(1)

import scipy.interpolate
from scipy import interpolate


class PDF():
	def __init__(self, ):
		pass

#==============================================================================
#Tool Part
#==============================================================================

	def GetWeights(self, mc, gamma, ):
		"""Returns a weights (array) given by
		OneWeights * Energy_weights
		"""
		#Uses numexpr for faster processing
		#lint:disable
		ow = mc['ow']
		trueE = mc['trueE']
		weights = numexpr.evaluate('ow * trueE **(-gamma)')
		#lint:enable
		return weights

	def _around(self, value):
		"""Produces an array in which the prescision of the value
		is rounded to the nearest integer. This is then multiplied
		by the precision, and the array is returned.
		"""
		return np.around(float(value) / self.precision) * self.precision

	def EmptyCache(self, ):
		"""Empties the cache by setting _g1/_w_cache to default value (nan)
		"""
		self._g1 = np.nan
		self._w_cache = np.nan

	def SelectEventsInBand(self, source, data, ):
		"""Selects events in a given sin declination band
		The width of the band is defined as 20 degrees.
		This function, if followed by SelectEventsInBox(),
		can be used to produce a square on the celestial sphere.
		"""
		dec_bandwidth = np.deg2rad(10.)
		min_dec = max(-np.pi / 2., source['dec'] - dec_bandwidth)
		max_dec = min(np.pi / 2., source['dec'] + dec_bandwidth)
		band_mask = ((data['sinDec']) > np.sin(min_dec)) & \
			((data['sinDec']) < np.sin(max_dec))
		data_band = data[band_mask]
		N_all = len(data)
		n_select = len(data_band)
		self.N_all = N_all
		self.n_select = n_select
		return data_band

	def SelectEventsInBox(self, source, data, ):
		"""Selects events in a given ra band.
		In combination with the declination band of SelectEventsInBand(),
		this produces a square on the surface of the celestial sphere.
		The width of the ra band is variable to account for sperical geometry
		"""
		dec_bandwidth = np.deg2rad(10.)
		min_dec = max(-np.pi / 2., source['dec'] - dec_bandwidth)
		max_dec = min(np.pi / 2., source['dec'] + dec_bandwidth)
		cosFact = np.amin(np.cos([min_dec, max_dec]))
		dPhi = np.amin([2. * np.pi, 2. * dec_bandwidth / cosFact])
		ra_dist = np.fabs((data["ra"] - source['ra'] + np.pi)
			% (2. * np.pi) - np.pi)
		mask = ra_dist < dPhi / 2.
		data_box = data[mask]
		return data_box

#==============================================================================
#Energy Part
#==============================================================================

	def weight(self, data, exp, mc, gamma):
		"""Calculates the Signal/Background function.
		Uses Finite Difference Methods to calculate 1st/2nd derivatives.
		Uses a Taylor series to estimate S(gamma).
		"""
		#Sets g (gamma), and the precision/stepsize (dg or h)
		g1 = self._around(gamma)
		dg = self.precision

		if g1 in self.w_cache.keys():
			#If derivatives already calculated, proceed.
			S1 = self.w_cache[g1]["S1"]
			a = self.w_cache[g1]["a"]
			b = self.w_cache[g1]["b"]
		else:
			#Add/subtract a step from g1 (g+h) and (g-h)
			g0 = self._around(g1 - dg)
			g2 = self._around(g1 + dg)

			#Calculate Spline stuff (sig/background 2D histograms)
			#Gives y(g), y(g+h), y(g-h)
			S0 = self.GenerateSplineStuff(exp, mc, gamma=g0).ev(
				data['logE'], data['sinDec'])
			S1 = self.GenerateSplineStuff(exp, mc, gamma=g1).ev(
				data['logE'], data['sinDec'])
			S2 = self.GenerateSplineStuff(exp, mc, gamma=g2).ev(
				data['logE'], data['sinDec'])

			#Calculates Derivatives Using Finite Difference Methods (1st/2nd)
			a = (S0 - 2. * S1 + S2) / (2. * dg ** 2)
			b = (S2 - S0) / (2. * dg)

			#Stores information
			self._g1 = g1
			self._w_cache = np.zeros((len(data),), dtype=[("S1", np.float),
				("a", np.float), ("b", np.float)])
			self._w_cache["S1"] = S1
			self._w_cache["a"] = a
			self._w_cache["b"] = b

			self.w_cache[g1] = dict()
			self.w_cache[g1]["S1"] = S1
			self.w_cache[g1]["a"] = a
			self.w_cache[g1]["b"] = b

		#Returns the Taylor series value for S(gamma) around g1
		val = np.exp(a * (gamma - g1) ** 2 + b * (gamma - g1) + S1)
		return val

	def weightFAST(self, gamma, w_cache):
		"""Quickly estimates the value of Signal/Background for Gamma.
		Uses precalculated values for first and second derivatives.
		Uses a Taylor series to estimate S(gamma).
		"""
		g1 = self._around(gamma)
		dg = self.precision

		g0 = self._around(g1 - dg)
		g2 = self._around(g1 + dg)

		#Uses Numexpr to quickly estimate S(gamma)
		#lint:disable
		S0 = w_cache[g0]
		S1 = w_cache[g1]
		S2 = w_cache[g2]

		val = numexpr.evaluate("exp((S0 - 2.*S1 + S2) / (2. * dg**2) *" +
			"(gamma - g1)**2 + (S2 - S0) / (2. * dg) * (gamma - g1) + S1)")
		#lint:enable
		return val

	def weightFAST2(self, gamma, w_cache):
		"""Calculates the Signal/Background function.
		Using precalculated values for the S(gamma) at different points.
		Uses Finite Difference Methods to calculate 1st/2nd derivatives.
		Uses a Taylor series to estimate S(gamma).
		"""
		g1 = self._around(gamma)
		dg = self.precision

		g0 = self._around(g1 - dg)
		g2 = self._around(g1 + dg)

		S0 = w_cache[g0]
		S1 = w_cache[g1]
		S2 = w_cache[g2]

		a = (S0 - 2. * S1 + S2) / (2. * dg ** 2)
		b = (S2 - S0) / (2. * dg)

		val = np.exp(a * (gamma - g1) ** 2 + b * (gamma - g1) + S1)
		return val

	def GenerateSplineDictForAllGamma(self, exp, mc, ):
		"""Loops over "GenerateSplineStuff" for each Gamma point.
		Returns the interpolated 2D spline of the log(sig/bkg) ratios.
		"""
		for gamma in self.GammaSupportPoints:
			if gamma in self.spline_cache.keys():
				pass
			else:
				self.spline_cache[gamma] = self.GenerateSplineStuff(exp,
					mc, gamma=gamma)

#******************************************************************************************************
	def GenerateBGWeightDictForAllGamma(self, exp):
		"""Loops over each gamma points, and check if the point has been
		added to the Background Weight Cache. If not, loads the associated
		2D spline function. Using the array of experiental data, it converts
		each event to a Log(sig/Bkg) value, using Log(Energy) and sin(Dec)
		with the fitted function.
		"""
		for gamma in self.GammaSupportPoints:
			if gamma in self.w_cache_BG.keys():
				pass
			else:
				g1 = self._around(gamma)
				self.w_cache_BG[gamma] = self.spline_cache[g1].ev(exp['logE'],
					exp['sinDec'])

#*********************************************************************************************************
	def GenerateSigWeightDictForAllGamma(self, sig_events):
		"""Loops over each gamma points, and check if the point has been
		added to the Signal Weight Cache. If not, loads the associated
		2D spline function. Using the array of experiental data, it converts
		each event to a Log(sig/Bkg) value, using Log(Energy) and sin(Dec)
		with the fitted function.
		"""
		for gamma in self.GammaSupportPoints:
			if gamma in self.w_cache_Sig.keys():
				pass
			else:
				g1 = self._around(gamma)
				self.w_cache_Sig[gamma] = self.spline_cache[g1].ev(
					sig_events['logE'], sig_events['sinDec'])

	def GenerateSplineStuff(self, exp, mc, gamma=2.):
		"""Produces normailsed histograms for signal and background.
		Produces a histogram of the ratio of sig/bkg.
		Interpolates, with a spline function,
		using the log of this 2D distribution.
		"""
		#Produces arrays from histograms for sin(dec) vs logE distribution
		#The first array is experimental (to simulate bkg)
		#The second array is Monte Carlo (to simulate signal)
		HBG = self.CreateEnergyHistForSplines(exp['sinDec'],
			exp['logE'], weights=np.ones_like(exp['logE']))
		HSig = self.CreateEnergyHistForSplines(np.sin(mc['trueDec']),
			mc['logE'], weights=self.GetWeights(mc, gamma))

		#Produces an array containing True if x> 0, False otherwise
		DomainBG = HBG > 0.
		DomainSig = HSig > 0.

		#Creates an array of ones as the default ratio
		ratio = np.ones_like(HBG, dtype=np.float)
		#Bitwise Addition giving an Truth Array
		#Returns True if True in both Sig and Bkg, otherwise False
		mask = DomainBG & DomainSig
		#Calculates the ratio sig/bkg for those entries True in Mask array
		ratio[mask] = (HSig[mask] / HBG[mask])

		#Finds the minimum ratio
		min_ratio = np.amax(ratio)
		#Where true in sig and false in bkg, sets ratio to minimum ratio
		np.copyto(ratio, min_ratio, where=DomainSig & ~DomainBG)

		#Sets Bin centers, and order of spline (for x and y)
		sinDecBinCenter = (self.sinDecBins[:-1] + self.sinDecBins[1:]) / 2.
		logEBinCenter = (self.EnergyBins[:-1] + self.EnergyBins[1:]) / 2.
		logE_order = 2

		#Fits a 2D spline function to the log of ratio array
		#This is 2nd order in both dimensions
		spline = scipy.interpolate.RectBivariateSpline(logEBinCenter,
			sinDecBinCenter, np.log(ratio), kx=logE_order, ky=logE_order, s=0)
		self.EnergyPDFSplineInterpolated = spline
		return spline

	def CreateEnergyHistForSplines(self, sinDec, logE, weights, ):
		"""Produces normalised histograms of distribution of variables
		Sin(Dec) and Log(Energy), for use in Spline interpolation.
		"""
		EnergyBins = self.EnergyBins
		sinDec_bins = self.sinDecBins
		#Produces the histogrmas
		h, binedges = np.histogramdd((logE, sinDec),
			bins=(EnergyBins, sinDec_bins),
			weights=weights)
		ndim = h.ndim
		#Normalises histogram
		norms = np.sum(h, axis=ndim - 2)
		norms[norms == 0.] = 1.
		h /= norms
		return h

#==============================================================================
#Background Functions
#==============================================================================

	def create_space_BG_pdf(self, exp, ):
		"""Creates the spatial PDF for background.
		Generates a histogram for the exp. distribution in sin declination.
		Fits a spline function to the distribution, giving a spatial PDF.
		Returns this spatial PDF.
		"""
		sinDec_bins = self.sinDecBins
		sinDec_range = (np.min(sinDec_bins), np.max(sinDec_bins))
		hist, bins = np.histogram(exp['sinDec'], density=True,
			bins=sinDec_bins, range=sinDec_range)
		bins = np.concatenate([bins[:1], bins, bins[-1:]])
		hist = np.concatenate([hist[:1], hist, hist[-1:]])
		bckg_spline = scipy.interpolate.InterpolatedUnivariateSpline(
								(bins[1:] + bins[:-1]) / 2.,
								np.log(hist), k=2)
		return bckg_spline

	def EvaluateB(self, ):
		data = self._ev
		bckg_spline = self.bckg_spline_space
		B_space = 1. / 2. / np.pi * np.exp(bckg_spline(data["sinDec"]))
		if self.UseTime is True:
			B_time = np.ones_like(B_space) / (self.DataEnd - self.DataStart)
		else:
			B_time = np.ones_like(B_space)
		self._ev_B = B_space * B_time

#==============================================================================
#Signal Functions
#==============================================================================

	def SpacePDFSignal(self, source, data, ):
		distance = astro.angular_distance(data['ra'],
			data['dec'], source['ra'], source['dec'])
		SpaceTerm = 1. / ((2. * np.pi * data['sigma'] ** 2.) *
			np.exp(-0.5 * (distance / data['sigma']) ** 2.))
		return SpaceTerm

	def TimePDFSignal(self, source, data, ):
		TimeTerm = self.SingleTimePDF(data['timeMJD'], source)
		return TimeTerm

	def S_source(self, source, data, ):
		SpaceTerm = self.SpacePDFSignal(source, data, )
		if self.UseTime is True:
			TimeTerm = self.TimePDFSignal(source, data, )
		else:
			TimeTerm = np.ones_like(SpaceTerm)
		return SpaceTerm * TimeTerm

	def EvaluateS(self):
		data = self._ev
		sources = self.sources
		NSources = len(sources)
		NData = len(data)
		self.SoB = np.zeros([NSources, NData])
		for i in range(len(sources)):
			self.SoB[i] = self.S_source(sources[i], data, )

#==============================================================================
#Time PDFs
#==============================================================================

	def AnalyticTimePDF(self, t, t_pp, DecayModelLenght):
		t_max = DecayModelLenght
		t = np.asarray(t)
		r = np.zeros_like(t)
		Norm = t_pp * (np.log(t_pp + t_max) - np.log(t_pp))

		mask = np.logical_and(t >= 0., t < t_max)
		r[mask] = (1. + t[mask] / t_pp) ** -1.
		r[mask] = r[mask] / Norm
		return r

	def StepFunc(self, x, x_0):
		return 0.5 * (np.sign(x - x_0) + 1.)

	def BoxFunc(self, x, deltaX):
		Norm = deltaX
		return (self.StepFunc(x, 0.) - self.StepFunc(x, deltaX)) / Norm

	def ComputeSourceWeightsTime(self, ):
		for source in self.sources:
			if self.TimeModel == 'Box':
				t_start = min(max(self.DataStart, source['discoverydate_mjd']),
					self.DataEnd)
				t_end = min(max(self.DataStart, source['discoverydate_mjd']
					+ self.TimeBoxLenght), self.DataEnd)
				TimeLengthInSeasons = t_end - t_start
				TotNorm = self.TimeBoxLenght
				source['weight_time'] = max(TimeLengthInSeasons /
					self.SeasonTimeSpan, 0.)
				source['TimeNorm'] = max(TimeLengthInSeasons / TotNorm, 0.)
			if self.TimeModel == 'BoxPre':
				t_start = min(max(self.DataStart, source['discoverydate_mjd']
					- self.TimeBoxLenght), self.DataEnd)
				t_end = min(max(self.DataStart, source['discoverydate_mjd']),
					self.DataEnd)
				TimeLengthInSeasons = t_end - t_start
				TotNorm = self.TimeBoxLenght
				source['weight_time'] = max(TimeLengthInSeasons /
					self.SeasonTimeSpan, 0.)
				source['TimeNorm'] = max(TimeLengthInSeasons /
					TotNorm, 0.)
			if self.TimeModel == 'Decay':
				t_pp = self.Model_tpp
				t_start = max(0., self.DataStart
					- (source['discoverydate_mjd']))
				t_end = min(self.DecayModelLenght,
					max((self.DataEnd - source['discoverydate_mjd']), 0.))
				TotNorm = t_pp * (np.log(self.DecayModelLenght + t_pp)
					- np.log(0. + t_pp))
				SeasonNorm = t_pp * (np.log(t_end + t_pp)
					- np.log(t_start + t_pp))
				source['TimeNorm'] = SeasonNorm / TotNorm
				source['weight_time'] = ((source['TimeNorm'] * TotNorm) /
					self.SeasonTimeSpan)

	def SingleTimePDF(self, t, source):
		t = np.asarray(t)
		r = np.zeros_like(t)
		if source['TimeNorm'] > 0.:
			Norm = source['TimeNorm']
		else:
			Norm = 1.
		mask = np.logical_and(t > self.DataStart, t < self.DataEnd)
		t = t[mask] - source['discoverydate_mjd']
		r[mask] = self.NuLightCurveFunc(t) / Norm
		return r

	def great_circle_distance(self, ra_1, dec_1, ra_2, dec_2, ):
		delta_dec = np.abs(dec_1 - dec_2)
		delta_ra = np.abs(ra_1 - ra_2)
		x = (np.cos(dec_1) * np.cos(dec_2) * (np.sin(delta_ra / 2.)) ** 2. +
			(np.sin(delta_dec / 2.)) ** 2.)
		return 2. * np.arcsin(np.sqrt(x))
