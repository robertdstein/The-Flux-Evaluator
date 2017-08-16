import numpy as np
import copy
import os
import cPickle as pickle
import time
import glob
import ConfigParser

from sys import stdout

import resource

from scipy.optimize import curve_fit
import scipy.optimize
import scipy as scp

from scripts.LLh import LLh


class GenerationControl(object, ):
	def __init__(self, settings=dict()):
		self.settings = settings

	def GenerateTestStatistics(self, k=0, n_trials=10, path='test.pkl',
		seed=np.nan, **kwargs):
		"""Creates the custom path for the given seed and k value
		Generates n_trials and wries results.
		"""
		#***********************************************************************
		#Redundant If statement, potentially
		if True:
			np.random.seed()
			if np.isnan(seed):
				self.seed = np.random.randint(0, 4294967295)
			else:
				self.seed = seed
			np.random.seed(self.seed)
			if self.settings["RunFast"]:
				path += "_Fast"
			else:
				path += "_Full"
			path += '_' + str(k) + '_' + str(self.seed) + '.pkl'
		test_stats = self.GenerateTrials(n_trials, k=k, )
		self.write_result_to_file(path, k, [], test_stats)

	def memory_usage_ps(self, ):
		"""Returns the memory usage.
		"""
		import subprocess
		out = subprocess.Popen(['ps', 'v', '-p', str(os.getpid())],
			stdout=subprocess.PIPE).communicate()[0].split(b'\n')
		vsz_index = out[0].split().index(b'RSS')
		mem = float(out[1].split()[vsz_index]) / 1024 / 1.e3
		return mem

	def GenerateTrials(self, n_trials, k=0., ):
		"""Generates trials.
		Takes as arguments n_trials, and k.
		"""
		test_stats = []

		#***********************************************************************
		#Sets tmp directory (currently Alex's)
		try:
			tmpdir = os.environ['TMPDIR']
		except KeyError:
#			tmpdir ="/afs/ifh.de/user/s/steinrob/scratch/PS_Data"
			tmpdir = "/afs/ifh.de/user/a/astasik/scratch/PS_Data/"

		#Gives the memory usage
		print 0, 'Memory usage: %s (Gb)' % self.memory_usage_ps()
		MemUse = str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1.e6)
		print 'Memory usage max: %s (Gb)' % MemUse

		#Reads in the config variables for each season of data
		data_conf = ConfigParser.ConfigParser()
		if self.settings["RunFast"]:
			data_conf.read("data_config_fast.ini")
		else:
			data_conf.read("data_config.ini")

		#Loops over each season, creating an LLh Object.
		#The data is randomised, and the LLh function is spline-fitted.
		self.seasons = []
		for season in data_conf._sections:
			new_instance = LLh(
				ExpPath=tmpdir + data_conf.get(season, "exp_path"),
				MCPath=tmpdir + data_conf.get(season, "mc_path"),
				AcceptanceWeightPath2=tmpdir + data_conf.get(season, "aw_path"),
				Livetime=eval(data_conf.get(season, "livetime")),
				StartDataTakingMJD=float(data_conf.get(season, "start_mjd")),
				EndDataTankingMJD=float(data_conf.get(season, "end_mjd")),
				**self.settings)
			new_instance.InitEverythingForMultipleTrials()
			self.seasons.append(new_instance)

		#Loops over n_trials
		for counter in range(n_trials):
			#Prints the progress
			self.print_progress(counter, n_trials)

			funcs = []
			for season in self.seasons:
				season.PrepareFakeDataSetAndEvalautePDF(k, )
				f = season.ProduceLLhFunction()
				funcs.append(f)

			def f_final(x):
				"""The likelihood function to be minimised,
				given the Ice Cube Dataset (specified in data_config).
				"""
				#Loops over each season of data
				for season in self.seasons:
					#If both "UseEnergy" and "FitGamma" are true, takes gamma from
					#Otherwise sets Gamma to be the true value
					if np.logical_and(self.settings['UseEnergy'] is True,
						self.settings['FitGamma'] is True):
						gamma = x[-1]
					else:
						gamma = season.InjectionGamma
					for source in season.sources:
						source['weight_acceptance'] = season.AcceptanceFitFunc(
							source['dec'], gamma)

				#****************************************************************
				#Loop over seasons again? Why not joined to previous?
				for season in self.seasons:
					season.sources['weight_distance'] = season.sources['distance'] ** (-2.)
					#What is this for?
					season.sources['weight_acceptance'] = season.sources['weight_acceptance']
					#Either accounts for  time variability, or weights seasons equally
					if self.settings['UseTime']:
						season.sources['weight_time'] = season.sources['weight_time']
					else:
						season.sources['weight_time'] = np.ones_like(season.sources['weight_time'])
					#Assigns an overall source weight
					season.sources['weight'] = (season.sources['weight_distance'] *
												season.sources['weight_time'] *
												season.sources['weight_acceptance'])

				#Puts the weights into a matrix
				WeightMatrix = np.zeros((len(self.seasons), len(self.seasons[0].sources)), dtype=float)
				for i, season in enumerate(self.seasons):
					WeightMatrix[i] = season.sources['weight']

				#***************************************************************************************************************
				#What does 'fitweights' means? Seems the opposite of what I'd expect
				#If the weights are not to be fitted,
				if (self.settings['FitWeights'] is False):
					norm = np.sum(WeightMatrix, axis=1)[:, None]
					for i, n in enumerate(norm):
						if n == 0:
							norm[i] = 1

					SourceWeights = WeightMatrix / norm
					SeasonWeights = np.sum(WeightMatrix, axis=1) / np.sum(np.sum(WeightMatrix, axis=0))

				#If weights are to be fitted,
				if (self.settings['FitWeights'] is True):
					SourceWeights = np.ones_like(WeightMatrix)
					SeasonWeights = WeightMatrix / np.sum(WeightMatrix, axis=0)

				#Loops over seasons to assign weights
				for i, season in enumerate(self.seasons):
					season.SeasonWeight = SeasonWeights[i]
					season.sources['weight'] = SourceWeights[i]

				value = 0.0
				for f in funcs:
					value += f(x)
				return value

			test_stat_res = self.MinimizeLLh(f_final)
			test_stats.append(test_stat_res)
			del(f_final)

		for season in self.seasons:
			del(season.SoB)
		return np.array(test_stats)

	def MinimizeLLh(self, f_final):
		#Gives number of sources FOR GIVEN TRIAL?
		nSources = len(np.load(self.settings['SourcePath']))

		#Checks whether to fit weights and gamma, and whether to use energy
		#Then minimises for the chosen configuration
		#Provides appropriate bounds the chosen configuration

		#If fit weights loops over a bound for each source!
		#If not, then gives only one set of bound

		if self.settings['FitWeights'] is True:
			if self.settings['UseEnergy'] is True:

				if self.settings['FitGamma'] is True:
					seed = np.append(np.ones(nSources), 2.)
					bounds = [(0., 1000.) for i in range(nSources)] + [(1., 4.)]
					res = scp.optimize.fmin_l_bfgs_b(f_final, seed,
						bounds=bounds, approx_grad=True)

				if self.settings['FitGamma'] is False:
					seed = np.ones(nSources)
					bounds = [(0., 1000.) for i in range(nSources)]
					res = scp.optimize.fmin_l_bfgs_b(f_final, seed,
						bounds=bounds, approx_grad=True)

			if self.settings['UseEnergy'] is False:
				seed = np.zeros(nSources)
				bounds = [(0., 1000.) for i in range(nSources)]
				res = scp.optimize.fmin_l_bfgs_b(f_final, seed,
					bounds=bounds, approx_grad=True)

		if self.settings['FitWeights'] is False:
			if self.settings['UseEnergy'] is True:

				if self.settings['FitGamma'] is True:
					rranges = (slice(0., 4., 0.5), slice(2., 4., 0.25))
					seed = scipy.optimize.brute(f_final, rranges, finish=None)
					res = scp.optimize.fmin_l_bfgs_b(f_final, seed,
						bounds=[(0., 1000.), (1., 4.)], approx_grad=True)

				if self.settings['FitGamma'] is False:
					res = scp.optimize.fmin_l_bfgs_b(f_final, 1.,
						bounds=[(0., 1000.)], approx_grad=True)

			if self.settings['UseEnergy'] is False:
				res = scp.optimize.fmin_l_bfgs_b(f_final, 1.,
					bounds=[(0., 1000.)], approx_grad=True)

		#Returns a value for LLH at minimum
		test_stat_res = -res[1]
		return test_stat_res

	def write_result_to_file(self, path, k, n_fit, test_stats):
		"""Adds information to a pickle file.
		"""
		#Checks if Pickle file exists
		self.check_if_pkl_file_exists(path)
		#******************************************************************
		#Opens the pickle file and creates a copy of the data?
		pkl_file = open(path, 'rb')
		test_stat_results = pickle.load(pkl_file)
		pkl_file.close()

		if k in test_stat_results.keys():
			test_stat_results[k] = np.append(test_stat_results[k], test_stats)
		else:
			test_stat_results[k] = test_stats
		#Saves the updated pickle file
		pkl_file = open(path, 'wb')
		pickle.dump(test_stat_results, pkl_file)
		pkl_file.close()

	def check_if_pkl_file_exists(self, path, ):
		"""Checks if the pickle file exists.
		If not, creates an empty pickle file.
		"""
		if os.path.isfile(path):
			pass
		else:
			test_stat_results = {}
			pkl_file = open(path, 'wb')
			pickle.dump(test_stat_results, pkl_file)
			pkl_file.close()
			print path, ' created'

	def print_generation_overview(self, path):
		"""Loads the pickle results file.
		Prints the number of entries for each flux.
		"""
		self.check_if_pkl_file_exists(path)
		pkl_file = open(path, 'rb')
		test_stat_results = pickle.load(pkl_file)
		pkl_file.close()
		print ''
		for key in np.sort(test_stat_results.keys()):
			print key, len(test_stat_results[key])
		print ''
		del(test_stat_results)

	def print_progress(self, counter, n_trials):
		"""Prints the progress (fraction)
		"""
		stdout.write("\r%.1f %%" % (float(counter + 1.) / n_trials * 100.))
		stdout.flush()

	def MergeTestResultPickles(self, DataPath='test_stat_results/test/test',
			OutPutPath="test_stat_results/test", RunFast=False):
		"""Searches for alll pickle files matching the DataPath path.
		Merges these into a single pickle file, and saves it to OutPutPath.
		Prints the combined results.
		"""
		if RunFast:
			DataPath += "_Fast"
			OutPutPath += "_Fast.pkl"
		else:
			DataPath += "_Full"
			OutPutPath += "_Full.pkl"

		#Creates an empty Pickle Path for results
		if os.path.isfile(OutPutPath):
			os.remove(OutPutPath)

		#returns a list of all files matching the path given
		FileList = glob.glob(DataPath + '_*.pkl')

		test_stat_results = {}

		#Loops over each file and adds it to the combined File
		for SingleFile in FileList:
			pkl_file = open(SingleFile, 'rb')
			SingleResult = pickle.load(pkl_file)
			pkl_file.close()
			for k in SingleResult.keys():
				if k in test_stat_results.keys():
					test_stat_results[k] = np.append(test_stat_results[k], SingleResult[k])
				else:
					test_stat_results[k] = SingleResult[k]

		#Saves merged file and prints results
		pkl_file = open(OutPutPath, 'wb')
		pickle.dump(test_stat_results, pkl_file)
		pkl_file.close()
		self.print_generation_overview(OutPutPath)

	##*******************************************************************************************************************************
	##Would this just crash? x is undefined!
	##I think this function is not used
	#def WeighterFunction(self, seasons, ):
		#for season in seasons:
			#if np.logical_and(self.settings['UseEnergy']==True, self.settings['FitGamma']==True):
				#gamma = x[-1]
			#else:
				#gamma = season.InjectionGamma
			#for source in season.sources:
				#source['weight_acceptance'] = season.AcceptanceFitFunc(source['dec'], gamma)

			#for season in seasons:
				#season.sources['weight_distance'] = season.sources['distance']**(-2.)
				#season.sources['weight_acceptance'] = season.sources['weight_acceptance']
				#if self.settings['UseTime']:
					#season.sources['weight_time'] = season.sources['weight_time']
				#else:
					#season.sources['weight_time'] = np.ones_like(season.sources['weight_time'])
						#
				#season.sources['weight'] = (season.sources['weight_distance']*
											#season.sources['weight_time']*
											#season.sources['weight_acceptance'])

				#print np.sum(season.sources['weight'])
			#print '=============================================='

			#WeightMatrix = np.zeros((len(seasons), len(seasons[0].sources)), dtype=float)
			#for i, season in enumerate(seasons):
				#WeightMatrix[i] = season.sources['weight']

			#if (self.settings['FitWeights']==False):
				#norm = np.sum(WeightMatrix, axis=1)[:,None]
				#mask = np.argwhere(norm==0.)
				#norm[mask] = np.ones_like(norm[mask])
				#SourceWeights = WeightMatrix / norm
				#SeasonWeights = np.sum(WeightMatrix, axis=1) / np.sum(np.sum(WeightMatrix, axis=0))
			#if (self.settings['FitWeights']==True):
				#SourceWeights = np.ones_like(WeightMatrix)
				#SeasonWeights = WeightMatrix / np.sum(WeightMatrix, axis=0)

			#for i, season in enumerate(seasons):
				#season.SeasonWeight = SeasonWeights[i]
				#season.sources['weight'] = SourceWeights[i]

			##This doesn't run, as far as I can tell
			#if False:
				#print('')
				#print('')
				#print WeightMatrix
				#print('')
				#for season in seasons:
					#print season.SeasonWeight#, season.sources['weight']