import numpy as np
import copy
import os
import cPickle as pickle
import time
import glob
import sys

from sys import stdout

import resource

from scipy.optimize import curve_fit
import scipy.optimize
import scipy as scp

from scripts.LLh import LLh

class GenerationControl(object, ):
	
	def __init__(self, settings=dict()):
		self.settings = settings
		
	def GenerateTestStatistics(self, k=0, n_trials=10, path='test.pkl', seed=np.nan, **kwargs):
		"""Creates the custom path for the given seed and k value
		Generates n_trials and wries results.
		"""
		#**********************************************************************************************************
		#Redundant If statement, potentially
		if True:
			np.random.seed()
			if np.isnan(seed):
				self.seed = np.random.randint(0, 4294967295)
			else:
				self.seed=seed
			np.random.seed(self.seed)
			path = path+'_'+str(k)+'_'+str(self.seed)+'.pkl'
		test_stats = self.GenerateTrials(n_trials, k=k, )
		self.write_result_to_file(path, k, [], test_stats)
		
	def memory_usage_ps(self, ):
		"""Returns the memory usage.
		"""
		import subprocess
		import os
		out = subprocess.Popen(['ps', 'v', '-p', str(os.getpid())],stdout=subprocess.PIPE).communicate()[0].split(b'\n')
		vsz_index = out[0].split().index(b'RSS')
		mem = float(out[1].split()[vsz_index]) / 1024 / 1.e3
		return mem
		
		
	def GenerateTrials(self, n_trials, k=0., ):
		"""Generates trials. 
		Takes as arguments n_trials, and k.
		"""

		test_stats = []
		
		#*************************************************************************************************************************
		#Sets tmp directory (currently Alex's)
		try:
			tmpdir=os.environ['TMPDIR']
		except KeyError:
#			tmpdir ="/afs/ifh.de/user/s/steinrob/scratch/PS_Data"
			tmpdir="/afs/ifh.de/user/a/astasik/scratch/PS_Data/"
			
		#Gives the memory usage
		print 0, 'Memory usage: %s (Gb)' % self.memory_usage_ps()
		MemUse = str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1.e6)
		print 'Memory usage max: %s (Gb)' % MemUse
		
		#Gives the Log Likelihood for each Ice Cube configuration
		#Should try with only one
		
#		self.LLh_instanceIC40 = LLh(
#				   ExpPath=tmpdir+'/FinalSample/IC40/exp/IC40_exp_corrected.npy',
#				   MCPath=tmpdir+'/FinalSample/IC40/mc/IC40_nugen_corrected.npy',
#				   AcceptanceWeightPath2 = tmpdir+'/DeclinationAcceptance/IC40',
#				   Livetime= 375.539,
#				   StartDataTakingMJD=54561.4746759,
#				   EndDataTankingMJD=54964.1892245,
#				   **self.settings)   
#		
#		self.LLh_instanceIC59 = LLh(
#				   ExpPath=tmpdir+'/FinalSample/IC59/exp/IC59_exp_corrected.npy',
#				   MCPath=tmpdir+'/FinalSample/IC59/mc/IC59_nugen_corrected.npy',
#				   AcceptanceWeightPath2 = tmpdir+'/DeclinationAcceptance/IC59',
#				   Livetime= 348.138,
#				   StartDataTakingMJD=54964.1892245,
#				   EndDataTankingMJD=55347.2862153,
#				   **self.settings)
#		
#		self.LLh_instanceIC79 = LLh(
#				   ExpPath=tmpdir+'/FinalSample/IC79/exp/IC79_exp_corrected.npy',
#				   MCPath=tmpdir+'/FinalSample/IC79/mc/IC79_nugen_corrected.npy',
#				   AcceptanceWeightPath2 = tmpdir+'/DeclinationAcceptance/IC79',
#				   Livetime= 315.506,
#				   StartDataTakingMJD=55347.2862153,
#				   EndDataTankingMJD=55694.4164699,
#				   **self.settings)

		self.LLh_instanceIC86_1 = LLh(
				   ExpPath= tmpdir+'/FinalSample/IC86_1/exp/IC86_1_exp_corrected.npy',
				   MCPath= tmpdir+'/FinalSample/IC86_1/mc/IC86_1_nugen_corrected.npy',
				   AcceptanceWeightPath2 = tmpdir+'/DeclinationAcceptance/IC86_1',
				   Livetime=332.61,
				   StartDataTakingMJD=55694.4164699,
				   EndDataTankingMJD=56062.420706,
				   **self.settings)          
			   
#		self.LLh_instanceIC86_2AndFollowing = LLh(
#				   ExpPath=tmpdir+'/FinalSample/IC86_2AndFollowing/exp/IC86_2_3_4_exp_corrected.npy',
#				   MCPath=tmpdir+'/FinalSample/IC86_2/mc/IC86_2_nugen_corrected.npy',
#				   AcceptanceWeightPath2 = tmpdir+'/DeclinationAcceptance/IC86_2AndFollowing',
#				   Livetime=330.38 + 359.95 + 367.21,
#				   StartDataTakingMJD=56062.420706,
#				   EndDataTankingMJD=57160.0440856,
#				   **self.settings)
							
		#Runs "InitEverything" for each configuration
							
#		self.LLh_instanceIC40.InitEverythingForMultipleTrials()
#		self.LLh_instanceIC59.InitEverythingForMultipleTrials()
#		self.LLh_instanceIC79.InitEverythingForMultipleTrials()
		self.LLh_instanceIC86_1.InitEverythingForMultipleTrials()
#		self.LLh_instanceIC86_2AndFollowing.InitEverythingForMultipleTrials()
			
		#Loops over n_trials
		for counter in range(n_trials):
			#Prints the progress			
			self.print_progress(counter, n_trials)

#			self.LLh_instanceIC40.PrepareFakeDataSetAndEvalautePDF(k, )
#			self.f_IC40 = self.LLh_instanceIC40.ProduceLLhFunction()
#
#			self.LLh_instanceIC59.PrepareFakeDataSetAndEvalautePDF(k, )
#			self.f_IC59 = self.LLh_instanceIC59.ProduceLLhFunction()
#
#			self.LLh_instanceIC79.PrepareFakeDataSetAndEvalautePDF(k, )
#			self.f_IC79 = self.LLh_instanceIC79.ProduceLLhFunction()

			#Apparently the best dataset

			self.LLh_instanceIC86_1.PrepareFakeDataSetAndEvalautePDF(k, )
			self.f_IC86_1 = self.LLh_instanceIC86_1.ProduceLLhFunction()
			
#			self.LLh_instanceIC86_2AndFollowing.PrepareFakeDataSetAndEvalautePDF(k, )
#			self.f_IC86_2AndFollowing = self.LLh_instanceIC86_2AndFollowing.ProduceLLhFunction()

			def f_final(x):
				"""The likelihood function to be minimised, given the Ice Cube Data
				"""
				#*****************************************************************************************************************************************************
				#Is this used?
				NSignalTotal = 0.
				
				#Loops over seasons				
				
#				seasons = [self.LLh_instanceIC40, self.LLh_instanceIC59, self.LLh_instanceIC79, self.LLh_instanceIC86_1, self.LLh_instanceIC86_2AndFollowing]
				seasons = [self.LLh_instanceIC86_1]
				
				#Loops over each season of data
				for season in seasons:
					#If both "UseEnergy" and "FitGamma" are true, takes gamma from 
					#Otherwise sets Gamma to be the true value
					if np.logical_and(self.settings['UseEnergy']==True, self.settings['FitGamma']==True):
						gamma = x[-1] 					
					else:
						gamma = season.InjectionGamma
					for source in season.sources:
						source['weight_acceptance'] = season.AcceptanceFitFunc(source['dec'], gamma)
				
				#***************************************************************************************************************************************************
				#Loop over seasons again? Why not joined to previous?
				for season in seasons:
					season.sources['weight_distance'] = season.sources['distance']**(-2.)
					#What is this for?
					season.sources['weight_acceptance'] = season.sources['weight_acceptance']
					#Either accounts for  time variability, or weights seasons equally
					if self.settings['UseTime']:
						season.sources['weight_time'] = season.sources['weight_time']
					else:
						season.sources['weight_time'] = np.ones_like(season.sources['weight_time'])
					#Assigns an overall source weight	
					season.sources['weight'] = (season.sources['weight_distance']*
												season.sources['weight_time']*
												season.sources['weight_acceptance'])
				
				#Puts the weights into a matrix
				WeightMatrix = np.zeros((len(seasons), len(seasons[0].sources)), dtype=float)
				for i, season in enumerate(seasons):
					WeightMatrix[i] = season.sources['weight']
					
				#***************************************************************************************************************
				#What does 'fitweights' means? Seems the opposite of what I'd expect
				#If the weights are not to be fitted, 							   
				if (self.settings['FitWeights']==False):
					norm = np.sum(WeightMatrix, axis=1)[:,None]
					for i, n in enumerate(norm):
						if n==0:
							norm[i]=1

					SourceWeights = WeightMatrix / norm
					SeasonWeights = np.sum(WeightMatrix, axis=1) / np.sum(np.sum(WeightMatrix, axis=0))
				
				#If weights are to be fitted, 				
				if (self.settings['FitWeights']==True):
					SourceWeights = np.ones_like(WeightMatrix)
					SeasonWeights = WeightMatrix / np.sum(WeightMatrix, axis=0)
				
				#Loops over seasons to assign weights
				for i, season in enumerate(seasons):
					season.SeasonWeight = SeasonWeights[i]
					season.sources['weight'] = SourceWeights[i]

				#***********************************************************************************************************************
				#Doesn't run 
				if False:
					print('')
					print('')
					print('')    
					for season in seasons:
						print season.SeasonWeight, season.sources['weight']       
			
#				return self.f_IC40(x)+self.f_IC59(x)+self.f_IC79(x)+self.f_IC86_1(x)+self.f_IC86_2AndFollowing(x)
				return self.f_IC86_1(x)
			
			
			test_stat_res = self.MinimizeLLh(f_final)
			test_stats.append(test_stat_res)
			del(f_final)
		
		del(self.LLh_instanceIC86_1.SoB)
		#*************************************************************************************************************************
		#Is this actually used before being reinitialised?
		MemUse = str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1.e6)
		return np.array(test_stats)
		
			
	def MinimizeLLh(self, f_final):
		#Gives number of sources FOR GIVEN TRIAL?
		nSources = len(np.load(self.settings['SourcePath']))
		
		#Checks whether to fit weights and gamma, and whether to use energy
		#Then minimises for the chosen configuration 
		#Provides appropriate bounds the chosen configuration
		
		#If fit weights loops over a bound for each source!
		#If not, then gives only one set of bound
		
		if self.settings['FitWeights'] == True:
			if self.settings['UseEnergy']==True:
				
				if self.settings['FitGamma']==True:
					seed = np.append(np.ones(nSources), 2.)
					bounds = [(0., 1000.) for i in range(nSources)] + [(1., 4.)]
					res = scp.optimize.fmin_l_bfgs_b(f_final, seed, bounds=bounds, approx_grad=True)
					
				if self.settings['FitGamma']==False:
					seed = np.ones(nSources)                 
					bounds = [(0., 1000.) for i in range(nSources)]
					res = scp.optimize.fmin_l_bfgs_b(f_final, seed, bounds=bounds, approx_grad=True)
					
			if self.settings['UseEnergy']==False:
				seed = np.zeros(nSources)
				bounds = [(0., 1000.) for i in range(nSources)]
				res = scp.optimize.fmin_l_bfgs_b(f_final, seed, bounds=bounds, approx_grad=True)
				
		if self.settings['FitWeights'] == False:
			if self.settings['UseEnergy'] == True:
				
				if self.settings['FitGamma']==True:
					rranges = (slice(0., 4., 0.5), slice(2., 4., 0.25))
					seed = scipy.optimize.brute(f_final, rranges, finish=None)
					res = scp.optimize.fmin_l_bfgs_b(f_final, seed, bounds=[(0., 1000.), (1., 4.)], approx_grad=True)
					
				if self.settings['FitGamma']==False:
					res = scp.optimize.fmin_l_bfgs_b(f_final, 1., bounds=[(0., 1000.)], approx_grad=True)
					
			if self.settings['UseEnergy'] == False:
				res = scp.optimize.fmin_l_bfgs_b(f_final, 1., bounds=[(0., 1000.)], approx_grad=True)           

		#****************************************************************************************************************************8
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
		#*******************************************************************************8
		#What triggers this 'else'?
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
		"""Loads the pickle results file, and prints the results
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
		stdout.write("\r%.1f %%" % (float(counter+1.)/n_trials*100.))
		stdout.flush()

		
	def MergeTestResultPickles(self, DataPath='test_stat_results/test/test', OutPutPath = 'test_stat_results/test.pkl'):
		"""Searches for alll pickle files matching the DataPath path.
		Merges these into a single pickle file, and saves it to OutPutPath.
		Prints the combined results.
		"""
		#Creates an empty Pickle Path for results
		if os.path.isfile(OutPutPath):
			os.remove(OutPutPath)
#		self.check_if_pkl_file_exists(OutPutPath)
#		
		#returns a list of all files matching the path given
		FileList = glob.glob(DataPath+'_*.pkl')
#		
		#**************************************************************************************
#		#Opens the empty file, and basically gives {} ?
#		pkl_file = open(OutPutPath, 'rb')
#		test_stat_results = pickle.load(pkl_file)
#		pkl_file.close()
		
		test_stat_results={}

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
		
	#*******************************************************************************************************************************
	#Would this just crash? x is undefined!
	#I think this function is not used
	def WeighterFunction(self, seasons, ):
		for season in seasons:
			if np.logical_and(self.settings['UseEnergy']==True, self.settings['FitGamma']==True):
				gamma = x[-1] 
			else:
				gamma = season.InjectionGamma
			for source in season.sources:
				source['weight_acceptance'] = season.AcceptanceFitFunc(source['dec'], gamma)
				
			for season in seasons:
				season.sources['weight_distance'] = season.sources['distance']**(-2.)
				season.sources['weight_acceptance'] = season.sources['weight_acceptance']
				if self.settings['UseTime']:
					season.sources['weight_time'] = season.sources['weight_time']
				else:
					season.sources['weight_time'] = np.ones_like(season.sources['weight_time'])
						
				season.sources['weight'] = (season.sources['weight_distance']*
											season.sources['weight_time']*
											season.sources['weight_acceptance'])
				
				print np.sum(season.sources['weight'])
			print '=============================================='
				
			WeightMatrix = np.zeros((len(seasons), len(seasons[0].sources)), dtype=float)
			for i, season in enumerate(seasons):
				WeightMatrix[i] = season.sources['weight']
							   
			if (self.settings['FitWeights']==False):
				norm = np.sum(WeightMatrix, axis=1)[:,None]
				mask = np.argwhere(norm==0.)
				norm[mask] = np.ones_like(norm[mask])
				SourceWeights = WeightMatrix / norm
				SeasonWeights = np.sum(WeightMatrix, axis=1) / np.sum(np.sum(WeightMatrix, axis=0))
			if (self.settings['FitWeights']==True):
				SourceWeights = np.ones_like(WeightMatrix)
				SeasonWeights = WeightMatrix / np.sum(WeightMatrix, axis=0)
					
			for i, season in enumerate(seasons):
				season.SeasonWeight = SeasonWeights[i]
				season.sources['weight'] = SourceWeights[i]

			#This doesn't run, as far as I can tell
			if False:
				print('')
				print('')    
				print WeightMatrix
				print('')    
				for season in seasons:
					print season.SeasonWeight#, season.sources['weight']
