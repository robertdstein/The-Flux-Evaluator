import numpy as np
import argparse

#from scripts.LLh import LLh
from scripts.GenerationControl import GenerationControl
#import os
import ConfigParser

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", default="Fast_with_fit")

cfg = parser.parse_args()

root = "/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/"

conf = ConfigParser.ConfigParser()
conf.read(root + "config.ini")

if cfg.config not in conf._sections:
	print "Searching for config section", cfg.config, "in", conf.sections()
	raise Exception("Config file not found.")

else:
	cat_dir = "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/"

	seed = np.random.randint(0, 4294967295)
	np.random.seed(seed)

	settings = {'UseEnergy': eval(conf.get(cfg.config, "UseEnergy")),
					'FitGamma': eval(conf.get(cfg.config, "FitGamma")),
					'FixedGamma': float(conf.get(cfg.config, "FixedGamma")),
					'UseTime': eval(conf.get(cfg.config, "UseTime")),
					'UseBoxFunction': eval(conf.get(
						cfg.config, "BoxFunction")),
					'TimeModel': conf.get(cfg.config, "TimeModel"),
					'TimeBoxLenght': float(conf.get(
						cfg.config, "TimeBoxLength")),
					'UseDecayModel': eval(conf.get(
						cfg.config, "UseDecayModel")),
					'Model_tpp': float(conf.get(cfg.config, "Model_tpp")),
					'DecayModelLenght': eval(conf.get(
						cfg.config, "DecayModelLength")),
					'FitWeights': eval(conf.get(cfg.config, "FitWeights")),
					'UseBox': eval(conf.get(cfg.config, "UseBox")),
					'SmearInjection': float(conf.get(
						cfg.config, "SmearInjection")),
					'MissTiming': eval(conf.get(cfg.config, "MissTiming")),
					'SourcePath': cat_dir + conf.get(cfg.config, "CatName"),
					'DataConfig': conf.get(cfg.config, "DataConfig"),
					'ConfigName': cfg.config
					}

	GenerationControlInstance = GenerationControl(settings)

	path = root + "results/results"

	k_values = np.array([0., 1.0, 2.0])

	np.random.shuffle(k_values)

	n_trials = 10

	for i in k_values:
		print i
		seed = seed + 1
		GenerationControlInstance.GenerateTestStatistics(k=i, n_trials=n_trials,
			path=path, seed=seed)
		print ''
		print '---------------------'