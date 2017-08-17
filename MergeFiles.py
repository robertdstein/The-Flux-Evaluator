# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:13:22 2017

@author: steinrob
"""
from scripts.GenerationControl import GenerationControl
from scripts.StatisticalEvaluation import sensitivity

import ConfigParser
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", default="Fast_with_fit")

cfg = parser.parse_args()

conf = ConfigParser.ConfigParser()
conf.read("config.ini")

if cfg.config not in conf._sections:
	print "Searching for config section", cfg.config, "in", conf.sections()
	raise Exception("Config file not found.")

else:
	GenerationControlInstance = GenerationControl()

	root = "/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/"
	path = root + "results/results"
	OutPutPath = root + "merged/results"
	PlotPath = root + "plots/"

	GenerationControlInstance.MergeTestResultPickles(DataPath=path,
		OutPutPath=OutPutPath, ConfigName=cfg.config)

	sens = sensitivity(path=OutPutPath, plot_path=PlotPath,
		plotting=True, ConfigName=cfg.config)
	sens.CreateSensitivyAllInOne()