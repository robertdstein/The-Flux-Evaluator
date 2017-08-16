# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:13:22 2017

@author: steinrob
"""
from scripts.GenerationControl import GenerationControl
from scripts.StatisticalEvaluation import sensitivity

GenerationControlInstance = GenerationControl()

root = "/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/"
path = root + "results/results"
OutPutPath = root + "merged/results"
PlotPath = root + "plots/"
RunFast = True

GenerationControlInstance.MergeTestResultPickles(DataPath=path,
	OutPutPath=OutPutPath, RunFast=RunFast)

sens = sensitivity(path=OutPutPath, plot_path=PlotPath,
	plotting=True, RunFast=RunFast)
sens.CreateSensitivyAllInOne()