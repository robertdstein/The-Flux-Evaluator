# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:13:22 2017

@author: steinrob
"""
from Alex_scripts.GenerationControl import GenerationControl

GenerationControlInstance = GenerationControl()

path = '/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/results/results'
OutPutPath = '/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/scripts/merged/results.pkl'

GenerationControlInstance.MergeTestResultPickles(DataPath=path, OutPutPath = OutPutPath)
#GenerationControlInstance.print_generation_overview(OutPutPath)