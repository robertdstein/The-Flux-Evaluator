import numpy as np
import sys
#from scripts.LLh import LLh
from scripts.GenerationControl import GenerationControl
#import os

src_dir = "/afs/ifh.de/user/a/astasik/scratch/PS_Data/"

seed = np.nan
if len(sys.argv) > 1:
	seed = int(sys.argv[1])
	seed = np.random.randint(0, 4294967295)
	np.random.seed(seed)

settings = {'UseEnergy': True,
				'FitGamma': True,
				'FixedGamma': 2.,
				'UseTime': False,
				'UseBoxFunction': True,
				'TimeModel': 'BoxPre',
				'TimeBoxLenght': 368.00423609999416,
				'UseDecayModel': False,
				'Model_tpp': 1.e-2,
				'DecayModelLenght': 30. * 365.,
				'FitWeights': True,
				'UseBox': False,
				'SmearInjection': 0.0,
				'MissTiming': 1. / 1.,
				'SourcePath': src_dir + '/Catalog/catalog_stack10.npy',
			}

GenerationControlInstance = GenerationControl(settings)

path = '/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/results/results'

k_values = np.array([0., 1., 2.0, 3.0, 4.0, 5.0])

np.random.shuffle(k_values)

n_trials = 10

for i in k_values:
	print i
	seed = seed + 1
	GenerationControlInstance.GenerateTestStatistics(k=i, n_trials=n_trials,
		path=path, seed=seed)
	print ''
	print '---------------------'
