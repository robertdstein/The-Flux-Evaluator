"""
Code to test impact of sensitivity for misalignment
"""
import os
import numpy as np
import ConfigParser
import argparse
import MergeFiles as MF
import RunCluster as RC
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import cPickle as pickle
import sys
sys.path.append('..')
from scripts.utils import coenders_7year_sensitivity, flux_to_k, k_to_flux
from common import plot_path, source_path, tde_pickle

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

file_name = "analysis_config/individual_TDEs.ini"

test_configs_file = source_path + file_name

Config = ConfigParser.ConfigParser()
Config.read(file_name)

if cfg.submit:
    for section in Config.sections():
        # print section

        # os.system(
        #     "python " + source_path + "RunLocal.py" +
        #     " -c " + section + " -f " + file_name + " -n 30 -s 5")

        RC.submit_to_cluster(30, section, file_name, ntrials=20, steps=10)
    RC.wait_for_cluster()

results = dict()

for section in Config.sections():
    try:
        fits = MF.run(section)

        print fits

        if fits is not None:

            [name, gamma] = section.split("/")[1].split("_")

            if name not in results.keys():
                results[name] = dict()

            print name, gamma
            datapoints = {
                "interpolation": fits["interpolation"],
                "polynom": fits["polynom"],
                "analytic": fits["mine"]
            }

            results[name][gamma] = datapoints
            print results
    except:
        pass

with open(tde_pickle, "wb") as f:
    pickle.dump(results, f)
