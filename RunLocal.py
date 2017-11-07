"""
Script to run the stacking code locally. Can be called by shell scripts for
running on cluster.

Through use of argparse, a given configuration for the code can be selected.
This can be given from the command line, in the form:

python RunLocal.py -c Desired_Configuration_Name

Each available configuration must be listed in "config.ini", and controls
options for fitting, such as which catalogue is to be used, and which seasons
of data should be included.

The season dataset are defined in separate data_config files under the
subdirectory data_configs/. In each data_config.ini file, there is a list of
each season to be included, as well as all necessary paths to find the data.

Any arbitrary season combination can be used through creation of a new .ini
data_config file, in which each required season is included. The associated
configuration, which is selected for use this script, must simply have the
name of the data_config file included in the DataConfig field.

"""
import numpy as np
import os
import argparse
import ConfigParser
from scripts.GenerationControl import GenerationControl

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", default="Fast_with_fit")
parser.add_argument("-f", "--conf_file", default="config.ini")
parser.add_argument("-n", "--ntrials", default=10, help="Number of trials")
parser.add_argument("-s", "--step", default=15, help="Number of flux steps")

cfg = parser.parse_args()

# Sets the root directory containing scripts and config files/subdirectories
root = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"

conf = ConfigParser.ConfigParser()
config_path = root + cfg.conf_file

if os.path.isfile(config_path):
    conf.read(config_path)
else:
    raise Exception(
        "Unable to file list of configs. No file found at " + config_path)

# If selected configuration is not found in config.ini, raise error
if cfg.config not in conf.sections():
    print "Searching for config section", cfg.config, "in", conf.sections()
    raise Exception("Config file not found.")

else:
    seed = np.random.randint(0, 4294967295)
    np.random.seed(seed)

    # Finds, from the given configuration, which data_config file is to be used
    data_conf_dir = root + "data_configs/" + conf.get(cfg.config, "DataConfig")
    print data_conf_dir

    settings = {'UseEnergy': eval(conf.get(cfg.config, "UseEnergy")),
                'FitGamma': eval(conf.get(cfg.config, "FitGamma")),
                'FixedGamma': float(conf.get(cfg.config, "FixedGamma")),
                'UseTime': eval(conf.get(cfg.config, "UseTime")),
                'SimTimeModel': conf.get(cfg.config, "SimTimeModel"),
                'SimTimeParameters': eval(conf.get(cfg.config,
                                                   "SimTimeParameters")),
                'ReconTimeModel': conf.get(cfg.config, "ReconTimeModel"),
                'ReconTimeParameters': eval(conf.get(cfg.config,
                                                     "ReconTimeParameters")),
                'FitWeights': eval(conf.get(cfg.config, "FitWeights")),
                'UseBox': eval(conf.get(cfg.config, "UseBox")),
                'SourcePath': conf.get(cfg.config, "CatName"),
                'DataConfig': data_conf_dir,
                'ConfigName': cfg.config
                }

    GenerationControlInstance = GenerationControl(settings)

    # Sets path to save output pickle files
    path = "/afs/ifh.de/user/s/steinrob/scratch/stacking_dump/results"

    # Selects a range of flux scales, for later use in Sensitivity graphs
    # The maximum value should be given in the ConfigFile as MaxK
    k_values = np.linspace(0.0, float(conf.get(cfg.config, "MaxK")), cfg.step)

    # Ensures that, if run over a cluster that can crash, a uniform
    # distribution over k values is sampled
    np.random.shuffle(k_values)

    # Sets number of times each k flux scale is looped over
    n_trials = int(cfg.ntrials)

    # Loops over k, and produces individual pickle result files for each
    for i in k_values:
        print i
        seed += 1
        GenerationControlInstance.generate_test_statistics(
            k=i, n_trials=n_trials, path=path, seed=seed)
        print ''
        print '---------------------'
