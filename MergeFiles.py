# -*- coding: utf-8 -*-
"""
Script to merge results of stacking code.

Through use of argparse, a given configuration for the code can be selected.
This can be given from the command line, in the form:

python MergeFiles.py -c Desired_Configuration_Name

Each available configuration must be listed in "config.ini", and controls
options for fitting, such as which catalogue is to be used, and which seasons
of data should be included.

Instead of being run from the command line, run() can also be called from
RunCluster.py.

"""
from scripts.GenerationControl import GenerationControl
from scripts.StatisticalEvaluation import Sensitivity

import ConfigParser
import argparse


def run(config_name,
        root="/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/"):
    """Checks if given configuration exists in the config.ini file, with a
    name matching the one given. If so, merges the test results into one
    pickle file, and produces Sensitivity graphs for the given configuration

    :param config_name: Name of config to be used, must be in config.ini
    :param root: path to directory containing results/, merged/, and plot/
    directories.
    """
    conf = ConfigParser.ConfigParser()
    conf.read("config.ini")

    if config_name not in conf.sections():
        print "Searching for config section", config_name, "in", conf.sections()
        raise Exception("Config file not found.")

    else:
        generation_control_instance = GenerationControl()

        path = root + "results/results"
        out_put_path = root + "merged/results"
        plot_path = root + "plots/"

        generation_control_instance.merge_test_result_pickles(
            data_path=path, output_path=out_put_path, config_name=config_name)

        sens = Sensitivity(path=out_put_path, plot_path=plot_path,
                           plotting=True, config=config_name)

        sens.CreateSensitivyAllInOne()

# If script is run from command line, automatically uses run()
# If imported into another script, run() must be explicitly called
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default="Fast_with_fit")

    cfg = parser.parse_args()
    run(cfg.config)
