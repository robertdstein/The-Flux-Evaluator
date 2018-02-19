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
from common import results_path, plot_path, merge_path
import ConfigParser
import argparse
import os


def run(config_name):
    """Checks if given configuration exists in the config.ini file, with a
    name matching the one given. If so, merges the test results into one
    pickle file, and produces Sensitivity graphs for the given configuration

    :param config_name: Name of config to be used, must be in config.ini
    """

    generation_control_instance = GenerationControl()

    path = results_path
    out_put_path = merge_path
    out_put_dir = os.path.dirname(out_put_path)
    if not os.path.isdir(out_put_dir):
        os.makedirs(out_put_dir)
        print "Making", out_put_dir

    generation_control_instance.merge_test_result_pickles(
        data_root=path, output_root=out_put_path, config_name=config_name)

    sens = Sensitivity(path=out_put_path, plot_path=plot_path,
                       plotting=True, config=config_name)

    fits = sens.CreateSensitivyAllInOne()
    return fits

# If script is run from command line, automatically uses run()
# If imported into another script, run() must be explicitly called
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default="Fast_with_fit")
    parser.add_argument("-f", "--conf_file", default="config.ini")

    # Sets the source_path directory containing scripts and config files/subdirectories
    root = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"

    cfg = parser.parse_args()

    config_path = root + "analysis_config/" + cfg.conf_file

    conf = ConfigParser.ConfigParser()
    conf.read(config_path)

    if cfg.config not in conf.sections():
        print "Searching for config section", cfg.config, "in", conf.sections()
        raise Exception("Config file not found.")

    else:
        run(cfg.config)
