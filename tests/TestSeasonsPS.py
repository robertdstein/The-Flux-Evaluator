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
import sys
sys.path.append('..')
from scripts.utils import coenders_7year_sensitivity, flux_to_k, k_to_flux
from common import plot_path, source_path, cat_path, log_path

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

file_name = "analysis_config/test_PS_seasons.ini"
pickle_names = "test_PS_seasons/dec_"

sindecs = np.linspace(0.9, -0.9, 13)

seasons = ["IC40", "IC59", "IC79", "IC86_1", "IC86_2"]
seasons = ["IC40", "IC59"]

test_configs_file = source_path + file_name

Config = ConfigParser.ConfigParser()

with open(test_configs_file, "w") as f:
    for season in seasons:
        for x in sindecs:
            name = pickle_names + season + "_" + "{0:.2f}".format(x) + "/"
            Config.add_section(name)
            Config.set(name, "UseEnergy", True)
            Config.set(name, "FitGamma", True)
            Config.set(name, "FixedGamma", 2)
            Config.set(name, "UseTime", False)
            Config.set(name, "SimTimeModel", "Box")
            Config.set(name, "SimTimeParameters", {"t0": 0, "length": 368})
            Config.set(name, "ReconTimeModel", "Box")
            Config.set(name, "ReconTimeParameters", {"t0": 0, "length": 368})
            Config.set(name, "FitWeights", False)
            Config.set(name, "UseBox", False)
            Config.set(name, "CatName", cat_path + "single_source_dec_" +
                       '{0:.2f}'.format(x) + ".npy")
            Config.set(name, "DataConfig", season)
            coenders_sens = coenders_7year_sensitivity(x)
            coenders_k_sens = flux_to_k(coenders_sens)
            # print name, coenders_sens, coenders_k_sens, 2 * coenders_k_sens

            Config.set(name, "MaxK", 5 * coenders_k_sens)

    Config.write(f)

if cfg.submit:
    for section in Config.sections():
        # os.system(
        #     "python " + source_path + "RunLocal.py" +
        #     " -c " + section + " -f " + file_name + " -n 1 -s 5")

        RC.submit_to_cluster(200, section, file_name, ntrials=1, steps=5)
    RC.wait_for_cluster()

allfits = []

datapoints = {
    "sindec": [],
    "interpolation": [],
    "polynom": [],
    "analytic": []
}

for season in seasons:
    for x in sindecs:
        name = pickle_names + "{0:.2f}".format(x) + "/"
        print name
        try:
            fits = MF.run(name)
            if fits is not None:
                datapoints["sindec"].append(x)
                datapoints["interpolation"].append(fits["interpolation"])
                datapoints["polynom"].append(fits["polynom"])
        except:
            pass

    datapoints["polynom_sens"] = k_to_flux(np.array(datapoints["polynom"]))
    datapoints["polynom_sens"]

    plot_range = np.linspace(-0.99, 0.99, 1000)

    plt.figure()
    ax1 = plt.subplot2grid((4, 1), (0, 0), colspan=3, rowspan=3)
    ax1.plot(plot_range, coenders_7year_sensitivity(plot_range),
             label=r"7 year PS analysis")

    ax1.scatter(
        datapoints["sindec"], datapoints["polynom_sens"], color='black',
        label='This code')

    ax1.plot(
        datapoints["sindec"], datapoints["polynom_sens"], color='black',
        linestyle="--")

    ax1.set_xlim(xmin=-1., xmax=1.)
    ax1.set_ylim(ymin=1.e-13, ymax=1.e-10)
    ax1.legend(loc='upper right', fancybox=True, framealpha=1.)
    ax1.grid(True, which='both')
    ax1.semilogy(nonposy='clip')
    ax1.set_ylabel(r"$E^2 \mathrm{d}N /\mathrm{d}E$ [ TeV cm$^{-2}$ s$^{-1}$ ]",
                   fontsize=12)

    plt.title('7 years PS Sample Sensitivity')

    ax2 = plt.subplot2grid((4, 1), (3, 0), colspan=3, rowspan=1, sharex=ax1)

    ratios = datapoints["polynom_sens"] /\
             coenders_7year_sensitivity(datapoints["sindec"])

    print ratios

    ax2.scatter(datapoints["sindec"], ratios)
    ax2.plot(datapoints["sindec"], ratios, linestyle="--")
    ax2.set_ylabel(r"ratio", fontsize=12)
    ax2.set_xlabel(r"sin($\delta$)", fontsize=12)
    #
    ax1.set_xlim(xmin=-1.0, xmax=1.0)
    ax2.set_ylim(ymin=0.5, ymax=1.5)
    ax2.grid(True)
    xticklabels = ax1.get_xticklabels()
    plt.setp(xticklabels, visible=False)
    plt.subplots_adjust(hspace=0.001)

    # yticks = ax2.yaxis.get_major_ticks()
    # yticks[-1].label1.set_visible(False)


    plt.savefig(plot_path + "SeasonsPS.pdf")
    plt.close()
