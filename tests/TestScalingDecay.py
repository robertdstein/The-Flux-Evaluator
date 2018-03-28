"""
Code to test impact of sensitivity for misalignment
"""
import os
import numpy as np
import math
import ConfigParser
import argparse
import MergeFiles as MF
import RunCluster as RC
from scripts import time_models as tm
from matplotlib import pyplot as plt
import sys
sys.path.append('..')
from scripts.utils import coenders_7year_sensitivity, flux_to_k, k_to_flux
from common import plot_path, source_path, cat_path, log_path

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

file_name = "analysis_config/test_scaling_decay.ini"
pickle_names = "test_scaling_decay/dec_"

test_configs_file = source_path + file_name

Config = ConfigParser.ConfigParser()

sindecs = [-0.5, 0.00, 0.5]

lengths = np.linspace(10, 300, 30)
sim_length = 100

with open(test_configs_file, "w") as f:
    for x in sindecs:
        for i in lengths:
            name = pickle_names + "{0:.2f}".format(x) + "/length_" + str(i)
            Config.add_section(name)
            sim_params = {"t0": -100, "length": 100, "t_pp": 1.0}
            Config.set(name, "UseEnergy", False)
            Config.set(name, "FitGamma", True)
            Config.set(name, "FixedGamma", 2)
            Config.set(name, "UseTime", True)
            Config.set(name, "SimTimeModel", "Decay")
            Config.set(name, "SimTimeParameters", sim_params)
            Config.set(name, "ReconTimeModel", "Box")
            Config.set(name, "ReconTimeParameters", {"t0": -100, "length": i})
            Config.set(name, "FitWeights", True)
            Config.set(name, "UseBox", True)
            Config.set(name, "CatName", cat_path + "single_source_dec_" +
                       '{0:.2f}'.format(x) + ".npy")
            Config.set(name, "DataConfig", "fast.ini")

            sim = tm.return_integrated_light_curve(sim_params["t0"],
                                                   "Decay", sim_params)

            alt = tm.return_integrated_light_curve(
                i + sim_params["t0"], "Decay", sim_params)

            maxk = 15 * ((1. / ((1 - sim) * alt)) + math.sqrt(i / 100))
            coenders_sens = coenders_7year_sensitivity(x)
            coenders_k_sens = flux_to_k(coenders_sens)

            Config.set(name, "MaxK", coenders_k_sens * maxk)

    Config.write(f)

if cfg.submit:
    for section in Config.sections():
        os.system(
            "python " + source_path + "RunLocal.py" +
            " -c " + section + " -f " + file_name + " -n 50 -s 5")

    #     RC.submit_to_cluster(10, section, file_name, ntrials=100, steps=20)
    # RC.wait_for_cluster()

plt.figure()
ax1 = plt.subplot(111)

for x in sindecs:
    allfits = []

    label = "SinDec = " + str(x)

    datapoints = {
        "length": [],
        "interpolation": [],
        "polynom": [],
        "analytic": []
    }

    for i in lengths:
        name = pickle_names + "{0:.2f}".format(x) + "/length_" + str(i)
        fits = MF.run(name)
        if fits is not None:
            datapoints["length"].append(i)
            datapoints["interpolation"].append(fits["interpolation"])
            datapoints["polynom"].append(fits["polynom"])

    datapoints["polynom_sens"] = k_to_flux(
        np.array(datapoints["polynom"]))

    ax1.plot(datapoints["length"], datapoints["polynom_sens"],
             label=label)

ax1.set_ylim(ymin=1.e-12, ymax=1.e-9)
ax1.legend(loc='upper right', fancybox=True, framealpha=1.)
ax1.grid(True, which='both')
ax1.semilogy(nonposy='clip')
ax1.set_ylabel(r"$E^2 \mathrm{d}N /\mathrm{d}E$ [ TeV cm$^{-2}$ s$^{-1}$ ]",
               fontsize=12)
ax1.set_xlabel("Length (Days)", fontsize=12)
plt.legend()
plt.tight_layout()
save_path = plot_path + "combined_test_scaling_decay.pdf"
plt.savefig(save_path)
