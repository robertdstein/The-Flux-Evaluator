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
from matplotlib import pyplot as plt
import sys
sys.path.append('..')
from scripts.utils import coenders_7year_sensitivity, flux_to_k, k_to_flux
from common import plot_path, source_path, cat_path, log_path


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

file_name = "analysis_config/test_scaling_box.ini"
pickle_names = "test_scaling_box/dec_"

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
            Config.set(name, "UseEnergy", False)
            Config.set(name, "FitGamma", True)
            Config.set(name, "FixedGamma", 2)
            Config.set(name, "UseTime", True)
            Config.set(name, "SimTimeModel", "Box")
            Config.set(name, "SimTimeParameters", {
                "t0": -100, "length": sim_length})
            Config.set(name, "ReconTimeModel", "Box")
            Config.set(name, "ReconTimeParameters", {"t0": -100, "length": i})
            Config.set(name, "FitWeights", True)
            Config.set(name, "UseBox", True)
            Config.set(name, "CatName", cat_path + "single_source_dec_" +
                       '{0:.2f}'.format(x) + ".npy")
            Config.set(name, "DataConfig", "fast.ini")
            maxk = 30 * ((100. / i)**0.01 + 1 * np.abs(100./float(i) - 1))

            coenders_sens = coenders_7year_sensitivity(x)
            coenders_k_sens = flux_to_k(coenders_sens)

            Config.set(name, "MaxK", coenders_k_sens * maxk)

    Config.write(f)


if cfg.submit:
    for section in Config.sections():
        os.system(
            "python " + source_path + "RunLocal.py" +
            " -c " + section + " -f " + file_name + " -n 200 -s 5")

    #     RC.submit_to_cluster(10, section, file_name, ntrials=200, steps=20)
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
save_path = plot_path + "combined_test_scaling_box.pdf"
plt.savefig(save_path)

print "Saving to", save_path

# plt.figure()
# plt.plot(datapoints["length"], datapoints["interpolation"],
#          label="interpolation")
# plt.plot(datapoints["length"], datapoints["polynom"],
#          label="ploynomial")
#
# n_points = 1
#
# min_val = min(datapoints["polynom"])
# min_index = datapoints["polynom"].index(min_val)
# low_index = max(min_index-n_points, 0)
# high_index = min(min_index + n_points + 1, len(datapoints["polynom"]))
#
# dataset = datapoints["polynom"][low_index: high_index]
# min_offsets = datapoints["length"][low_index: high_index]
#
# parabola = np.polyfit(min_offsets, dataset, 2)
#
# def f(x):
#     return (parabola[0] * (x ** 2)) + (parabola[1] * x) + parabola[2]
#
# fitted_offset = -parabola[1] / (2 * parabola[0])
# fitted_sens = f(fitted_offset)
#
# message = "Minimum occurs at " + str.format('{0:.1f}', fitted_offset) + " \n"
# message += "Minimum sensitivity is " + str.format('{0:.2f}', fitted_sens)
#
# print message
#
# plt.annotate(message, xy=(0.3, 0.9), xycoords="axes fraction")
#
# mask = f(np.array(datapoints["length"])) < max(datapoints["polynom"])
#
# plt.plot(
#     np.array(datapoints["length"])[mask],
#     f(np.array(datapoints["length"])[mask]),
#     linestyle="--", label="Parabola")
#
# plt.xlabel("Reconstruction Length (Days)")
# plt.ylabel("Sensitivity (Arbitrary Units)")
# plt.legend()
# plt.tight_layout()
# plot_path = user_dir + "plots/combined_test_scaling.pdf"
# plt.savefig(plot_path)
#
# print "Saving to", plot_path
