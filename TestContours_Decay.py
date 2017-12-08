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
import matplotlib.cm as cm
from scripts import time_models as tm

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

user_dir = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"
save_dir = "/afs/ifh.de/user/s/steinrob/scratch/The-Flux-Evaluator__Output/"
file_name = "test_contours_decay.ini"
pickle_names = "test_contours_decay/results_"

test_configs_file = user_dir + file_name

Config = ConfigParser.ConfigParser()

lower_length = 5
upper_length = 150

lengths = np.linspace(lower_length, upper_length, 30)
# lengths = [5.]
sim_length = 50

offsets = np.linspace(-95, 55, 31)

# lengths = np.linspace(10, 90, 5)
# offsets = np.linspace(-75, 50, 6)

source_time = 50

with open(test_configs_file, "w") as f:
    for length in lengths:
        for offset in offsets:
            if ((offset + length) > 0) and (sim_length > offset):
                name = pickle_names + str(length) + "/" + str(offset)
                sim_params = {"t0": source_time, "length": sim_length,
                              "t_pp": 1.0}
                Config.add_section(name)
                Config.set(name, "UseEnergy", False)
                Config.set(name, "FitGamma", True)
                Config.set(name, "FixedGamma", 2)
                Config.set(name, "UseTime", True)
                Config.set(name, "SimTimeModel", "Decay")
                Config.set(name, "SimTimeParameters", sim_params)
                Config.set(name, "ReconTimeModel", "Box")
                Config.set(name, "ReconTimeParameters", {
                    "t0": source_time + offset, "length": length})
                Config.set(name, "FitWeights", True)
                Config.set(name, "UseBox", True)
                Config.set(name, "CatName",
                           "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/" +
                           "Catalogue/catalogue00.npy")
                Config.set(name, "DataConfig", "fast.ini")

                sim = tm.return_integrated_light_curve(source_time + offset,
                                                       "Decay", sim_params)

                alt = tm.return_integrated_light_curve(
                    source_time + offset + sim_length, "Decay", sim_params)

                maxk = 80 * (((1. / (alt - sim)) +
                             math.sqrt(length / 100.)))

                N = alt - sim

                bkg = length

                Nsum = N + math.sqrt(length)

                scale = 20

                print length, offset, maxk, sim, alt

                Config.set(name, "MaxK", maxk)
        # raw_input("prompt")
    Config.write(f)

# raw_input("prompt")

if cfg.submit:

    os.system("rm " + user_dir + "logs/*")

    for section in Config.sections():
         # os.system(
         #     "python " +
         #     "/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/RunLocal.py" +
             # " -c " + section + " -f " + file_name + " -n 100 -s 10")

        RC.submit_to_cluster(10, section, file_name, ntrials=100, steps=20)
    RC.wait_for_cluster()

all_fitted_lengths = []
all_fitted_sens = []

sens_map = np.ones((len(lengths), len(offsets))) * np.nan
log_sens_map = np.ones((len(lengths), len(offsets))) * np.nan

for j, offset in enumerate(offsets):

    allfits = []

    datapoints = {
        "length": [],
        "interpolation": [],
        "polynom": [],
        "analytic": []
    }

    for i, length in enumerate(lengths):
        if ((offset + length) > 0) and (sim_length > offset):
            name = pickle_names + str(length) + "/" + str(offset)
            try:
                fits = MF.run(name)
                if fits is not None:
                    datapoints["length"].append(length)
                    datapoints["interpolation"].append(fits["interpolation"])
                    datapoints["polynom"].append(fits["polynom"])
                    datapoints["analytic"].append(fits["mine"])
                    sens_map[i][j] = fits["polynom"]
                    log_sens_map[i][j] = np.log(fits["polynom"])
            except ValueError:
                pass

    try:
        plt.figure()
        plt.plot(datapoints["length"], datapoints["interpolation"],
                 label="interpolation")
        plt.plot(datapoints["length"], datapoints["polynom"],
                 label="ploynomial")

        n_points = 1

        min_val = min(datapoints["polynom"])
        min_index = datapoints["polynom"].index(min_val)
        low_index = max(min_index-n_points, 0)
        high_index = min(min_index + n_points + 1, len(datapoints["polynom"]))

        dataset = datapoints["polynom"][low_index: high_index]
        min_lengths = datapoints["length"][low_index: high_index]

        parabola = np.polyfit(min_lengths, dataset, 2)

        def f(x):
            return (parabola[0] * (x ** 2)) + (parabola[1] * x) + parabola[2]

        fitted_length = -parabola[1] / (2 * parabola[0])
        fitted_sens = f(fitted_length)

        if (min(datapoints["length"]) < fitted_length) and (
            max(datapoints["length"]) > fitted_length) and (parabola[0] > 0):

            all_fitted_lengths.append(fitted_length)
            all_fitted_sens.append(fitted_sens)

            message = "Minimum occurs at " + str.format('{0:.1f}', fitted_length) + " \n"
            message += "Minimum sensitivity is " + str.format('{0:.2f}', fitted_sens)

            print message

            plt.annotate(message, xy=(0.3, 0.9), xycoords="axes fraction")

            mask = f(np.array(datapoints["length"])) < max(datapoints["polynom"])

            plt.plot(
                np.array(datapoints["length"])[mask],
                f(np.array(datapoints["length"])[mask]),
                linestyle="--", label="Parabola")
            
        else:
            all_fitted_lengths.append(np.nan)
            all_fitted_sens.append(np.nan)

        plt.xlabel("Reconstruction Length (Days)")
        plt.ylabel("Sensitivity (Arbitrary Units)")
        plt.legend()
        plt.tight_layout()
        plot_path = save_dir + "plots/test_contours_decay/offset_" + \
                    str(offset) + ".pdf"
        output_dir = os.path.dirname(plot_path)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        print "Saving to", plot_path
        plt.savefig(plot_path)
        plt.close()

    except ValueError:
        all_fitted_lengths.append(np.nan)
        all_fitted_sens.append(np.nan)

plt.figure()
plt.plot(offsets, all_fitted_sens)
plt.xlabel("Offset")
plt.ylabel("Best Sensitivity")
plt.savefig(save_dir + "plots/combined_contour_decay_sens.pdf")
plt.close()

plt.figure()
plt.plot(offsets, all_fitted_lengths)
plt.xlabel("Offset")
plt.ylabel("Reconstruction Length at Best Sensitivity")
plt.savefig(save_dir + "plots/combined_contour_decay_lengths.pdf")
plt.close()

for i, map_to_plot in enumerate([sens_map, log_sens_map]):

    cmap = cm.jet_r

    plt.figure()
    plt.imshow(
        map_to_plot, aspect="auto", interpolation="nearest", cmap=cmap,
        extent=(offsets[0], offsets[-1], lengths[-1], lengths[0])
        )
    cbar = plt.colorbar()
    cbar.set_label(["Sensitivity", "Log(Sensitivity)"][i])

    mask = np.isnan(sens_map)
    min_val = min(sens_map[~mask])
    loc = np.where(sens_map == min_val)
    plt.scatter(x=offsets[loc[1]], y=lengths[loc[0]], marker="*", color="white")
    plt.xlabel("Offset (Days)")
    plt.ylabel("Reconstruction Length (Days)")

    plt.savefig(save_dir + "plots/combined_contour_decay_map" + ["", "_log"][i]
                + ".pdf")
    plt.close()

for i, length in enumerate(lengths):
    print i, length

    raw_y = sens_map[i]

    mask = ~np.isnan(raw_y)
    x = np.array(offsets)[mask]
    y = raw_y[mask]

    print x
    print y

    plt.figure()
    plt.plot(x, y, label="ploynomial")

    n_points = 1

    min_val = min(y)
    min_index = list(y).index(min_val)
    low_index = max(min_index - n_points, 0)
    high_index = min(min_index + n_points + 1, len(y))

    dataset = y[low_index: high_index]
    min_lengths = y[low_index: high_index]

    parabola = np.polyfit(min_lengths, dataset, 2)

    def f(x):
        return (parabola[0] * (x ** 2)) + (parabola[1] * x) + parabola[2]

    fitted_length = -parabola[1] / (2 * parabola[0])
    fitted_sens = f(fitted_length)

    if (min(x) < fitted_length) and (
                max(x) > fitted_length) and (
        parabola[0] > 0):

        message = "Minimum occurs at " + str.format('{0:.1f}',
                                                    fitted_length) + " \n"
        message += "Minimum sensitivity is " + str.format('{0:.2f}',
                                                          fitted_sens)

        print message

        plt.annotate(message, xy=(0.3, 0.9), xycoords="axes fraction")

        mask = f(x) < max(y)

        plt.plot(
            x[mask],
            f(x)[mask],
            linestyle="--", label="Parabola")

    plt.xlabel("Offset (Days)")
    plt.ylabel("Sensitivity (Arbitrary Units)")
    plt.legend()
    plt.tight_layout()
    plot_path = save_dir + "plots/test_contours_decay/length_" + str(length) + \
                ".pdf"
    output_dir = os.path.dirname(plot_path)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    print "Saving to", plot_path
    plt.savefig(plot_path)
    plt.close()
