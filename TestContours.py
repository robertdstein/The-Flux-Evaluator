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

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

user_dir = "/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/"
file_name = "test_contours.ini"
pickle_names = "test_contours/results_"

test_configs_file = user_dir + file_name

Config = ConfigParser.ConfigParser()

lower_length = 10
upper_length = 300

lengths = np.linspace(lower_length, upper_length, 30)
sim_length = 100

offsets = np.linspace(0, 90, 10)

with open(test_configs_file, "w") as f:
    for length in lengths:
        for offset in offsets:
            if length > offset:
                name = pickle_names + str(length) + "/" + str(offset)
                Config.add_section(name)
                Config.set(name, "UseEnergy", False)
                Config.set(name, "FitGamma", True)
                Config.set(name, "FixedGamma", 2)
                Config.set(name, "UseTime", True)
                Config.set(name, "SimTimeModel", "Box")
                Config.set(name, "SimTimeParameters", {
                    "t0": -100, "length": sim_length})
                Config.set(name, "ReconTimeModel", "Box")
                Config.set(name, "ReconTimeParameters", {"t0": -100 + offset, "length":
                    length})
                Config.set(name, "FitWeights", True)
                Config.set(name, "UseBox", True)
                Config.set(name, "CatName",
                           "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/" +
                           "Catalogue/catalogue00.npy")
                Config.set(name, "DataConfig", "fast.ini")
                maxk = 30 * ((100. / length) ** 0.01 + 1 * np.abs(100. / float(length) - 1)) * \
                       (100. / (100. - float(np.abs(offset))))
                Config.set(
                    name, "MaxK", maxk)

    Config.write(f)

if cfg.submit:
    for section in Config.sections():
         os.system(
             "python " +
             "/afs/ifh.de/user/s/steinrob/Desktop/python/stacking/RunLocal.py" +
             " -c " + section + " -f " + file_name + " -n 100 -s 10")

    #     RC.submit_to_cluster(10, section, file_name, ntrials=200, steps=20)
    # RC.wait_for_cluster()

all_fitted_lengths = []
all_fitted_sens = []

sens_map = np.ones((len(lengths), len(offsets))) * np.nan

print sens_map

for j, offset in enumerate(offsets):

    allfits = []

    datapoints = {
        "length": [],
        "interpolation": [],
        "polynom": [],
        "analytic": []
    }

    print length, offset

    for i, length in enumerate(lengths):
        if length > offset:
            name = pickle_names + str(length) + "/" + str(offset)
            fits = MF.run(name)
            if fits is not None:
                datapoints["length"].append(length)
                datapoints["interpolation"].append(fits["interpolation"])
                datapoints["polynom"].append(fits["polynom"])
                datapoints["analytic"].append(fits["mine"])
                sens_map[i][j] = fits["polynom"]

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

    plt.xlabel("Reconstruction Length (Days)")
    plt.ylabel("Sensitivity (Arbitrary Units)")
    plt.legend()
    plt.tight_layout()
    plot_path = user_dir + "plots/test_contours/offset_" + str(offset) + ".pdf"
    output_dir = os.path.dirname(plot_path)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    print "Saving to", plot_path
    plt.savefig(plot_path)
    plt.close()

plt.figure()
plt.plot(offsets, all_fitted_sens)
plt.xlabel("Offset")
plt.ylabel("Best Sensitivity")
plt.savefig(user_dir + "plots/combined_contour_sens.pdf")
plt.close()

plt.figure()
plt.plot(offsets, all_fitted_lengths)
plt.xlabel("Offset")
plt.ylabel("Reconstruction Length at Best Sensitivity")
plt.savefig(user_dir + "plots/combined_contour_lengths.pdf")
plt.close()

cmap = cm.jet_r

plt.figure()
plt.imshow(
    sens_map, aspect="auto", interpolation="bilinear", cmap=cmap,
    extent=(offsets[0], offsets[-1], lengths[-1], lengths[0])
    )

plt.colorbar()
plt.savefig(user_dir + "plots/combined_contour_map.pdf")
plt.close()


