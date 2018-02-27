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
from common import root_path, source_path, cat_path, log_path, output_path

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

file_name = "analysis_config/test_daifang_misalignment.ini"
pickle_names = "test_misalignment_daifang/offset_"

test_configs_file = source_path + file_name

Config = ConfigParser.ConfigParser()

offsets = np.linspace(0, 90, 1)

with open(test_configs_file, "w") as f:
    for i in offsets:
        name = pickle_names + str(i)
        Config.add_section(name)
        Config.set(name, "UseEnergy", False)
        Config.set(name, "FitGamma", True)
        Config.set(name, "FixedGamma", 2)
        Config.set(name, "UseTime", False)
        Config.set(name, "SimTimeModel", "Box")
        Config.set(name, "SimTimeParameters", {"t0": 0, "length": 368})
        Config.set(name, "ReconTimeModel", "Box")
        Config.set(name, "ReconTimeParameters", {"t0": i -0, "length": 368})
        Config.set(name, "FitWeights", False)
        Config.set(name, "UseBox", False)
        Config.set(name, "CatName", cat_path + "catalogue00.npy")
        Config.set(name, "DataConfig", "fast.ini")
        Config.set(name, "MaxK", 100. * (100. / (100. - float(np.abs(i)))))

    Config.write(f)

if cfg.submit:
    for section in Config.sections():
        os.system(
            "python " + source_path + "RunLocal.py" +
            " -c " + section + " -f " + file_name + " -n 100 -s 5")

    #     RC.submit_to_cluster(10, section, file_name, ntrials=200, steps=20)
    # RC.wait_for_cluster()

allfits = []

datapoints = {
    "offset": [],
    "interpolation": [],
    "polynom": [],
    "analytic": []
}

for i in offsets:
    name = pickle_names + str(i)
    fits = MF.run(name)
    if fits is not None:
        datapoints["offset"].append(i)
        datapoints["interpolation"].append(fits["interpolation"])
        datapoints["polynom"].append(fits["polynom"])

plt.figure()
plt.plot(datapoints["offset"], datapoints["interpolation"],
         label="Interpolation")
plt.plot(datapoints["offset"], datapoints["polynom"],
         label="Polynomial")

n_points = 1

min_val = min(datapoints["polynom"])
min_index = datapoints["polynom"].index(min_val)
low_index = max(min_index-n_points, 0)
high_index = min(min_index + n_points + 1, len(datapoints["polynom"]))

dataset = datapoints["polynom"][low_index: high_index]
min_offsets = datapoints["offset"][low_index: high_index]

parabola = np.polyfit(min_offsets, dataset, 2)

def f(x):
    return (parabola[0] * (x ** 2)) + (parabola[1] * x) + parabola[2]

fitted_offset = -parabola[1] / (2 * parabola[0])
fitted_sens = f(fitted_offset)

message = "Minimum occurs at " + str.format('{0:.1f}', fitted_offset) + " \n"
message += "Minimum sensitivity is " + str.format('{0:.2f}', fitted_sens)

print message

plt.annotate(message, xy=(0.3, 0.9), xycoords="axes fraction")

mask = f(np.array(datapoints["offset"])) < max(datapoints["polynom"])

plt.plot(np.array(datapoints["offset"])[mask],
         f(np.array(datapoints["offset"])[mask]),
         linestyle="--", label="Parabola")

plt.xlabel("Offset (Days)")
plt.ylabel("Sensitivity (Arbitrary Units)")
plt.legend()
plt.tight_layout()
plot_path = output_path + "plots/combined_test_daifang_misalignment.pdf"
plt.savefig(plot_path)

print "Saving to", plot_path
