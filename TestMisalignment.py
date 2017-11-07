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

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

user_dir = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"
file_name = "test_misalignment.ini"
pickle_names = "test_misalignment/offset_"

test_configs_file = user_dir + file_name

Config = ConfigParser.ConfigParser()

offsets = np.linspace(-90, 90, 37)

with open(test_configs_file, "w") as f:
    for i in offsets:
        name = pickle_names + str(i)
        Config.add_section(name)
        Config.set(name, "UseEnergy", False)
        Config.set(name, "FitGamma", True)
        Config.set(name, "FixedGamma", 2)
        Config.set(name, "UseTime", True)
        Config.set(name, "SimTimeModel", "Box")
        Config.set(name, "SimTimeParameters", {"t0": 0, "length": 100})
        Config.set(name, "ReconTimeModel", "Box")
        Config.set(name, "ReconTimeParameters", {"t0": i, "length": 100})
        Config.set(name, "FitWeights", True)
        Config.set(name, "UseBox", True)
        Config.set(name, "CatName",
                   "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/" +
                   "Catalogue/catalogue00.npy")
        Config.set(name, "DataConfig", "fast.ini")
        Config.set(name, "MaxK", 30 * (100. / (100. - float(np.abs(i)))))

    Config.write(f)

if cfg.submit:
    for section in Config.sections():
        # os.system(
        #     "python " + user_dir + "RunLocal.py" +
        #     " -c " + section + " -f " + file_name + " -n 100 -s 10")

        RC.submit_to_cluster(10, section, file_name, ntrials=200, steps=20)
    RC.wait_for_cluster()

allfits = []

datapoints = {
    "offset": [],
    "interpolation": [],
    "polynom": [],
    "analytic": []
}

for i in offsets:
    name = pickle_names + str(i)
    fits = MF.run(name, user_dir)
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

fitted_offset = parabola[1] / (2 * parabola[0])
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
plot_path = user_dir + "plots/combined_test_misalignment.pdf"
plt.savefig(plot_path)

print "Saving to", plot_path
