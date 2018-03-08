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

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action="store_true")
cfg = parser.parse_args()

user_dir = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"
save_dir = "/afs/ifh.de/user/s/steinrob/scratch/The-Flux-Evaluator__Output/"
file_name = "analysis_config/test_scaling_decay.ini"
pickle_names = "test_scaling_decay/length_"

test_configs_file = user_dir + file_name

Config = ConfigParser.ConfigParser()

lengths = np.linspace(5, 300, 60)
sim_length = 100

with open(test_configs_file, "w") as f:
    for i in lengths:
        name = pickle_names + str(i)
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
        Config.set(name, "CatName",
                   "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/" +
                   "Catalogue/catalogue00.npy")
        Config.set(name, "DataConfig", "fast.ini")

        sim = tm.return_integrated_light_curve(sim_params["t0"],
                                               "Decay", sim_params)

        alt = tm.return_integrated_light_curve(
            i + sim_params["t0"], "Decay", sim_params)

        maxk = 15 * ((1. / ((1 - sim) * alt)) + math.sqrt(i / 100))
        print i, maxk, sim, alt
        Config.set(name, "MaxK", maxk)

    Config.write(f)

# raw_input("prompt")

os.system("rm " + user_dir + "logs/*")

if cfg.submit:
    for section in Config.sections():
        # os.system(
        #     "python " + source_path + "RunLocal.py" +
        #     " -c " + section + " -f " + file_name + " -n 100 -s 5")

        RC.submit_to_cluster(10, section, file_name, ntrials=100, steps=20)
    RC.wait_for_cluster()

allfits = []

datapoints = {
    "length": [],
    "interpolation": [],
    "polynom": [],
    "analytic": []
}

for i in lengths:
    name = pickle_names + str(i)
    fits = MF.run(name, save_dir)
    if fits is not None:
        datapoints["length"].append(i)
        datapoints["interpolation"].append(fits["interpolation"])
        datapoints["polynom"].append(fits["polynom"])
        datapoints["analytic"].append(fits["mine"])

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
min_offsets = datapoints["length"][low_index: high_index]

parabola = np.polyfit(min_offsets, dataset, 2)

def f(x):
    return (parabola[0] * (x ** 2)) + (parabola[1] * x) + parabola[2]

fitted_offset = -parabola[1] / (2 * parabola[0])
fitted_sens = f(fitted_offset)

message = "Minimum occurs at " + str.format('{0:.1f}', fitted_offset) + " \n"
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
plt.gca().set_ylim(bottom=0)
plt.legend()
plt.tight_layout()
plot_path = save_dir + "plots/combined_test_scaling_decay.pdf"
plt.savefig(plot_path)

print "Saving to", plot_path
