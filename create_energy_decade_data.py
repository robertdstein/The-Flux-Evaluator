import ConfigParser
import numpy as np
import os
from common import source_path, e_decades_dir

data_path = source_path + "data_configs/IC86.ini"

data_conf = ConfigParser.ConfigParser()
data_conf.read(data_path)

step = 0.25

# Loops over each season, creating an LLh Object.
# The data is randomised, and the LLh function is spline-fitted.
for season in data_conf.sections():

    exp_path = data_conf.get(season, "data_dir") + \
               data_conf.get(season, "exp_path")
    mc_path = data_conf.get(season, "data_dir") + \
              data_conf.get(season, "mc_path")

    exp = np.load(exp_path)
    mc = np.load(mc_path)

    print exp.dtype.names

    logE = exp["logE"]

    min_logE = int(np.min(logE))
    max_logE = int(np.max(logE)) + 1

    logE_bands = np.arange(min_logE, max_logE, step=step)

    print len(logE), np.min(logE), np.max(logE)
    print logE_bands

    total = 0

    veto = []

    for lower_lim in logE_bands:
        upper_lim = lower_lim + step

        mask = (logE >= lower_lim) & (logE < upper_lim)

        data = logE[mask]

        n = len(data)

        print "Between", lower_lim, "and", upper_lim, "we have", n

        total += n

        if n < 1000:
            print "Removing band", lower_lim, "due to insufficient statistics."
            veto.append(lower_lim)

    logE_bands = [x for x in logE_bands if x not in veto]

    print "In total, that's", total
    print "We are left with the following bands:"
    print logE_bands

    aw_dir = e_decades_dir + os.path.dirname(data_conf.get(season, "aw_path"))

    print aw_dir

    if not os.path.isdir(aw_dir):
        os.makedirs(aw_dir)

    print e_decades_dir

    cmd = "cp " + data_conf.get(season, "data_dir") +\
          data_conf.get(season, "aw_path") + "* " + aw_dir

    print os.system(cmd)

    for lower_lim in logE_bands:
        upper_lim = lower_lim + step

        name = str(lower_lim) + "<logE<" + str(upper_lim)

        new_paths = []

        for path in [exp_path, mc_path]:

            full_data = np.load(path)
            logE = full_data["logE"]
            mask = (logE >= lower_lim) & (logE < upper_lim)
            cut_data = full_data[mask]

            base = os.path.basename(path)
            new = "SUBSAMPLE_" + name + "_" + base

            output = e_decades_dir + new
            np.save(output, cut_data)

            print output

            new_paths.append(new)

        # aw = np.load()

        new_config_path = source_path + "data_configs/SUBSAMPLE_" + season + \
                          "_" + name + ".ini"

        new = ConfigParser.ConfigParser()

        section_name = season + "_" + name
        new.add_section(section_name)

        for (var, val) in data_conf.items(season):
            new.set(section_name, var, val)

        new.set(section_name, "exp_path", new_paths[0])
        new.set(section_name, "mc_path", new_paths[1])
        new.set(section_name, "data_dir", e_decades_dir)

        with open(new_config_path, "wb") as f:
            new.write(f)




