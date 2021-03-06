"""
Contains all information regarding time PDF models. New Models can be created
by including their functions in this script.
"""
import numpy as np


def step_func(x, x_0):
    return 0.5 * (np.sign(x - x_0) + 1.)


def box_func_dict(parameters):
    """Unpacks parameters, and creates a dictionary containing the  model
    start/end times, as well as the total normalisation

    :param parameters: Parameters for box model ("t0" is offset of neutrino
    lightcurve relative to discovery date, and "length" is width of emission
    box)
    :return: Dictionary
    """
    box_dict = dict()
    box_dict["model_start"] = 0. + float(parameters["t0"])
    box_dict["model_end"] = (
        float(parameters["length"]) + float(parameters["t0"]))
    box_dict["tot_norm"] = float(parameters["length"])
    return box_dict


def box_func(t, parameters):
    """Returns a Box Function, between the start and end time, relative to the
    discovery date.

    :param t: time displacement from discovery
    :param parameters: Parameters for box model ("t0" is offset of neutrino
    lightcurve relative to discovery date, and "length" is width of emission
    box)
    :return: Box Function value at x
    """
    box_dict = box_func_dict(parameters)
    val = (
        (step_func(t, box_dict["model_start"])
         - step_func(t, box_dict["model_end"])
         )
        / box_dict["tot_norm"])
    return val


def box_func_integrated(t, parameters):
    """Returns the ratio of the integral of the Box Function from -infinity
    to t, and -infinity to +infinity.

    :param t: time
    :param parameters: Parameters for box model ("t0" is offset of neutrino
    lightcurve relative to discovery date, and "length" is width of emission
    box)
    :return: Fraction of total integral covered up to t
    """
    t = np.asarray(t)
    box_dict = box_func_dict(parameters)
    norm = box_dict["tot_norm"]

    r = np.zeros_like(t)
    r[t > box_dict["model_end"]] = np.ones_like(t[t > box_dict["model_end"]])
    mask = np.logical_and(
        t >= box_dict["model_start"],
        t <= box_dict["model_end"])

    # r = np.zeros_like(t)
    # r[t > parameters["length"]] = np.ones_like(t[t > parameters["length"]])
    # mask = np.logical_and(
    #     t >= 0,
    #     t <= parameters["length"])

    r[mask] = (t[mask]-box_dict["model_start"]) / norm

    # print t, t-box_dict["model_start"], parameters["length"], r, mask

    return r


def box_func_overlap(data_start, data_end, discovery_date, parameters):
    """Calculates the overlap of a box function light curve and a given
    detector run time.

    :param data_start: Start of data taking (MJD)
    :param data_end: End of data taking (MJD)
    :param discovery_date: Date of Source Discovery (MJD)
    :param parameters: Parameters for box model ("t0" is offset of neutrino
    lightcurve relative to discovery date, and "length" is width of emission
    box)
    :return: season_norm (Overlap of datataking and given Box Model) and also
    tot_norm (the overall width of the given box model)

    """
    box_dict = box_func_dict(parameters)
    t_start = min(
        max(data_start, discovery_date + box_dict["model_start"]), data_end)
    t_end = min(
        max(data_start, discovery_date + box_dict["model_end"]), data_end)

    season_norm = t_end-t_start
    return season_norm, box_dict["tot_norm"]


def analytic_dict(parameters):
    """Unpacks parameters, and creates a dictionary containing the  model
    start/end times, as well as the total normalisation

    :param parameters: Parameters for analytic model ("t0" is offset of neutrino
    lightcurve relative to discovery date, "length" is width of emission
    box and "t_pp" is the time constant.)
    :return: Dictionary
    """
    a_dict = dict()
    a_dict["model_start"] = 0 + parameters["t0"]
    a_dict["model_end"] = parameters["length"] + parameters["t0"]
    a_dict["tot_norm"] = parameters["t_pp"] * (
        np.log(1 + parameters["length"] / parameters["t_pp"]))

    return a_dict


def analytic_func(t, parameters):
    """Analytic exponential Decay Model, beginning at t0, ending at
    t0+length, to describe the neutrino emission.

    :param t: Time
    :param parameters: Parameters for analytic model ("t0" is offset of neutrino
    lightcurve relative to discovery date, "length" is width of emission
    box and "t_pp" is the time constant.)
    :return: Value of f at t
    """
    a_dict = analytic_dict(parameters)

    t = np.asarray(t)
    r = np.zeros_like(t)

    mask = np.logical_and(
        t >= a_dict["model_start"], t < a_dict["model_end"])
    r[mask] = (1. + (t[mask] - parameters["t0"]) / parameters["t_pp"]) ** -1.
    r[mask] = r[mask] / a_dict["tot_norm"]
    return r

def analytic_integrated(t, parameters):
    """Returns the ratio of the integral of the Analytic expoential decay model
    from -infinity to t, and -infinity to +infinity.

    :param t: Time
    :param parameters: Parameters for analytic model ("t0" is offset of neutrino
    lightcurve relative to discovery date, "length" is width of emission
    box and "t_pp" is the time constant.)
    :return: Fraction of total integral covered up to t
    """
    a_dict = analytic_dict(parameters)
    t = np.asarray(t)
    r = np.ones_like(t)

    mask = np.logical_and(
        t >= a_dict["model_start"], t < a_dict["model_end"])
    r[mask] = np.log(1. + (t[mask] - parameters["t0"]) / parameters["t_pp"])
    r[mask] = r[mask] / np.log(1 + parameters["length"]/parameters["t_pp"])
    r[t < a_dict["model_start"]] = 0.
    return r


def analytic_overlap(data_start, data_end, discovery_date, parameters):
    a_dict = analytic_dict(parameters)
    t_start = min(max(data_start, discovery_date + a_dict["model_start"]),
                  data_end)
    t_end = min(max(data_start, discovery_date + a_dict["model_end"]),
                  data_end)

    season_norm = a_dict["t_pp"] * (
        np.log(t_end + a_dict["t_pp"]) - np.log(t_start + a_dict["t_pp"]))

    return season_norm, a_dict["tot_norm"]


def return_light_curve(t, name, parameters):
    """Returns the light curve function corresponding to a given name

    :param t: Time
    :param name: Name of model
    :param parameters: Dictionary containing all necessary parameters for
    the given model
    :return: The light curve function
    """
    if name == "Box":
        return box_func(t, parameters)
    elif name == 'Decay':
        return analytic_func(t, parameters)
    else:
        raise Exception("Model not found!")


def return_integrated_light_curve(t, name, parameters):
    if name == "Box":
        return box_func_integrated(t, parameters)
    elif name == "Decay":
        return analytic_integrated(t, parameters)
    else:
        raise Exception("Model not found!")


def return_norms(name, data_start, data_end, discovery_date, parameters):
    if name == "Box":
        return box_func_overlap(data_start, data_end, discovery_date, parameters)
    elif name == "Decay":
        return analytic_overlap(data_start, data_end, discovery_date, parameters)
    else:
        raise Exception("Model not found!")

if __name__ == '__main__':
    # If this script is simply run directly, it will produce plots of both
    # the time PDFs and their integral functions

    import os

    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure()

    ax1 = plt.subplot(311)
    ax1.set_title("Time PDF")
    ax2 = plt.subplot(312)
    ax2.set_title("Normalised Integral function of Time PDF")
    ax3 = plt.subplot(313)
    ax3.set_title("Approximated Integral of Time PDF")

    t_range = np.linspace(-150, 150, 301)
    parameters ={
        "t0": -100,
        "length" : 100,
        "t_pp" : 1.0,
    }

    for model in ["Box", "Decay"]:
        ax1.plot(
            t_range,
            return_light_curve(t_range, model, parameters),
            label=model)
        ax2.plot(
            t_range,
            return_integrated_light_curve(t_range, model, parameters),
            label=model)

        approx_integral = []
        val = 0.0

        for t in t_range:
            val += return_light_curve(t, model, parameters)
            approx_integral.append(val)

        approx_integral = np.array(approx_integral)
        approx_integral = approx_integral/approx_integral[-1]

        ax3.plot(t_range, approx_integral, label="model")

    fig.set_size_inches(8, 12)

    root = "/afs/ifh.de/user/s/steinrob/Desktop/python/The-Flux-Evaluator/"
    path = root + "plots/time_PDFs/all_PDFs.pdf"

    if not os.path.isdir(os.path.dirname(path)):
        os.mkdir(os.path.dirname(path))

    plt.savefig(root + "plots/time_PDFs/all_PDFs.pdf")
    plt.close()

    graphs = {
        "perfect": {
            "t0": -100,
            "length": 100,
            "t_pp": 2.0 },
        "misaligned": {
            "t0": -50,
            "length": 100,
            "t_pp": 2.0 },
        "scaled": {
            "t0": -100,
            "length": 50,
            "t_pp": 2.0},
        "both": {
            "t0": -50,
            "length": 200,
            "t_pp": 2.0 }
    }

    for sim_model in ["Box", "Decay"]:

        for recon_model in ["Box", "Decay"]:

            for name, parameters in graphs.iteritems():
                fig = plt.figure()

                ax1 = plt.subplot(211)
                ax2 = plt.subplot(212, sharex=ax1)

                x = np.linspace(-150, 150, 3001)

                y1 = return_light_curve(t_range, sim_model, graphs["perfect"])
                y2 = return_light_curve(t_range, recon_model, parameters)

                ax1.fill_between(t_range, 0, y1, facecolor="blue")

                ax2.fill_between(t_range, 0, y2, facecolor="red")

                ax1.set_ylim(0, 1.2 * max(y1 + y2))
                ax2.set_ylim(0, 1.2 * max(y1 + y2))
                fig.subplots_adjust(hspace=0)
                ax2.invert_yaxis()
                plt.xlabel("Time (days)")
                ax2.set_ylabel("Recon PDF")
                ax1.set_ylabel("Emission PDF")

                ax1.get_xaxis().set_visible(False)

                fig.set_size_inches(8, 6)
                plt.savefig(
                    root + "plots/time_PDFs/" + sim_model + "_" +
                    recon_model + "_" + name + ".png")
                plt.close()

