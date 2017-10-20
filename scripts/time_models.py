"""
Contains all information regarding time PDF models. New MOdels can be created
by including their functions in this script.
"""
import numpy as np


def StepFunc(x, x_0):
    return 0.5 * (np.sign(x - x_0) + 1.)


def box_func_dict(parameters):
    box_dict = dict()
    box_dict["model_start"] = 0 + parameters["t0"]
    box_dict["model_end"] = parameters["length"] + parameters["t0"]
    box_dict["tot_norm"] = parameters["length"]
    return box_dict


def box_func(x, parameters):
    """Returns a Box Function, between the start and end time, relative to the

    :param x: time displacement from discovery
    :param parameters: Parameter Array
    :return: Box Function
    """
    box_dict = box_func_dict(parameters)

    f = (
        (StepFunc(x, box_dict["model_start"])
         - StepFunc(x, box_dict["model_end"])) / box_dict["tot_norm"])
    return f


def box_func_overlap(data_start, data_end, discovery_date, parameters):
    """Calculates the overlap of a box function light curve and a given
    detector run time

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
    a_dict = dict()
    a_dict["model_start"] = 0 + parameters["t0"]
    a_dict["model_end"] = parameters["length"] + parameters["t0"]
    a_dict["tot_norm"] = parameters["t_pp"] * (
        np.log(parameters["t_pp"] + a_dict["model_end"])
        - np.log(parameters["t_pp"] + a_dict["model_start"]))

    return a_dict


def AnalyticTimePDF(x,  parameters):
    a_dict = analytic_dict(parameters)

    t_max = parameters["length"]
    t = np.asarray(x)
    r = np.zeros_like(x)

    mask = np.logical_and(t >= 0., t < t_max)
    r[mask] = (1. + t[mask] / parameters["t_pp"]) ** -1.
    r[mask] = r[mask] / a_dict["tot_norm"]
    return r


def analytic_overlap(data_start, data_end, discovery_date, parameters):
    a_dict = dict()
    t_start = min(max(data_start, discovery_date + a_dict["model_start"]),
                  data_end)
    t_end = min(max(data_start, discovery_date + a_dict["model_end"]),
                  data_end)

    season_norm = a_dict["t_pp"] * (
        np.log(t_end + a_dict["t_pp"]) - np.log(t_start + a_dict["t_pp"]))

    return season_norm, a_dict["tot_norm"]


def return_light_curve(name, x, parameters):
    if name == "Box":
        return box_func(x, parameters)
    elif name == 'Decay':
        return AnalyticTimePDF(x, parameters)
    else:
        raise Exception("Model not found!")

def return_norms(name, data_start, data_end, discovery_date, parameters):
    if name == "Box":
        return box_func_overlap(data_start, data_end, discovery_date, parameters)
    elif name == "Decay":
        return analytic_overlap(data_start, data_end, discovery_date, parameters)
    else:
        raise Exception("Model not found!")

