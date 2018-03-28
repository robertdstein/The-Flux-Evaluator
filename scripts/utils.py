import numpy as np
import sys
from scipy.interpolate import interp1d
sys.path.append('..')
from common import coenders_7year_sens_path

# Factor to convert k-scale to flux and vice versa
factor = 1.e-12


def k_to_flux(x):
    return x * factor


def flux_to_k(x):
    return x / factor


def coenders_7year_sensitivity(sindec=0.0):
    """Interpolates between the saved values of the Stefan Coenders 7 year PS
    analysis sensitivity. Then converts given values for sin(declination to
    the equivalent skylab sensitivity. Adds values for Sindec = +/- 1,
    equal to nearest known value.

    :param sindec: Sin(declination)
    :return: 7 year PS sensitivity at sindec
    """
    data = np.load(coenders_7year_sens_path)
    sindecs = np.sin(np.array([x[0] for x in data]))
    sens = np.array([x[2] for x in data])

    # Extend range of sensitivity to +/- 1 through approximation,
    # by 1d-extrapolation of first/last pair

    sindecs = np.append(-1, sindecs)
    sindecs = np.append(sindecs, 1)

    lower_diff = sens[0] - sens[1]

    upper_diff = sens[-1] - sens[-2]

    sens = np.append(sens[0] + lower_diff, sens)
    sens = np.append(sens, sens[-1] + upper_diff)

    sens_ref = interp1d(sindecs, sens)

    return sens_ref(sindec)

# def TFE_IC40_sensitivity(sindec=0.0):
#     """Interpolates between the saved values of the point source sensitivity,
#     calculated with this code, using only IC40 data. Then converts given values
#     for sin(declination) to the equivalent IC40 sensitivity.
#
#     :param sindec: Sin(declination)
#     :return: IC40 PS sensitivity at sindec
#     """
#     data = np.load(IC40_sens_path)
#     decs = np.array([x[0] for x in data])
#     sens = np.array([x[2] for x in data])
#     sens_ref = interp1d(np.sin(decs), sens)
#     return sens_ref(sindec)