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
    the equivalent skylab sensitivity.

    :param sindec: Sin(declination)
    :return: 7 year PS sensitivity at sindec
    """
    data = np.load(coenders_7year_sens_path)
    decs = np.array([x[0] for x in data])
    sens = np.array([x[2] for x in data])
    sens_ref = interp1d(np.sin(decs), sens)
    return sens_ref(sindec)

coenders_7year_sensitivity()