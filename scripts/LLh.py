import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as scp
import scipy.optimize
import copy

import numexpr
numexpr.set_num_threads(1)

import time
from scipy import interpolate

import numpy.lib.recfunctions

from scripts.Injector import Injector
from scripts.PDF import PDF
from scripts.RandomTools import RandomTools

import resource


class LLh(PDF, Injector):

    """A class for the Log Likelihood. Handles loading the data, preparing
    randomised datasets based on the data, and calculates the Log Likelihood
    for a given parameter set.

    """

    def __init__(self, **kwargs):
        # Bins for energy (Tev?)
        self.EnergyBins = np.linspace(1., 10., 40 + 1)
        # Bins for sin declination (not evenly spaced)
        self.sinDecBins = np.unique(np.concatenate([
                                np.linspace(-1., -0.9, 2 + 1),
                                np.linspace(-0.9, -0.2, 8 + 1),
                                np.linspace(-0.2, 0.2, 15 + 1),
                                np.linspace(0.2, 0.9, 12 + 1),
                                np.linspace(0.9, 1., 2 + 1),
                                ]))

        # Provides grid values for spectral index (gamma)
        self.GammaGrid = np.linspace(0, 5., 20 + 1)
        # Sets precision
        self.precision = .1

        # Sets weights to default value of nan
        self.SeasonWeight = np.nan
        self.AcceptanceWeight = np.nan

        # Sets "_..." to False
        self._ReturnInjectorNExp = False

        # Initialises empty caches
        self._g1 = np.nan
        self._w_cache = np.nan
        self.w_cache = dict()
        self.spline_cache = dict()
        self.w_cache_BG = dict()
        self.w_cache_Sig = dict()

        # Produces a set (i.e no duplicates) of datapoints for gamma
        # This is best on 33 points between 0.9 and 4.1
        # Each point is modified by _around(i)
        # Useful for different precisions, where rounding errors might
        # otherwise lead to duplicates in set
        self.GammaSupportPoints = set(
            [self._around(i) for i in np.linspace(0.9, 4.1, 30 + 3)])

        self.N_all = np.nan
        self.n_select = np.nan
        self._exp = np.nan

        # Sets start and end date for data taking (typo irrelevant)
        self.DataStart = kwargs['StartDataTakingMJD']
        self.DataEnd = kwargs['EndDataTakingMJD']

        # Sets the livetime, and the total Season length
        self.Livetime = kwargs['Livetime']
        self.SeasonTimeSpan = self.DataEnd - self.DataStart

        # Loads "Experimental" and "Monte Carlo" data sets (in .npy format)
        self.exp = np.load(kwargs['ExpPath'])
        self.mc = np.load(kwargs['MCPath'])

        # ***********************************************************************
        # Loads sources from the given source path
        # Adds field to numpy array
        self._sources2 = np.load(kwargs['SourcePath'])
        self._sources = np.lib.recfunctions.append_fields(
            self._sources2, 'TimeNorm', data=np.ones_like(self._sources2['ra']))
        del self._sources2

        # Interpolates (in 2D) a function for the Gamma and declination
        try:
            dec_bins = np.load(
                kwargs['AcceptanceWeightPath'] + '_bins_dec.npy')
            gamma_bins = np.load(
                kwargs['AcceptanceWeightPath'] + '_bins_gamma.npy')
            values = np.load(kwargs['AcceptanceWeightPath'] + '_values.npy')
            f = interpolate.interp2d(
                dec_bins, gamma_bins, values, kind='linear')
            self.AcceptanceFitFunc = f
            del f
        except:
            print('No Acceptance Files')

        # Sets Use Energy, Fit Gamma and Fixed Gamma
        self.UseEnergy = kwargs['UseEnergy']
        self.FitGamma = kwargs['FitGamma']
        self.FixedGamma = kwargs['FixedGamma']

        # Toggles using a box for sources to reduce runtime
        self.UseBox = kwargs['UseBox']

        # Sets Fit Weights (True/False for fitting the weights of each)
        self.FitWeights = kwargs['FitWeights']

        # Sets UseTime and Time Model (Box/BoxPre/Decay)
        self.UseTime = kwargs['UseTime']
        self.TimeModel = kwargs['TimeModel']

        # If fitting for time
        if self.UseTime:
            # **************************************************************************************************************************
            # Adds a 'Discovery Delay' (10 days), which is subtracted from
            # discovery date
            self.DiscDelay = 10.
            self.sources['discoverydate_mjd'] = (
                self.sources['discoverydate_mjd'] - self.DiscDelay)

            if self.TimeModel == 'Box':
                self.TimeBoxLength = kwargs['TimeBoxLength']

            if self.TimeModel == 'BoxPre':
                self.TimeBoxLength = kwargs['TimeBoxLength']

            if self.TimeModel == 'Decay':
                # Sets Model Timescale (p-p), and decay mode length
                self.Model_tpp = kwargs['Model_tpp']
                self.DecayModelLength = kwargs['DecayModelLength']
                self.Model_Length = self.DecayModelLength

        # ****************************************************************************************************************************************
        # Sets smear Injection and MissTiming
        # Inject with wrong weight
        # MT Inject with scaled time window (longer/shorter
        self.SmearInjection = kwargs['SmearInjection']
        self.MissTiming = kwargs['MissTiming']

        # Creates empty dictionary, filled by Injector.find_and_apply_band_mask
        self.InjectionBandMask = dict()

        # If Injection Gamma is not given, sets it to 2.
        # Then gives Injection Weights. Otherwise loads injection gamma.
        if 'InjectionGamma' in kwargs.keys():
            self.InjectionGamma = kwargs['InjectionGamma']
        else:
            self.InjectionGamma = 2.
            self.WeightsInject = self.get_weights(self.mc, self.InjectionGamma)
        print 'Injection Spectrum', self.InjectionGamma

        # Creates a weight cache with default value nan
        self.WeightCache = np.nan

        # Makes sure unblinding is false
        self.Unblind = False

        # If not unblinded, times are shuffled. Otherwise gives a warning.
        if self.Unblind is False:
            np.random.shuffle(self.exp['timeMJD'])
        else:
            print('WARNING: Running in Unblinding Mode')

        # ***********************************************************************
        # Sets a default of Bootstrap=False, if not provided in init
        # Statistical Method to evaluate uncertanties
        # Test sensitivity
        self.BootStrap = False
        try:
            self.BootStrap = kwargs['BootStrap']
        except:
            pass

# ==============================================================================
# Running Part
# ==============================================================================
    # Designed for backwards compatibility,
    # so that self.sources returns self._sources
    @property
    def sources(self, ):
        """Gives self._sources"""
        return self._sources

    def prepare_fake_data_set_and_evaluate_pdf(self, k, ):
        """Prepares a scrambled/randomised dataset.

        :param k: Flux scale
        :return: Scrambled dataset
        """
        # Initialises empty dictionaries and nan-values
        self.w_cache = dict()
        self.w_cache_Sig = dict()
        self._ev = np.nan
        self._ev_B = np.nan
        self._ev_S = np.nan
        self.EnergyWeightCache = np.nan
        self.SoB = np.nan

        # If UseBox, only selects events in Box,
        # In either case, 'blinds' the experimental data
        if self.UseBox is True:
            self.w_cache_BG = dict()
            self._ev = self.selects_events_in_band(self.sources[0], self.exp)
            self._ev = self.scramble_exp_data(self._ev)
            self._ev = self.select_events_in_box(self.sources[0], self._ev, )
        if self.UseBox is False:
            self._ev = self.scramble_exp_data(self.exp)

        # Produces signal events
        self.sig_events = self.generate_sig_events(self.sources, self.mc, k, )

        # ***********************************************************************
        # If UseEnergy, then generates energy weights?
        if self.UseEnergy is True:
            self.generate_sig_weight_dict_for_all_gamma(self.sig_events)
            self.generate_bkg_weight_dict_for_all_gamma(self._ev)

        # Creates a combined dataset
        self._ev = np.concatenate((self._ev, self.sig_events))
        self.N_all = len(self.exp) + len(self.sig_events)
        self.N = len(self._ev)

        # ***********************************************************************
        # What does this do?
        self.evaluate_bkg()
        self.evaluate_sig()

        del self._ev
        self.SoB = self.SoB / self._ev_B
        del self._ev_B

        # If UseEnergy but not FitGamma, set gamma to FixedGam
        if self.UseEnergy is True:
            if self.FitGamma is False:
                gamma = self.FixedGamma
                self.EnergyWeightCache = np.append(
                    self.weight_fast(gamma, self.w_cache_BG),
                    self.weight_fast(gamma, self.w_cache_Sig))

    def init_everything_for_multiple_trials(self, ):
        """Initialises all attributes which need to be calculated once for
        all trials. Scrambles/randomises the experimental data, creates a
        spatial pdf for the background and assigns weight distances. Also
        produces optional Energy PDFs and Time weights, if required.
        """
        # Produces scrambled experimental data set
        self._ev = self.exp
        self._ev = self.scramble_exp_data(self._ev)

        # Finds a spatial PDF for the background, based on the experimental
        # Sin Declination distribution
        bckg_spline_space = self.create_space_bkg_pdf(self._ev)
        self.bckg_spline_space = bckg_spline_space

        # Assigns a weight to each source, equal to 1/(r^2) for distance r
        self.sources['weight_distance'] = self.sources['distance']**(-2.)

        # If accounting for energy, produces Energy PDFs
        if self.UseEnergy is True:
            print('Initialising Energy PDFs')
            self.generate_spline_dict_for_all_gamma(self.exp, self.mc)
            self.generate_bkg_weight_dict_for_all_gamma(self._ev)

        # If using time, calculates Time weights for the source
        if self.UseTime is True:
            self.compute_source_weights_time()
            self.init_random_generator_pdf()

    def ProduceLLhFunction(self, ):
        """ Creates a function to give f(x) = -1 * Test Statistic (x),
        and returns it.
        
        :return: function
        """
        f = lambda x: - self.test_stat_new_llh_fit_weights(x, )

        return f

# ==============================================================================
# End Running Part
# ==============================================================================

    def test_stat_new_llh_fit_weights(self, params, ):
        """Calculates the test statistic, given the parameters. Uses numexpr
        for faster calculations.

        :param params: Parameters from minimisation
        :return: Test Statistic
        """
        N_all = self.N_all
        N = self.N
        SoB = self.SoB

        # If using Energy, finds gamma and the remaining parameters
        # Calculates the Signal/background ratio for
        if self.UseEnergy is True:
            n = params[:-1]
            gamma = params[-1]
            w = np.append(self.weight_fast(gamma, self.w_cache_BG),
                          (self.weight_fast(gamma, self.w_cache_Sig)))

        if self.UseEnergy is False:
            n = params[:]
            w = 1.

        n *= self.SeasonWeight
        b = self.sources['weight'][:, None]

        # If not fitting weights, use array. Otherwise, use matrix.
        if self.FitWeights is False:
            y = np.sum(numexpr.evaluate('(b * SoB)'), axis=0)
            a = numexpr.evaluate('n * (w*y-1.)')
            del w
        if self.FitWeights is True:
            y = numexpr.evaluate('SoB')
            c = np.matrix(numexpr.evaluate('w*SoB-1.'))
            del w
            a = np.squeeze(np.asarray(n * c))
            del c

        llh_value = numexpr.evaluate('log1p(a / N_all)')

        # Make sure to delete variables from memory
        del SoB
        del y
        del a

        if self.UseBox is True:
            llh_value = llh_value.sum() + (N_all - N) * np.log1p(-n / N_all)
        else:
            llh_value = llh_value.sum()

        # Definition of test statistic
        return 2. * llh_value