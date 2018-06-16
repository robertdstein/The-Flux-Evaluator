import numpy as np
from icecube import astro

import scipy

import numexpr
numexpr.set_num_threads(1)

import scipy.interpolate

import time_models as tm

class PDF():

    """Class for handling PDFs and weighting.

    """

    def __init__(self, ):
        pass

# ==============================================================================
# Tool Part
# ==============================================================================

    def get_weights(self, mc, gamma, ):
        """Returns a weight array given by
        OneWeights * Energy_weights.

        :param mc: Monte Carlo
        :param gamma: Spectral Index
        :return: Weights
        """
        # Uses numexpr for faster processing
        ow = mc['ow']
        trueE = mc['trueE']
        weights = numexpr.evaluate('ow * trueE **(-gamma)')
        return weights

    def _around(self, value):
        """Produces an array in which the precision of the value
        is rounded to the nearest integer. This is then multiplied
        by the precision, and the new value is returned.

        :param value: value to be processed
        :return: value after processed
        """
        return np.around(float(value) / self.precision) * self.precision

    def empty_cache(self, ):
        """Empties the cache by setting _g1/_w_cache to default value (nan). """
        self._g1 = np.nan
        self._w_cache = np.nan

    def selects_events_in_band(self, source, data, ):
        """Selects events in a given sin declination band
        The width of the band is defined as 20 degrees.
        This function, if followed by select_events_in_box(),
        can be used to produce a square on the celestial sphere.

        :param source: Individual Source
        :param data: Data to be considered
        :return: data (modified)
        """
        dec_bandwidth = np.deg2rad(10.)
        min_dec = max(-np.pi / 2., source['dec'] - dec_bandwidth)
        max_dec = min(np.pi / 2., source['dec'] + dec_bandwidth)
        band_mask = ((data['sinDec']) > np.sin(min_dec)) & \
            ((data['sinDec']) < np.sin(max_dec))
        data_band = data[band_mask]
        N_all = len(data)
        n_select = len(data_band)
        self.N_all = N_all
        self.n_select = n_select
        return data_band

    def select_events_in_box(self, source, data, ):
        """Selects events in a given ra band.
        In combination with the declination band of selects_events_in_band(),
        this produces a square on the surface of the celestial sphere.
        The width of the ra band is variable to account for spherical geometry.

        :param source: Individual Source
        :param data: Data to be considered
        :return: data (modified)
        """
        dec_bandwidth = np.deg2rad(10.)
        min_dec = max(-np.pi / 2., source['dec'] - dec_bandwidth)
        max_dec = min(np.pi / 2., source['dec'] + dec_bandwidth)
        cosFact = np.amin(np.cos([min_dec, max_dec]))
        dPhi = np.amin([2. * np.pi, 2. * dec_bandwidth / cosFact])
        ra_dist = np.fabs(
            (data["ra"] - source['ra'] + np.pi) % (2. * np.pi) - np.pi)
        mask = ra_dist < dPhi / 2.
        data_box = data[mask]
        return data_box

# ==============================================================================
# Energy Part
# ==============================================================================

    def weight(self, data, exp, mc, gamma):
        """Calculates the Signal/Background function.
        Uses Finite Difference Methods to calculate 1st/2nd derivatives.
        Uses a Taylor series to estimate S(gamma).

        :param data:
        :param exp: Experimental data (background)
        :param mc: Monte Carlo data (signal)
        :param gamma: Spectral Index
        :return: Estimated value for S(gamma)
        """
        # Sets g (gamma), and the precision/stepsize (dg or h)
        g1 = self._around(gamma)
        dg = self.precision

        if g1 in self.w_cache.keys():
            # If derivatives already calculated, proceed.
            S1 = self.w_cache[g1]["S1"]
            a = self.w_cache[g1]["a"]
            b = self.w_cache[g1]["b"]
        else:
            # Add/subtract a step from g1 (g+h) and (g-h)
            g0 = self._around(g1 - dg)
            g2 = self._around(g1 + dg)

            # Calculate Spline stuff (sig/background 2D histograms)
            # Gives y(g), y(g+h), y(g-h)
            S0 = self.generate_spline_stuff(exp, mc, gamma=g0).ev(
                data['logE'], data['sinDec'])
            S1 = self.generate_spline_stuff(exp, mc, gamma=g1).ev(
                data['logE'], data['sinDec'])
            S2 = self.generate_spline_stuff(exp, mc, gamma=g2).ev(
                data['logE'], data['sinDec'])

            # Calculates Derivatives Using Finite Difference Methods (1st/2nd)
            a = (S0 - 2. * S1 + S2) / (2. * dg ** 2)
            b = (S2 - S0) / (2. * dg)

            # Stores information in weight cache
            self._g1 = g1
            self._w_cache = np.zeros(
                (len(data),), dtype=[("S1", np.float), ("a", np.float),
                                     ("b", np.float)])
            self._w_cache["S1"] = S1
            self._w_cache["a"] = a
            self._w_cache["b"] = b

            self.w_cache[g1] = dict()
            self.w_cache[g1]["S1"] = S1
            self.w_cache[g1]["a"] = a
            self.w_cache[g1]["b"] = b

        # Returns the Taylor series value for S(gamma) around g1
        val = np.exp(a * (gamma - g1) ** 2 + b * (gamma - g1) + S1)
        return val

    def weight_fast(self, gamma, w_cache):
        """Quickly estimates the value of Signal/Background for Gamma.
        Uses precalculated values for first and second derivatives.
        Uses a Taylor series to estimate S(gamma).

        :param gamma: Spectral Index
        :param w_cache: Weight cache
        :return: Estimated value for S(gamma)
        """
        g1 = self._around(gamma)
        dg = self.precision

        g0 = self._around(g1 - dg)
        g2 = self._around(g1 + dg)

        # Uses Numexpr to quickly estimate S(gamma)
        S0 = w_cache[g0]
        S1 = w_cache[g1]
        S2 = w_cache[g2]

        val = numexpr.evaluate(
            "exp((S0 - 2.*S1 + S2) / (2. * dg**2) * (gamma - g1)**2" + \
            " + (S2 -S0) / (2. * dg) * (gamma - g1) + S1)"
        )
        return val

    def weight_fast2(self, gamma, w_cache):
        """Calculates the Signal/Background function.
        Using precalculated values for the S(gamma) at different points.
        Uses Finite Difference Methods to calculate 1st/2nd derivatives.
        Uses a Taylor series to estimate S(gamma).

        :param gamma: Spectral Index
        :param w_cache: Weight cache
        :return: Estimated value for S(gamma)
        """
        g1 = self._around(gamma)
        dg = self.precision

        g0 = self._around(g1 - dg)
        g2 = self._around(g1 + dg)

        S0 = w_cache[g0]
        S1 = w_cache[g1]
        S2 = w_cache[g2]

        a = (S0 - 2. * S1 + S2) / (2. * dg ** 2)
        b = (S2 - S0) / (2. * dg)

        val = np.exp(a * (gamma - g1) ** 2 + b * (gamma - g1) + S1)
        return val

    def generate_spline_dict_for_all_gamma(self, exp, mc, ):
        """Loops over "generate_spline_stuff" for each Gamma point.
        Returns the interpolated 2D spline of the log(sig/bkg) ratios.

        :param exp: Experimental data (background)
        :param mc: Monte Carlo data (signal)
        """
        for gamma in self.GammaSupportPoints:
            if gamma in self.spline_cache.keys():
                pass
            else:
                self.spline_cache[gamma] = self.generate_spline_stuff(
                    exp, mc, gamma=gamma)

    def generate_bkg_weight_dict_for_all_gamma(self, exp):
        """Loops over each gamma points, and check if the point has been
        added to the Background Weight Cache. If not, loads the associated
        2D spline function. Using the array of experimental data, it converts
        each event to a Log(Sig/Bkg) value, using Log(Energy) and sin(Dec)
        with the fitted function.

        :param exp: Experimental data (Background)
        """
        for gamma in self.GammaSupportPoints:
            if gamma in self.w_cache_BG.keys():
                pass
            else:
                g1 = self._around(gamma)
                self.w_cache_BG[gamma] = self.spline_cache[g1].ev(
                    exp['logE'], exp['sinDec'])

    def generate_sig_weight_dict_for_all_gamma(self, sig_events):
        """Loops over each gamma points, and check if the point has been
        added to the Signal Weight Cache. If not, loads the associated
        2D spline function. Using the array of experimental data, it converts
        each event to a Log(sig/Bkg) value, using Log(Energy) and sin(Dec)
        with the fitted function.
        """
        for gamma in self.GammaSupportPoints:
            if gamma in self.w_cache_Sig.keys():
                pass
            else:
                g1 = self._around(gamma)
                self.w_cache_Sig[gamma] = self.spline_cache[g1].ev(
                    sig_events['logE'], sig_events['sinDec'])

    def generate_spline_stuff(self, exp, mc, gamma=2.):
        """Produces normalised histograms for signal and background.
        Produces a histogram of the ratio of sig/bkg. Interpolates, with a
        spline function, using the log of this 2D distribution.

        :param exp: Experimental data (Background)
        :param mc: Monte Carlo data (signal)
        :param gamma: Spectral Index
        :return: 2D Spline function
        """
        # Produces arrays from histograms for sin(dec) vs logE distribution
        # The first array is experimental (to simulate bkg)
        # The second array is Monte Carlo (to simulate signal)
        h_bkg = self.create_energy_hist_for_splines(
            exp['sinDec'], exp['logE'], weights=np.ones_like(exp['logE']))
        h_sig = self.create_energy_hist_for_splines(
            np.sin(mc['trueDec']), mc['logE'], weights=self.get_weights(mc,
                                                                        gamma))

        # Produces an array containing True if x> 0, False otherwise
        domain_bkg = h_bkg > 0.
        domain_sig = h_sig > 0.

        # Creates an array of ones as the default ratio
        ratio = np.ones_like(h_bkg, dtype=np.float)
        # Bitwise Addition giving an Truth Array
        # Returns True if True in both Sig and Bkg, otherwise False
        mask = domain_bkg & domain_sig
        # Calculates the ratio sig/bkg for those entries True in Mask array
        ratio[mask] = (h_sig[mask] / h_bkg[mask])

        # Finds the minimum ratio
        min_ratio = np.amax(ratio)
        # Where true in sig and false in bkg, sets ratio to minimum ratio
        np.copyto(ratio, min_ratio, where=domain_sig & ~domain_bkg)

        # Sets Bin centers, and order of spline (for x and y)
        sin_dec_bin_center = (self.sinDecBins[:-1] + self.sinDecBins[1:]) / 2.
        log_e_bin_center = (self.EnergyBins[:-1] + self.EnergyBins[1:]) / 2.
        log_e_order = 2

        # Fits a 2D spline function to the log of ratio array
        # This is 2nd order in both dimensions
        spline = scipy.interpolate.RectBivariateSpline(
            log_e_bin_center, sin_dec_bin_center, np.log(ratio),
            kx=log_e_order, ky=log_e_order, s=0)
        self.EnergyPDFSplineInterpolated = spline

        return spline

    def create_energy_hist_for_splines(self, sin_dec, log_e, weights, ):
        """Produces normalised histograms of distribution of variables
        Sin(Dec) and Log(Energy), for use in Spline interpolation.

        :param sin_dec: Sin Declination
        :param log_e: Log(Energy)
        :param weights: Weights for distribution
        :return: Normalised histogram values
        """
        energy_bins = self.EnergyBins
        sin_dec_bins = self.sinDecBins
        # Produces the histograms
        h, binedges = np.histogramdd(
            (log_e, sin_dec), bins=(energy_bins, sin_dec_bins), weights=weights)
        ndim = h.ndim

        # Normalises histogram
        norms = np.sum(h, axis=ndim - 2)
        norms[norms == 0.] = 1.
        h /= norms

        return h

# ==============================================================================
# Background Functions
# ==============================================================================

    def create_space_bkg_pdf(self, exp, ):
        """Creates the spatial PDF for background.
        Generates a histogram for the exp. distribution in sin declination.
        Fits a spline function to the distribution, giving a spatial PDF.
        Returns this spatial PDF.

        :param exp: Experimental data (background)
        :return: Background spline function
        """
        sin_dec_bins = self.sinDecBins
        sin_dec_range = (np.min(sin_dec_bins), np.max(sin_dec_bins))
        hist, bins = np.histogram(
            exp['sinDec'], density=True, bins=sin_dec_bins, range=sin_dec_range)

        bins = np.concatenate([bins[:1], bins, bins[-1:]])
        hist = np.concatenate([hist[:1], hist, hist[-1:]])

        #In case of low statistics, fill zero-bins with minimum value!!!

        histarr = np.array(hist)
        mask = histarr > 0
        # print histarr, np.log(histarr[mask]), np.sum(np.log(histarr[mask]))

        minval = np.min(histarr[mask])
        histarr[~mask] = minval

        # print histarr, np.log(histarr), np.sum(np.log(histarr))
        hist = list(histarr)
        # print hist
        #
        # raw_input("prompt")


        bkg_spline = scipy.interpolate.InterpolatedUnivariateSpline(
                                (bins[1:] + bins[:-1]) / 2.,
                                np.log(hist), k=2)
        return bkg_spline

    def evaluate_bkg(self, ):
        """Loads the background spline function, calculates the background
        PDF function, and stores it.
        """
        data = self._ev
        bckg_spline = self.bckg_spline_space
        B_space = 1. / 2. / np.pi * np.exp(bckg_spline(data["sinDec"]))

        # If UseTime, scale the background with the runtime
        if self.UseTime is True:
            B_time = np.ones_like(B_space) / (self.DataEnd - self.DataStart)
        else:
            B_time = np.ones_like(B_space)

        self._ev_B = B_space * B_time

# ==============================================================================
# Signal Functions
# ==============================================================================

    def space_pdf_signal(self, source, data, ):
        """Calculates the angular distance between the source_path and the data.
        Uses a Gaussian PDF function, centered on the source_path.
        Returns the value of the Gaussian at the given distances.

        :param source: Single Source
        :param data: Dataset to be compared
        :return: Array of Gaussian values
        """
        distance = astro.angular_distance(
            data['ra'], data['dec'], source['ra'], source['dec'])
        space_term = (1. / (2. * np.pi * data['sigma'] ** 2.) *
                      np.exp(-0.5 * (distance / data['sigma']) ** 2.))
        return space_term

    def time_pdf_signal(self, source, data, ):
        """Converts data times to be in units of days since sourcediscovery

        :param source: Reference Source
        :param data: Data
        :return: Array of converted times
        """
        time_term = self.single_time_pdf(data['timeMJD'], source)
        return time_term

    def S_source(self, source, data, ):
        """Calculates the spatial and time PDFs, and multiplies to return a
        combined PDF value.

        :param source: Source
        :param data: Data
        :return: Spatial * Time PDF value
        """
        space_term = self.space_pdf_signal(source, data, )
        if self.UseTime is True:
            time_term = self.time_pdf_signal(source, data, )
        else:
            time_term = np.ones_like(space_term)

        return space_term * time_term

    def evaluate_sig(self):
        """For the dataset, calculates the Signal PDF. This is stored in
        self.SoB to save memory, but is later replaced by the Signal over
        Background ratio value.
        """
        data = self._ev
        sources = self.sources
        NSources = len(sources)
        NData = len(data)
        self.SoB = np.zeros([NSources, NData])
        for i in range(len(sources)):
            # Stage 1 of calculating SoB, processed further later (saves memory)
            self.SoB[i] = self.S_source(sources[i], data, )

# ==============================================================================
# Time PDFs
# ==============================================================================



    def compute_source_weights_time(self, ):
        """Loops over sources. For any of the given Time Models, calculates
        the overlap between the associated time PDF and the data-taking
        period. Assigns time weights for the source_path, as well as Time
        Normalisation.
        """
        for source in self.sources:
            # Calculates for a box starting at discovery and lasting for
            # a time period TimeBoxLength
            season_norm, total_norm = tm.return_norms(
                self.ReconTimeModel, self.DataStart, self.DataEnd,
                source["discoverydate_mjd"],
                self.ReconTimeParameters)

            source['weight_time'] = max(season_norm / self.SeasonTimeSpan, 0.)
            source['TimeNorm'] = max(season_norm / total_norm, 0.)

            sim_season_norm, sim_total_norm = tm.return_norms(
                self.SimTimeModel, self.DataStart, self.DataEnd,
                source["discoverydate_mjd"],
                self.SimTimeParameters)

            source['sim_TimeNorm'] = max(sim_season_norm / sim_total_norm, 0.)

            # print source_path['weight_time'], source_path['TimeNorm'], season_norm

            # # Calculates for a box beginning TimeBoxLength before discovery
            # # and continuing until discovery.
            # if self.ReconTimeModel == 'BoxPre':
            #     t_start = min(
            #         max(self.DataStart,
            #             source_path['discoverydate_mjd'] - self.ReconTimeLength),
            #         self.DataEnd)
            #
            #     t_end = min(max(self.DataStart, source_path['discoverydate_mjd']),
            #         self.DataEnd)
            #
            #     time_length_in_seasons = t_end - t_start
            #     tot_norm = self.ReconTimeLength
            #     source_path['weight_time'] = max(
            #         time_length_in_seasons / self.SeasonTimeSpan, 0.)
            #     source_path['TimeNorm'] = max(time_length_in_seasons / tot_norm, 0.)

            # Calculates for for an arbitrary exponential decay model
            # if self.ReconTimeModel == 'Decay':
            #     t_pp = self.ReconModelParameters["t_pp"]
            #     t_start = max(
            #         0., self.DataStart - (source_path['discoverydate_mjd']))
            #     t_end = min(self.ReconTimeLength,
            #         max((self.DataEnd - source_path['discoverydate_mjd']), 0.))
            #     tot_norm = t_pp * (
            #         np.log(self.ReconTimeLength + t_pp) - np.log(0. + t_pp))
            #     SeasonNorm = t_pp * (
            #         np.log(t_end + t_pp) - np.log(t_start + t_pp))
            #     source_path['TimeNorm'] = SeasonNorm / tot_norm
            #     source_path['weight_time'] = (
            #         (source_path['TimeNorm'] * tot_norm) / self.SeasonTimeSpan)

    def single_time_pdf(self, t, source):
        """Processes the array of event times (in mjd). Checks if there is
        Time Normalisation, or sets it to 1. Checks if each event occurred in
        the given season. If a time passes, rescales the time (in mjd) to a
        time since source_path discovery.

        For these times, creates an array of zeroes and fills the
        corresponding entries with the normalised neutrino light curve values
        for the rescaled times, or leaves 0 if outside of data window.

        :param t: Array of event times
        :param source: Source
        :return: Array containing expected neutrinos for each event
        """
        t = np.asarray(t)
        r = np.zeros_like(t)
        if source['TimeNorm'] > 0.:
            norm = source['TimeNorm']
        else:
            norm = 1.
        mask = np.logical_and(t > self.DataStart, t < self.DataEnd)
        t = t[mask] - source['discoverydate_mjd']
        r[mask] = self.ReconNuLightCurveFunc(t) / norm
        return r

