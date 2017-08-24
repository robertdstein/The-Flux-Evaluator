import numpy as np
import copy
import Sphere
import RandomTools
from scipy import stats
from scipy.stats import norm
import healpy as hp

class Injector():

    def __init__(self, ):
        pass

    def generate_sig_events(self, sources, mc, k=0., ):
        """Produces signal events based on the Monte Carlo signal data.
        Smears/Randomises the Monte Carlo dataset in various ways.

        :param sources: list of sources
        :param mc: Monte Carlo data (signal)
        :param k: Flux scaling
        :return: Signal events
        """
        inject_sources = copy.deepcopy(sources)

        # Random Variates drawn from a normal distribution, with a width
        # scale equal to 'SmearInjection'.
        new_flux = inject_sources['flux'] * np.exp(stats.norm.rvs(
            loc=0., scale=self.SmearInjection, size=len(sources)))
        normalisation = np.sum(sources['flux']) / np.sum(new_flux)
        inject_sources['flux'] = new_flux * normalisation * k

        sig_events = self.extract_events_from_mc(inject_sources, mc,)
        return sig_events

    def find_and_apply_band_mask(self, source, mc, dec_bandwidth):
        """For a given source, creates a mask to include only Monte Carlo
        events which lie in a declination band, of width 10 degrees, centered
        on the source.

        :param source: Single source
        :param mc: Monte Carlo data (signal)
        :param dec_bandwidth: Width of declination band (in radians)
        :return: Returns Monte Carlo events in the band,
        the solid angle coverage, and the mask for the band.
        """
        # Sets a declination band 5 degrees above and below each source
        min_dec = max(-np.pi / 2., source['dec'] - dec_bandwidth)
        max_dec = min(np.pi / 2., source['dec'] + dec_bandwidth)
        # Gives the solid angle coverage of the sky for the band
        omega = 2. * np.pi * (np.sin(max_dec) - np.sin(min_dec))

        # Check for season weighting as a function of declination
        if self._ReturnInjectorNExp:
            self.InjectionBandMask = dict()

        # Checks if the source mask has already been evaluated.
        # If not, creates a mask for MC events lying in the band.
        if source['name'] in self.InjectionBandMask.keys():
            band_mask = self.InjectionBandMask[source['name']]
        else:
            band_mask = np.logical_and(np.greater(mc["trueDec"], min_dec),
                                       np.less(mc["trueDec"], max_dec))
            self.InjectionBandMask[source['name']] = band_mask
        return (mc[band_mask]), omega, band_mask

    def extract_events_from_mc(self, sources, mc, ):
        """Randomises a a set of Monte Carlo data for a given set of sources.
        For each source, selects only those events lying in a declination
        band centered on the source.

        :param sources: list of sources
        :param mc: Monte Carlo (signal)
        :return: signal events
        """
        sig_events = np.empty((0, ),
                              dtype=[("ra", np.float), ("sinDec", np.float),
                                     ("sigma", np.float), ("logE", np.float),
                                     ("dec", np.float), ('timeMJD', np.float),
                                     ])

        # Sets detector livetime in seconds
        livetime = self.Livetime * (60. * 60. * 24.)
        #**************************************************************************
        # Sometimes 10 degrees, sometimes 5? (PDF selects_events_in_band)
        # Conservative
        dec_bandwidth = np.deg2rad(5.)
        TotMuN = 0.

        # Loops over sources to add in expected number of neutrinos
        for source in sources:
            # Only includes events lying in a +/- 5 degree declination band
            SourceMC, omega, band_mask = self.find_and_apply_band_mask(
                source, mc, dec_bandwidth)

            # If using time, calculates the Fluence using the time model
            # Otherwise calculates Fluence as flux * livetime (i.e. constant)
            if self.UseTime is True:
                EfficencyFactor = livetime / (
                    (self.DataEnd-self.DataStart) * 24. * 60. * 60.)
                if self.TimeModel == 'Box':
                    TotalTime = self.TimeBoxLenght
                if self.TimeModel == 'BoxPre':
                    TotalTime = self.TimeBoxLenght
                if self.TimeModel == 'Decay':
                    TotalTime = self.Model_tpp * (
                        np.log(self.Model_Length + self.Model_tpp)
                        - np.log(0. + self.Model_tpp))
                TotalTimeDays = TotalTime
                TotalTime *= (60. * 60. * 24.)
                fluence = EfficencyFactor * TotalTime * source['flux']
            else:
                fluence = source['flux'] * livetime

            if self._ReturnInjectorNExp is True:
                source['weight_distance'] = 1.

            # **************************************************************************
            # Add in bracket?
            # Recalculates the one weights to account for the band mask?
            SourceMC['ow'] = (self.WeightsInject[band_mask] / omega) * source[
                'weight_distance'] * fluence

            # Expectation number of Neutrinos, equal to sum of one weights
            MuN = np.sum(SourceMC['ow'], dtype=np.float)

            if self._ReturnInjectorNExp is True:
                return MuN

            # If weighting for time, calculates weighted expectation value
            # for number of neutrinos, and adds it to total expectation value.
            if self.UseTime is True:
                TotMuN += (MuN * source['weight_time'] / (TotalTimeDays /
                                                          self.SeasonTimeSpan))
            # If not time weighting, simply adds expectation value to total.
            else:
                TotMuN += MuN

            # Draws random number from expectation value of neutrinos
            n_signal = np.random.poisson(MuN)
            if n_signal < 1:
                continue

            # Creates a normalised array of OneWeights
            p_select = SourceMC['ow'] / np.sum(SourceMC['ow'])

            # Creates an array with n_signal entries.
            # Each entry is a random integer between 0 and no. of sources.
            # The probability for each integer is equal to the OneWeight of
            # the corresponding source.
            ind = np.random.choice(len(SourceMC['ow']), size=n_signal,
                                   p=p_select)

            # Selects the sources corresponding to the random integer array
            sam_ev = SourceMC[ind]

            # Rotates the Monte Carlo events onto the source
            sam_ev = self.rotate_struct(sam_ev, source['ra'], source['dec'])

            # Generates random numbers according to Time profile
            if self.UseTime is True:
                sam_ev['timeMJD'] = (self.GenerateNRandomNumbers(n_signal) +
                                     source['discoverydate_mjd'])
                sam_ev = self.check_time_borders(sam_ev, )

            sig_events = np.concatenate((sig_events, sam_ev))
        return sig_events

    def check_time_borders(self, sam_ev, ):
        """Checks to ensure each event lies within the season.
        Removes those events which lie outside the window.

        :param sam_ev: events
        :return: events (modified)
        """
        mask = np.logical_and(sam_ev['timeMJD'] > self.DataStart,
                              sam_ev['timeMJD'] < self.DataEnd)
        return sam_ev[mask]

    def merge_struct_arrays(self, data1, data2):
        """Joins two datasets, using np.concatenate().

        :param data1: first dataset
        :param data2: second dataset
        :return: joined dataset
        """
        data_final = np.concatenate((data1, data2))
        return data_final

    def scramble_exp_data(self, exp):
        """Scrambles "Experimental" data set.
        If the data is unblinded, returns the dataset and a warning.
        If not unblinded, then produces a new dataset,
        drawing values randomly from the original dataset.

        :param exp: Experimental Data
        :return: Experimental Data (shuffled)
        """
        if self.Unblind is False:
            if self.BootStrap is True:
                N_tot = len(exp)
                # Produces N_tot random integers, each between 0 and N_tot
                x = np.random.choice(N_tot, N_tot)
                # Reassigns each entry of exp,
                # using a randomly selected value from the dataset
                exp = exp[x]
            # Assigns a flat random distribution for Right Ascension
            exp['ra'] = np.random.uniform(0, 2 * np.pi, size=len(exp))
            # Randomly reorders the times
            np.random.shuffle(exp['timeMJD'])
        else:
            print('WARNING: Running in Unblinding Mode')
        return exp

    def inject_signal_events(self, exp, mc, sources, k=0., ):
        """Adds signal events to to a dataset containing experimental data
        (background).

        :param exp: Experimental Data (background)
        :param mc: Monte Carlo Data (signal)
        :param sources: Set of sources
        :param k: Flux scale
        :return: Full dataset (original with added signal), and the new
        signal events set
        """
        sig_events = self.generate_sig_events(sources, mc, k)
        exp = self.scramble_exp_data(exp)
        data = self.merge_struct_arrays(exp, sig_events)
        return data, sig_events

# ****************************************************************************************************
    def rotate(self, ra1, dec1, ra2, dec2, ra3, dec3):
        """Rotate ra1 and dec1 in a way that ra2 and dec2 will exactly map
        onto ra3 and dec3, respectively. All angles are treated as radians.
        Essentially rotates the events, so that they behave as if they were
        originally incident on the source.

        :param ra1: Event Right Ascension
        :param dec1: Event Declination
        :param ra2: True Event Right Ascension
        :param dec2: True Event Declination
        :param ra3: Source Right Ascension
        :param dec3: Source Declination
        :return: Returns new Right Ascensions and Declinations
        """
        # Turns Right Ascension/Declination into Azimuth/Zenith for healpy
        phi1 = ra1 - np.pi
        zen1 = np.pi/2. - dec1
        phi2 = ra2 - np.pi
        zen2 = np.pi/2. - dec2
        phi3 = ra3 - np.pi
        zen3 = np.pi/2. - dec3

        # Rotate each ra1 and dec1 towards the pole?
        x = np.array([hp.rotator.rotateDirection(
            hp.rotator.get_rotation_matrix((dp, -dz, 0.))[0], z, p)
            for z, p, dz, dp in zip(zen1, phi1, zen2, phi2)])

        # Rotate **all** these vectors towards ra3, dec3 (source)
        zen, phi = hp.rotator.rotateDirection(np.dot(
            hp.rotator.get_rotation_matrix((-phi3, 0, 0))[0],
            hp.rotator.get_rotation_matrix((0, zen3, 0.))[0]), x[:, 0], x[:, 1])

        dec = np.pi/2. - zen
        ra = phi + np.pi
        return np.atleast_1d(ra), np.atleast_1d(dec)

    # **************************************************************************************
    def rotate_struct(self, ev, ra, dec):
        """Modifies the events by reassigning the ra of the events
        Also does something to declination?
        Removes the True Monte Carlo information from sampled events,
        (ra/dec/Energy/OneWeights)

        :param ev: events
        :param ra: Right Ascensions
        :param dec: Declinations
        :return: events (modified)
        """
        names = ev.dtype.names
        ev["ra"], rot_dec = self.rotate(ev["ra"], np.arcsin(ev["sinDec"]),
                                        ev["trueRa"], ev["trueDec"],
                                        ra, dec)

        # ???????? Why reassign sindec but not dec?
        if "dec" in names:
            ev["dec"] = rot_dec
        ev["sinDec"] = np.sin(rot_dec)
        # "delete" Monte Carlo information from sampled events
        non_mc = [name for name in names
                  if name not in ["trueRa", "trueDec", "trueE", "ow"]]
        ev = ev[non_mc].copy()

        return ev
