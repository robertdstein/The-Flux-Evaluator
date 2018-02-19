import numpy as np
import copy
import os
import cPickle as pickle
import time
import glob
import ConfigParser
import subprocess

from sys import stdout

import resource

from scipy.optimize import curve_fit
import scipy.optimize
import scipy as scp

from scripts.LLh import LLh


class GenerationControl(object, ):

    """Handles the creation of a Log Likelihood Function for the dataset,
    and the repeated minimisation of this function. Also manages the merging
    of test results.

    """

    def __init__(self, settings=dict()):
        self.settings = settings
        self.seasons = []

        # **********************************************************************
        # Sets tmp directory (currently Alex's)
        try:
            tmpdir = os.environ['TMPDIR']
        except KeyError:
            # tmpdir ="/afs/ifh.de/user/s/steinrob/scratch/PS_Data"
            tmpdir = "/afs/ifh.de/user/a/astasik/scratch/PS_Data/"

        try:
            # Reads in the config variables for each season of data
            data_conf = ConfigParser.ConfigParser()
            data_conf.read(self.settings["DataConfig"])

            # Loops over each season, creating an LLh Object.
            # The data is randomised, and the LLh function is spline-fitted.
            for season in data_conf.sections():
                new_instance = LLh(
                    ExpPath=tmpdir + data_conf.get(season, "exp_path"),
                    MCPath=tmpdir + data_conf.get(season, "mc_path"),
                    AcceptanceWeightPath=tmpdir + data_conf.get(season, "aw_path"),
                    Livetime=eval(data_conf.get(season, "livetime")),
                    StartDataTakingMJD=float(data_conf.get(season, "start_mjd")),
                    EndDataTakingMJD=float(data_conf.get(season, "end_mjd")),
                    **self.settings)
                new_instance.init_everything_for_multiple_trials()
                self.seasons.append(new_instance)

        except KeyError:
            pass

    def generate_test_statistics(self, k=0, n_trials=10, path='test.pkl',
                                 seed=np.nan, **kwargs):
        """Creates the custom path for the given seed and k value
        Generates n_trials and writes results to a pickle file.

        :param k: Flux Scale
        :param n_trials: Number of trials for each scale
        :param path: Path to save pickle file
        :param seed: Random seed
        """
        # **********************************************************************
        # Look into it
        np.random.seed()
        if np.isnan(seed):
            self.seed = np.random.randint(0, 4294967295)
        else:
            self.seed = seed
        np.random.seed(self.seed)

        # Sets the save path for the pickle file, Generates the trials for
        # the given flux scale, and saves the pickle results file
        path += "_" + self.settings["ConfigName"]
        path += '__' + str(k) + '_' + str(self.seed) + '.pkl'

        # Checks the Pickle Directory exists, and if not, creates it
        results_dir = os.path.dirname(path)
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)

        test_stats = self.generate_trials(n_trials, k=k, )
        self.write_result_to_file(path, k, test_stats)

    def memory_usage_ps(self, ):
        """Returns the memory usage."""
        out = subprocess.Popen(
            ['ps', 'v', '-p', str(os.getpid())],
            stdout=subprocess.PIPE).communicate()[0].split(b'\n')

        vsz_index = out[0].split().index(b'RSS')
        mem = float(out[1].split()[vsz_index]) / 1024 / 1.e3
        return mem

    def generate_trials(self, n_trials, k=0., ):
        """Generates trials, by lopping over each season of data and creating a
        combined Log Likelihood Function. Minimises this LLh function,
        returning the result.

        :param n_trials: Number of times the minimisation should be performed
        :param k: Flux scale
        :return: Test statistics result
        """
        test_stats = []

        # Gives the memory usage
        print 0, 'Memory usage: %s (Gb)' % self.memory_usage_ps()
        memory_use = str(float(resource.getrusage(
            resource.RUSAGE_SELF).ru_maxrss) / 1.e6)
        print 'Memory usage max: %s (Gb)' % memory_use

        # Loops over n_trials
        for counter in range(n_trials):
            # Prints the progress
            self.print_progress(counter, n_trials)

            # Loops over seasons to create each season's LLh function
            llh_funcs = []
            for season in self.seasons:
                season.prepare_fake_data_set_and_evaluate_pdf(k, )
                f = season.ProduceLLhFunction()
                llh_funcs.append(f)

            def f_final(x):
                """The likelihood function to be minimised, given the Ice
                Cube Dataset. The included seasons are listed in the
                data_config file, and the name of the data_config file to be
                used is listed in config.ini under the DataConfig field.

                :param x: parameter array for minimisation
                :return: The combined Log Likelihood for all included seasons,
                given parameter set x
                """

                NSignalTotal = 0.

                params = []

                # Loops over each season of data
                for season in self.seasons:
                    # If both "UseEnergy" and "FitGamma" are true,
                    # takes gamma from parameter array.
                    # Otherwise sets Gamma to be the correct value.
                    if np.logical_and(self.settings['UseEnergy'] is True,
                                      self.settings['FitGamma'] is True):
                        gamma = x[-1]
                    else:
                        gamma = season.InjectionGamma
                        p = np.array(x)
                        p = np.append(p, gamma)
                        params.append(p)

                    # Sets the source acceptance
                    for source in season.sources:
                        source['weight_acceptance'] = season.AcceptanceFitFunc(
                            source['dec'], gamma)

                    season.sources['weight_distance'] = (
                        season.sources['distance'] ** (-2.))

                    # Either accounts for  time variability, or
                    # weights seasons equally
                    if self.settings['UseTime']:
                        season.sources['weight_time'] = \
                            season.sources['weight_time']
                    else:
                        season.sources['weight_time'] = np.ones_like(
                            season.sources['weight_time'])

                    # Assigns an overall source weight
                    season.sources['weight'] = (
                        season.sources['weight_distance'] *
                        season.sources['weight_time'] *
                        season.sources['weight_acceptance'])

                    NSignalTotal += np.sum(season.sources['weight'])

                # print x

                # Puts the weights into a matrix
                WeightMatrix = np.zeros(
                    (len(self.seasons), len(self.seasons[0].sources)),
                    dtype=float)
                for i, season in enumerate(self.seasons):
                    WeightMatrix[i] = season.sources['weight']

                # If the weights are not to be fitted, assign a weight of 1
                # to events which lie outside the seasons, and renormalise
                if self.settings['FitWeights'] is False:
                    norm = np.sum(WeightMatrix, axis=1)[:, None]
                    for i, n in enumerate(norm):
                        if n == 0:
                            norm[i] = 1

                    SourceWeights = WeightMatrix / norm
                    SeasonWeights = np.sum(WeightMatrix, axis=1) / np.sum(
                        np.sum(WeightMatrix, axis=0))

                # If weights are to be fitted,
                if self.settings['FitWeights'] is True:
                    SourceWeights = np.ones_like(WeightMatrix)
                    norm = np.sum(WeightMatrix, axis=0)
                    altSeasonWeights = np.zeros_like(WeightMatrix)
                    for i, n in enumerate(norm):
                        if n > 0:
                            altSeasonWeights[:, i] = WeightMatrix[:, i] / n
                    # SeasonWeights = WeightMatrix / np.sum(WeightMatrix, axis=0)
                    SeasonWeights = altSeasonWeights

                # Loops over seasons to assign weights
                for i, season in enumerate(self.seasons):
                    season.SeasonWeight = SeasonWeights[i]
                    season.sources['weight'] = SourceWeights[i]

                # Loops over each LLh function to give overall LLh
                value = 0.0
                # print "LLH values:"
                for j, f in enumerate(llh_funcs):

                    if np.logical_and(self.settings['UseEnergy'] is True,
                                      self.settings['FitGamma'] is False):
                        y = np.array(params[j])
                    else:
                        y = np.array(x)

                    value += f(y)

                return value

            test_stat_res = self.minimise_llh(f_final)

            test_stats.append(test_stat_res)
            del f_final

        # Ensures cluster cache is not full
        for season in self.seasons:
            del season.SoB
        return np.array(test_stats)

    def minimise_llh(self, f_final):
        """Minimises the Log Likelihood Function using scipy optimize,
        taking account of the chosen settings (specified in the config.ini),
        such as whether to fit gamma.

        :param f_final: Combined LLh function
        :return: Value for LLh at minimum
        """
        # Gives number of sources in catalogue
        n_sources = len(np.load(self.settings['SourcePath']))

        # Checks whether to fit weights and gamma, and whether to use energy
        # Then minimises for the chosen configuration
        # Provides appropriate bounds the chosen configuration

        # If fit weights loops over a bound for each source!
        # If not, then gives only one set of bound

        if self.settings['FitWeights'] is True:
            if self.settings['UseEnergy'] is True:

                if self.settings['FitGamma'] is True:
                    seed = np.append(np.ones(n_sources), 2.)
                    bounds = [(0., 1000.) for i in range(n_sources)] + \
                             [(1., 4.)]
                    res = scp.optimize.fmin_l_bfgs_b(
                        f_final, seed, bounds=bounds, approx_grad=True,
                    )

                if self.settings['FitGamma'] is False:
                    seed = np.ones(n_sources)
                    bounds = [(0., 1000.) for i in range(n_sources)]
                    res = scp.optimize.fmin_l_bfgs_b(
                        f_final, seed, bounds=bounds, approx_grad=True)

            if self.settings['UseEnergy'] is False:
                seed = np.zeros(n_sources)
                bounds = [(0., 1000.) for i in range(n_sources)]
                res = scp.optimize.fmin_l_bfgs_b(
                    f_final, seed, bounds=bounds, approx_grad=True)

        if self.settings['FitWeights'] is False:
            if self.settings['UseEnergy'] is True:

                if self.settings['FitGamma'] is True:
                    rranges = (slice(0., 4., 0.5), slice(2., 4., 0.25))
                    seed = scipy.optimize.brute(f_final, rranges, finish=None)
                    res = scp.optimize.fmin_l_bfgs_b(
                        f_final, seed, bounds=[(0., 1000.), (1., 4.)],
                        approx_grad=True)

                if self.settings['FitGamma'] is False:
                    res = scp.optimize.fmin_l_bfgs_b(
                        f_final, [1.], bounds=[(0., 1000.)], approx_grad=True)

            if self.settings['UseEnergy'] is False:
                res = scp.optimize.fmin_l_bfgs_b(
                    f_final, [1.], bounds=[(0., 1000.)], approx_grad=True)

        # Returns a value for LLH at minimum
        test_stat_res = -res[1]

        return test_stat_res

    def write_result_to_file(self, path, k, test_stats):
        """Adds information to a pickle file, and saves it.

        :param path: Path to save pickle file
        :param k: Flux Scale
        :param test_stats: Array of test stat results
        """
        # Checks if Pickle file exists
        self.check_if_pkl_file_exists(path)

        # Opens the pickle file and creates a copy of the data
        pkl_file = open(path, 'rb')
        test_stat_results = pickle.load(pkl_file)
        pkl_file.close()

        # Either creates a new key, or adds results to existing key entry
        if k in test_stat_results.keys():
            test_stat_results[k] = np.append(test_stat_results[k], test_stats)
        else:
            test_stat_results[k] = test_stats

        # Saves the updated pickle file
        pkl_file = open(path, 'wb')
        pickle.dump(test_stat_results, pkl_file)
        pkl_file.close()

    def check_if_pkl_file_exists(self, path, ):
        """Checks if the pickle file exists.
        If not, creates an empty pickle file.

        :param path: pickle file path
        """
        if os.path.isfile(path):
            pass
        else:
            test_stat_results = {}
            with open(path, 'wb') as pkl_file:
                pickle.dump(test_stat_results, pkl_file)
            print path, ' created'

    def print_generation_overview(self, path):
        """Loads the pickle results file.
        Prints the number of entries for each flux.

        :param path: Pickle File path
        """
        self.check_if_pkl_file_exists(path)
        with open(path, 'rb') as pkl_file:
            test_stat_results = pickle.load(pkl_file)
        print ''
        for key in np.sort(test_stat_results.keys()):
            print key, len(test_stat_results[key])
        print ''
        del test_stat_results

    def print_progress(self, counter, n_trials):
        """Prints the progress of function, as a percentage

        :param counter: counter of trial
        :param n_trials: Total umber of trials
        """
        stdout.write("\r%.1f %%" % (float(counter + 1.) / n_trials * 100.))
        stdout.flush()

    def merge_test_result_pickles(self, data_path='test_stat_results/test/test',
                                  output_path="test_stat_results/test",
                                  config_name="Fast_with_fit"):
        """Searches for all pickle files matching the data_path path.
        Merges these into a single pickle file, and saves it to OutPutPath.
        Prints the combined results.

        :param data_path: Path of individual pickle files
        :param output_path: Path to save merged pickle file
        :param config_name: Name of configuration (used for naming pickle file)
        """

        data_path += "_" + config_name
        output_path += "_" + config_name + ".pkl"

        # Checks if results directory exists
        output_dir = os.path.dirname(output_path)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # Creates an empty Pickle Path for results
        if os.path.isfile(output_path):
            os.remove(output_path)

        # Returns a list of all files matching the path given
        file_list = glob.glob(data_path + '__*.pkl')

        test_stat_results = {}

        # Loops over each file and adds it to the combined File
        for singlefile in file_list:
            pkl_file = open(singlefile, 'rb')
            single_result = pickle.load(pkl_file)
            pkl_file.close()
            for k in single_result.keys():
                if k in test_stat_results.keys():
                    test_stat_results[k] = np.append(test_stat_results[k],
                                                     single_result[k])
                else:
                    test_stat_results[k] = single_result[k]

        # Saves merged file and prints results
        with open(output_path, 'wb') as pkl_file:
            pickle.dump(test_stat_results, pkl_file)
        self.print_generation_overview(output_path)
