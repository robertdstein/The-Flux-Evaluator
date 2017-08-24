import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# **********************************************************************************
# Is this even used????
class RandomTools(object, ):

    """Contains functions for randomness. called Init

    """

    def __init__(self, ):
        pass

    def reset_random_generator_pdf(self, ):
        """Resets the attributes InversInterpol and RandFuncNorm to their
        default values of nan.
        """
        self.InversInterpol = np.nan
        self.RandomFuncNorm = np.nan

    def init_random_generator_pdf(self, ):
        """Produces neutrino light curve for given time model.
        Also creates integrated light curve. Creates an inverse interpolation
        of the integrated light curve, so that a fraction between 0 and 1
        can be converted into a time.

        :return: Inverse interpolation of integrated light curve
        """
        t_min = 0.
        if self.TimeModel == 'Box':
            self.NuLightCurveFunc = lambda x: self.BoxFunc(x,
                                                           self.TimeBoxLength)
            t_max = self.TimeBoxLength
            norm = (t_max - t_min)
            IntegralNuLightCurveFunc = lambda x: x / norm
        if self.TimeModel == 'BoxPre':
            self.NuLightCurveFunc = lambda x: self.BoxFunc(x,
                                                           self.TimeBoxLength)
            t_max = self.TimeBoxLength
            norm = (t_max - t_min)
            IntegralNuLightCurveFunc = lambda x: x / norm
        if self.TimeModel == 'Decay':
            self.NuLightCurveFunc = lambda x: self.AnalyticTimePDF(
                x, self.Model_tpp, self.DecayModelLength)
            t_max = self.DecayModelLength
            norm = self.Model_tpp * (np.log(t_max + self.Model_tpp) -
                                     np.log(t_min + self.Model_tpp))
            IntegralNuLightCurveFunc = lambda x: self.Model_tpp * (
                np.log(x + self.Model_tpp) - np.log(self.Model_tpp)) / norm

        self.IntegralNuLightCurveFunc = np.vectorize(IntegralNuLightCurveFunc)
        t = np.linspace(t_min, t_max, 1.e4)
        y = self.IntegralNuLightCurveFunc(t)

        self.InversInterpol = interp1d(self.IntegralNuLightCurveFunc(t), t,
                                       kind='linear')
        return self.InversInterpol

    def generate_n_random_numbers(self, n_events, ):
        """Generates random numbers, and uses InversInterpol to convert these
        values back to a light curve.

        :param n_events: Number of points to convert
        :return: A set of light curve time values
        """
        values = self.InversInterpol(np.random.uniform(0., 1., n_events))
        return values

    # Not Called!
    def plot_test_and_pdf(self, n_events, t_min=0., t_max=1.,
                          path='TimePDFTest.pdf'):
        """Produces a graph to test the behaviour of the time models.
        Randomly chooses values with GeerateNRandomNumbers, and creates a
        histogram of the data. PLots the original PDF for comparison

        :param n_events: Number of events to simulate
        :param t_min: Minimum time value for curve
        :param t_max: Maximum time value for curve
        :param path: Path to save file
        """
        extra_plot_space = 100.
        values = self.generate_n_random_numbers(n_events)
        x = np.linspace(
            t_min - extra_plot_space, t_max + extra_plot_space, 1000 + 1)

        plt.figure()
        plt.hist(values, bins=50+1, histtype='stepfilled', lw=0, normed=True,
                 label='data', alpha=0.5)
        plt.plot(x, self.NuLightCurveFunc(x), lw=2, color='blue',
                 label=r'$\nu$ light curve')

        plt.grid()
        plt.xlim(t_min - extra_plot_space, t_max + extra_plot_space)
        plt.legend(loc='best', fancybox=True, framealpha=1.)
        plt.xlabel('Time [d]')
        plt.ylabel('PDF value [a.u.]')
        plt.savefig(path)
        plt.show()
