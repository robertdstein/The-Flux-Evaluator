import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import time_models as tm


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
        t_min = self.SimTimeParameters["t0"]
        # t_min = 0

        # -------- Simulation Section --------

        # Sets the simulated neutrino light curve

        self.NuLightCurveFunc = lambda t : tm.return_light_curve(t,
            self.SimTimeModel, self.SimTimeParameters)

        # Sets the integrated light curve function, normalised to yield 1 at
        # t = model_end

        IntegralNuLightCurveFunc = lambda t: tm.return_integrated_light_curve(t,
            self.SimTimeModel, self.SimTimeParameters)
        self.IntegralNuLightCurveFunc = np.vectorize(IntegralNuLightCurveFunc)

        # Sets the maximum t covered by the light curve, and generates a
        # range of t values. Linearly interpolates the integral function,
        # as a function of t

        t_max = self.SimTimeParameters["length"] + self.SimTimeParameters["t0"]
        # t_max = self.SimTimeParameters["length"]

        t = np.linspace(t_min, t_max, 1.e4)

        self.InversInterpol = interp1d(
            self.IntegralNuLightCurveFunc(t), t, kind='linear')

        # print t
        # print self.IntegralNuLightCurveFunc(t)
        # print IntegralNuLightCurveFunc(t)
        # raw_input("prompt")

        # -------- Reconstruction Section --------

        # Sets the recon light curve, which is independent of the simulation
        # model

        self.ReconNuLightCurveFunc = lambda t: tm.return_light_curve(t,
            self.ReconTimeModel, self.ReconTimeParameters)

        # Sets the recon integrated light curve function, normalised to yield
        # 1 at t = recon_model_end

        ReconIntegralNuLightCurveFunc = \
            lambda t: tm.return_integrated_light_curve(t,
            self.ReconTimeModel, self.ReconTimeParameters)

        self.ReconIntegralNuLightCurveFunc = np.vectorize(
            ReconIntegralNuLightCurveFunc)

        # Sets the maximum t covered by the recon light curve, and generates a
        # range of t values. Linearly interpolates the recon integral function,
        # as a function of t

        recon_t_min = self.ReconTimeParameters["t0"]

        recon_t_max = self.ReconTimeParameters["length"] +\
                      self.ReconTimeParameters["t0"]

        recon_t = np.linspace(recon_t_min, recon_t_max, 1.e4)

        self.ReconInversInterpol = interp1d(
            self.ReconIntegralNuLightCurveFunc(recon_t), recon_t, kind='linear')

        # raw_input("prompt")

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
