import numpy as np
import scipy.optimize
import scipy.stats


class Chi2_LeftTruncated(object):
    """ A class similar to the ones from scipy.stats
       allowing to fit left-truncated chi^2 distributions.
    """
   
    def __init__(self, data, cut=0., **kwargs):
        """ Fit the given ensemble of measurements with a chi^2 function.

        `data` is a list of test statistics values.
        `cut` defines where the distribution is truncated.
        """

        data_left  = data[data <= cut]
        data_right = data[data > cut]
        N_all      = len(data)
        N_left     = len(data_left)

        # three parameters will be fitted: dof, location, scale
        p_start  = [ 2., -1., 1. ]
        p_bounds = [ (0., None),    #      dof > 0
                    (None, 0.),    # location < 0 for 'truncated' effect
                    (1e-5, 1e5) ]  #    shape ~ free

        # define the fit function: likelihood for chi^2 distribution,
        # plus knowledge about the amount of truncated data
        def func(p):
           dist = scipy.stats.chi2(*p)
           loglh  = dist.logpdf(data_right).sum()
        #           loglh += N_left*dist.logcdf(cut)
           return -loglh

        res = scipy.optimize.minimize(func,x0=[2., 1., 1.])
        #                                     x0=p_start,)
        #                                     bounds=p_bounds)

        if not res.success:
           print 'Chi2 fit did not converge! Result is likely garbage.'

        self._q_left = N_left / float(N_all)
        self._cut = cut
        self._f   = scipy.stats.chi2(*res.x)
        self._ks  = scipy.stats.kstest(data_right, self._f.cdf)[0]


    def pdf(self, x):
        """ Probability density function.
        """
        m_flat = (0.<=x) & (x<=self._cut)
        m_chi2 = (x > self._cut)

        x = np.asarray(x)
        r = np.zeros_like(x)
        if self._cut == 0.:
           r[m_flat] = self._f.cdf(self._cut)
        else:
           r[m_flat] = self._f.cdf(self._cut) / self._cut
        r[m_chi2] = self._f.pdf(x[m_chi2])
        return r


    def cdf(self, x):
        """ Cumulative distribution function.
        """
        return self._f.cdf(x)


    def sf(self, x):
        """ Survival function.
        """
        return self._f.sf(x)


    def isf(self, x):
        """ Inverse survival function.
        """
        return self._f.isf(x)


    def ppf(self, x):
        """ Percent point function.
        """
        return self._f.ppf(x)


    def __str__(self):
        return 'Left-truncated Chi^2 Distribution:\n'+\
              '\t DoF      = {0:7.2f}\n'.format(self._f.args[0])+\
              '\t Location = {0:7.2f}\n'.format(self._f.args[1])+\
              '\t Scale    = {0:7.2f}\n'.format(self._f.args[2])+\
              '\t KS       = {0:7.2%}'.format(self._ks)




def fraction_above_background(ts, n_inj, n_eff, ts_bkg, error=False, ordered=False):
   """
   Calculate the fraction of trials above a given threshold (ts_bkg).
   The test statistic (ts) and number of (true) injected signal events
   for each trial is reweighted to an effective number of signal events
   (n_eff) in the Poisson process.

   If the fits are ordered by ascending TS, set `ordered` to True.

   The error calculation originates from SkyLab.
   """

   w = poisson_weight(n_eff, n_inj)
   if ordered:
       i_above = np.searchsorted(ts, ts_bkg, side='right')
       sum_above = w[i_above:].sum()
       sum_below = w[:i_above].sum()

       if error:
           err = np.sqrt(np.sum(w[i_above:]**2)) / (sum_below + sum_above)
           return sum_above / (sum_below + sum_above), err
       else:
           return sum_above / (sum_below + sum_above)
       # NOTE: The following should be even faster, provided that sum(w) ~= 1.
       # return w[i_above:].sum()
   else:
       above = ts > ts_bkg
       sum_above = w[above].sum()
       sum_below = w[~above].sum()

       if error:
           err = np.sqrt(np.sum(w[above]**2.)) / (sum_below + sum_above)
           return sum_above / (sum_below + sum_above), err
       else:
           return sum_above / (sum_below + sum_above)
