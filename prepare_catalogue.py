"""Script to produce catalogues for use in stacking analysis.

The catalogues themselves are randomly produced for the purpose of trialing
the code. Modification of variable n can produces a catalogue with an
arbitrary number of sources.

"""

import numpy as np

root = "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/Catalogue/"


def read_in_catalogue():
    """Produces a catalogue with a single source.

    :return: Source Array
    """
    sources = np.empty(
        1, dtype=[("ra", np.float), ("dec", np.float),
                  ("flux", np.float), ("n_exp", np.float),
                  ("weight", np.float), ("weight_acceptance", np.float),
                  ("weight_time", np.float),
                  ("weight_distance", np.float),
                  ("discoverydate_mjd", np.float),
                  ("distance", np.float), ('name', 'a30'),
                  ])

    sources['ra'] = np.array([np.deg2rad(180.)])
    sources['dec'] = np.arcsin(0.99)
    sources['flux'] = np.array([1.e-9])
    sources['weight'] = np.array([1.0])
    sources['distance'] = np.array([1.0])
    sources['discoverydate_mjd'] = (
        np.array([55800.4164699]) )
    sources['name'] = 'SN_01'

    return sources

# Saves a single-source catalogue
single_source_array = read_in_catalogue()
np.save(root + "catalogue00.npy", single_source_array)


def read_in_catalogue_stack(n_sources):
    """Produces a catalogue of n sources. Attributes are randomised within
    physical bounds.

    :param n_sources: Number of sources in catalogue
    :return: Source Array
    """

    sources = np.empty(
        n, dtype=[("ra", np.float), ("dec", np.float),
                  ("flux", np.float), ("n_exp", np.float),
                  ("weight", np.float), ("weight_acceptance", np.float),
                  ("weight_time", np.float),
                  ("weight_distance", np.float),
                  ("norm_time", np.float),
                  ("global_weight_norm_time", np.float),
                  ("discoverydate_mjd", np.float),
                  ("distance", np.float), ('name', 'a30'),
                  ])

    sources['ra'] = np.deg2rad(np.linspace(0., 360, n_sources + 1)[:-1])
    sources['dec'] = np.deg2rad(np.linspace(-90., 90., n_sources + 2)[1:-1])

    sources['flux'] = np.ones_like(sources['ra'])
    normalisation = n_sources * 1.e-9 / np.sum(sources['flux'])
    sources['flux'] *= normalisation
    sources['weight'] = np.ones_like(sources['ra'])
    sources['distance'] = np.ones_like(sources['ra'])
    sources['discoverydate_mjd'] = (
        55694.4164699 + (np.array(range(n_sources)) / float(n_sources)) *
        368.00423609999416)
    sources['name'] = ['SN' + str(i) for i in range(n_sources)]

    return sources

# Saves an n-source catalogue
n = 10
n_source_array = read_in_catalogue_stack(n)
np.save(root + "catalogue_stack" + str(n) + ".npy", n_source_array)
