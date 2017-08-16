
# coding: utf-8

# In[16]:

import numpy as np

root = "/afs/ifh.de/user/s/steinrob/scratch/PS_Data/Catalogue/"

def read_in_catalogue():
	sources = np.empty((1), dtype=[("ra", np.float), ("dec", np.float),
									("flux", np.float), ("n_exp", np.float),
									("weight", np.float),
									("weight_acceptance", np.float),
									("weight_time", np.float),
									("weight_distance", np.float),
									("norm_time", np.float),
									("discoverydate_mjd", np.float),
									("distance", np.float),
									('name', 'a30'),
								])
	sources['ra'] = np.array([np.deg2rad(180.)])
	sources['dec'] = np.arcsin(-0.4)
	sources['flux'] = np.array([1.e-9])
	sources['weight'] = np.array([1.0])
	sources['distance'] = np.array([1.0])
	sources['discoverydate_mjd'] = np.array([55694.4164699])# + 1.*368.00423609999416
	sources['name'] = 'SN_01'

	return sources

sources = read_in_catalogue()
np.save(root + 'catalogue00.npy', sources)

def read_in_catalogue():
	sources = np.empty((1), dtype=[("ra", np.float), ("dec", np.float),
									("flux", np.float), ("n_exp", np.float),
									("weight", np.float),
									("weight_acceptance", np.float),
									("weight_time", np.float),
									("weight_distance", np.float),
									("discoverydate_mjd", np.float),
									("distance", np.float),
									('name', 'a30'),
								])
	sources['ra'] = np.array([np.deg2rad(180.)])
	sources['dec'] = np.arcsin(0.99)
	sources['flux'] = np.array([1.e-9])
	sources['weight'] = np.array([1.0])
	sources['distance'] = np.array([1.0])
	sources['discoverydate_mjd'] = np.array([55694.4164699]) - 0.5*368.00423609999416
	sources['name'] = 'SN_01'

	return sources

sources = read_in_catalogue()
np.save(root+"catalogue00.npy", sources)


def read_in_catalogue_stack(n):
	n_sources = n

	sources = np.empty((n_sources), dtype=[("ra", np.float), ("dec", np.float),
									("flux", np.float), ("n_exp", np.float),
									("weight", np.float),
									("weight_acceptance", np.float),
									("weight_time", np.float),
									("weight_distance", np.float),
									("norm_time", np.float),
									("global_weight_norm_time", np.float),
									("discoverydate_mjd", np.float),
									("distance", np.float),
									('name', 'a30'),
									])
	sources['ra'] = np.deg2rad(np.linspace(0., 360, n_sources + 1)[:-1])
	sources['dec'] = np.deg2rad(np.linspace(-90., 90., n_sources + 2)[1:-1])

#    sources['dec'] = 0. * np.ones_like(sources['ra'])

	sources['flux'] = np.ones_like(sources['ra'])
	Norm = n_sources * 1.e-9 / np.sum(sources['flux'])
	sources['flux'] = sources['flux'] * Norm
	sources['weight'] = np.ones_like(sources['ra'])
	sources['distance'] = np.ones_like(sources['ra'])
	sources['discoverydate_mjd'] = np.ones_like(sources['ra']) * np.array([55694.4164699])

	sources['discoverydate_mjd'] = 55694.4164699 + (np.array(range(5)) / float(n)) * 368.00423609999416
	print sources['discoverydate_mjd']
	sources['name'] = ['SN' + str(i) for i in range(n_sources)]

	print Norm

	return sources

n = 5
sources = read_in_catalogue_stack(n)
np.save(root + "catalogue_stack'+str(n)+'.npy", sources)
np.save(root + "/catalogue_stack5.npy", sources)
