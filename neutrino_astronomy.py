import cPickle as pickle

from common import tde_pickle

with open(tde_pickle) as f:
    tde_dict = pickle.load(f)

for (tde, vals) in tde_dict.items():
    print tde, vals["polynom"]
