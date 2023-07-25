import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
matplotlib.use("agg")
import argparse
import subprocess
import os
from joblib import Parallel, delayed
import scipy.interpolate
import h5py

import matplotlib
matplotlib.use('Agg')
import os
import sys
import sphviewer
import pdb
from matplotlib.pyplot import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sphviewer.tools import QuickView
from glob import glob
import trident

import yt
import trident
from numpy import *
from yt import *
yt.enable_parallelism()
from mpl_toolkits.axes_grid1 import AxesGrid
#import glob
from yt import YTQuantity
from matplotlib.pyplot import *
from matplotlib.pyplot import cm


parser = argparse.ArgumentParser()
parser.add_argument('--param1', type=int)
parser.add_argument('--param2', type=str)
#parser.add_argument('--param3', type=str)

args = parser.parse_args()

jj = args.param1     ## 5 or 200  (kpc)
kk = args.param2     ## simulation-type
#ll = args.param3     ## different ions

cuts = np.array([5, 200])

simulations_name = np.array(['jet1', 'jet2', ...])

# ....    parametrically


NN = cuts[jj]
simulatoin = simulations_name[kk]


def Loop_snap(snapshot):
    ## Includes the yt, trident, the main body of your code ...
    ##ist of the ions l
    plt.savefig( '/figure' + str(snap_Num).zfill(4) + '.png', dpi=150)
                    

parallel(n_jobs=48)(delayed(Loop_snap)(mm) for mm in range(len(number_snapshot)))
