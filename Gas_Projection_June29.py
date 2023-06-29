import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
#import illustris_python as il
import h5py
#import illustris_python.snapshot as snapp
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

unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}


"""
Here we present some explanations.

Camera:
    The Camera class is a container that stores the camera parameters. The camera is an object that lives in the space and has spherical coordinates (r,theta,phi), centred around the location (x,y,z). Angles *theta* and *phi* are given in degrees, and enable to rotate the camera along the x-axis and y-axis, respectively. The *roll* angle induces rotations along the line-of-sight, i.e., the z-axis.
"""


### Interesting toturial: http://alejandrobll.github.io/py-sphviewer/content/tutorial_projections.html

### Here we study the density profile of the gas particles. We consider two different cases. 1) When we take into account all of the gas particles and monitor their density profile. 2) When we consider the star forming gas and see how do they look like!!!



#data.keys()
#data["PartType0"].keys()


"""
<KeysViewHDF5 ['CoolingRate', 'Coordinates', 'Density', 'DivBcleaningFunctionGradPhi', 'DivBcleaningFunctionPhi', 'DivergenceOfMagneticField', 'ElectronAbundance', 'HeatingRate', 'HydroHeatingRate', 'IMFFormationProperties', 'InternalEnergy', 'MagneticField', 'Masses', 'MetalCoolingRate', 'Metallicity', 'NetHeatingRateQ', 'NeutralHydrogenAbundance', 'ParticleChildIDsNumber', 'ParticleIDGenerationNumber', 'ParticleIDs', 'SmoothingLength', 'StarFormationRate', 'Velocities']>

"""
# you could put a for loop here that replaces 'm12_mcvt_m2_10000_tor4_pr45_100Myr' with whaterver directories are in /hernquist_lab/AGN_Feedback_Fire then goes to the 'output' folder in each 
Directory_path = '/n/holylfs05/LABS/hernquist_lab/AGN_Feedback_Fire/m12_mcvt_m2_10000_tor4_pr45_100Myr_lower/output/'

diir = glob(Directory_path + '*.hdf5')

List_diir = np.atleast_1d(diir)
snap_order = []
List_diir_image = []

# this calls all of the snapshots in the output folder
for ii in range(len(List_diir)): 
    AA = List_diir[ii]
    snap_ii = int(AA[AA.find('snapshot_')+len('snapshot_'):AA.find('.hdf5')])
    snap_order = np.append(snap_order, [snap_ii])
    List_diir_image = np.append(List_diir_image, List_diir[ii])
    # add txt conversion line later! 

# this puts the snapshots in the correct numerical order when they are runa nd saved 
int_snap = snap_order.astype(int) 
sorted = np.sort(int_snap)
List_diir = List_diir_image[np.argsort(int_snap)]

@derived_field(name= ('PartType0', 'number_density'), units="cm**(-3)", sampling_type="cell",force_override=True)
def _dinos(field, data):
    x_H=0.76
    proton_mass=YTArray([1.67262178*10**-24],'g')
    particlemass   = 4.0 / (3.0*x_H + 1.0 + 4.0*x_H*np.array(data[('PartType0','ElectronAbundance')]))*proton_mass
    den= data['PartType0','Density'] /particlemass
    return den

oo = 0
for fname in List_diir:
    data = h5py.File(fname,"r")
    BH_Center = data["PartType5"]["Coordinates"][:]
    print(BH_Center)
    BH_Center = data["PartType5"]["Coordinates"][:]
    Gas_location = data["PartType0"]["Coordinates"][:] - BH_Center
    Gas_mass = 1e10*data["PartType0"]["Masses"][:]
    Gas_Softening = data["PartType0"]["SmoothingLength"][:]
    
    # this is where I started adding additional parameters
    ds = yt.load(fname, unit_base=unit_base) 
    #ds = yt.load(Directory_path)
    ions_names=['H I','Mg II','C IV','N V','O VI', 'O VII', 'O VIII', 'Ne VIII']  ## Exploration!
    trident.add_ion_fields(ds, ions=ions_names, ftype="PartType0")
    # H1_Number = ds.all_data()[('PartType0', 'H_p1_number_density')]
    # O6_ion_Fraction = ds.all_data()[('gas', 'O_p6_ion_fraction')]
    # O6_ion_mass = ds.all_data()[('gas', 'O_p6_ion_mass')] # MIGHT NEED TO CHANGE THIS TO SAY MASS INSTEAD OF DENSITY
    # Oxygen6 = O6_ion_Fraction
    Oxygen6_mass = ds.all_data()[('gas', 'O_p6_ion_fraction')] # THIS LINE ISNT DONE!!!!
    
    NN = 100 # this defines the x, y, and z axis ranges for the plots 
    #mass = np.ones(len(gas_Loc[:,0]))
    hh = Gas_Softening
    mass = Gas_mass
    Particles = sphviewer.Particles(Gas_location, Oxygen6_mass, hh)

    Scene = sphviewer.Scene(Particles)
    
    fig = plt.figure(1,figsize=(15,5))
    
    #set a figure title on top
    fig.suptitle(r"Precessing kinetic jet with low energy flux: " + str(sorted[oo]+1) + " with Gas_Softening, Gas_mass, and O6_ion_Fraction", fontsize=17, x=0.5, y=1.4)
    
    plt.subplots_adjust(top =1.8, bottom=0.2,hspace=0.3, wspace=0.3)
    
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    extendd = [-NN,NN,-NN,NN]
    Scene.update_camera(r='infinity', t=0, p = 0, roll = 0, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent= extendd)

    Render = sphviewer.Render(Scene)
    Render.set_logscale()
    img1 = Render.get_image()
    extent1 = Render.get_extent()
    divider = make_axes_locatable(ax1)

    #ax1.imshow(img1, extent=extent1, origin='lower', cmap=plt.cm.jet, vmax= 3.0, rasterized=True)
    image1 = ax1.imshow(img1, extent=extent1, origin='lower', cmap=plt.cm.jet, rasterized=True)
    divider = make_axes_locatable(ax1)
    cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    fig.add_axes(cax)
    cb = fig.colorbar(image1, cax=cax, orientation="horizontal")
    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    cb.ax.tick_params(labelsize=15)

    #cb = fig.colorbar(img1, cax=cax, orientation="horizontal")
    #cb = fig.colorbar(img1)

    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    #cb.ax.tick_params(labelsize=15)
    #ax1.set_title(labb_tit[uu])
    ax1.set_xlabel('$X$(kpc)', size=12)
    ax1.set_ylabel('$Y$(kpc)', size=12)

    Scene.update_camera(r='infinity', t=-90, p = -90, roll = 0, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent= extendd)
    
    Render = sphviewer.Render(Scene)
    Render.set_logscale()
    img2 = Render.get_image()
    extent2 = Render.get_extent()
    #divider = make_axes_locatable(ax2)
    #ax2.imshow(img2, extent=extent2, origin='lower',cmap=plt.cm.jet, vmax= 3.0, rasterized=True)
    #image2 = ax2.imshow(img2/imag3, extent=extent2, origin='lower',cmap=plt.cm.jet, rasterized=True)
    image2 = ax2.imshow(img2, extent=extent2, origin='lower', cmap=plt.cm.jet, rasterized=True)

    divider = make_axes_locatable(ax2)
    cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    fig.add_axes(cax)
    cb = fig.colorbar(image2, cax=cax, orientation="horizontal")
    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    cb.ax.tick_params(labelsize=15)


    #cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    #cb = fig.colorbar(img2)

    #fig.add_axes(cax)

    #cb = fig.colorbar(img2, cax=cax, orientation="horizontal")
    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    #cb.ax.tick_params(labelsize=15)
    #ax1.set_title(labb_tit[uu])
    #ax2.set_title(labb_tit[uu])
    ax2.set_xlabel('$Y$(kpc)', size=12)
    ax2.set_ylabel('$Z$(kpc)', size=12)

    Scene.update_camera(r='infinity', t=90, p = 0, roll = -90, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent= extendd)

    Render = sphviewer.Render(Scene)
    Render.set_logscale()
    img3 = Render.get_image()
    extent3 = Render.get_extent()
    divider = make_axes_locatable(ax3)
    #ax3.imshow(img3, extent=extent3, origin='lower', cmap=plt.cm.jet, vmax= 3.0, rasterized=True)
    image3 = ax3.imshow(img3, extent=extent2, origin='lower',cmap=plt.cm.jet, rasterized=True)

    cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    fig.add_axes(cax)
    cb = fig.colorbar(image3, cax=cax, orientation="horizontal")
    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    cb.ax.tick_params(labelsize=15)


    #cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    #fig.add_axes(cax)
    #cb = fig.colorbar(img2)

    #cb = fig.colorbar(img3, cax=cax, orientation="horizontal")
    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    #cb.ax.tick_params(labelsize=15)
    #ax1.set_title(labb_tit[uu])
    #ax3.set_title(labb_tit[uu])
    ax3.set_xlabel('$Z$(kpc)', size=12)
    ax3.set_ylabel('$X$(kpc)', size=12)
    
    Default_dir = '/n/home13/nqutob/AGN_Feedback/ion_snapshots/' + 'm12_mcvt_m2_10000_tor4_pr45_100Myr_lower'
    Default_dir_pdf = '/n/home13/nqutob/AGN_Feedback/ion_snapshots/' + 'm12_mcvt_m2_10000_tor4_pr45_100Myr_lower_pdf'

    try:
        os.mkdir(Default_dir)
    except:
        pass

    try:
        os.mkdir(Default_dir_pdf)
    except:
        pass
    
    plt.savefig(Default_dir + '/Fire' + str(sorted[oo]).zfill(3) + '.png', dpi = 600, transparent = True, bbox_inches='tight')
    plt.savefig(Default_dir_pdf + '/Fire' + str(sorted[oo]).zfill(3) + '.pdf', dpi = 600, transparent = True, bbox_inches='tight')

    plt.close()
    oo += 1
    
    #pdb.set_trace()
