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


Dir_path = '/n/holylfs05/LABS/hernquist_lab/AGN_Feedback_Fire/m12_mcvt_m2_t95_3000_tor4/output/snapshot_057.hdf5'

data = h5py.File("/n/holylfs05/LABS/hernquist_lab/AGN_Feedback_Fire/m12_mcvt_m2_t95_3000_tor4/output/snapshot_057.hdf5","r")

#data.keys()

#data["PartType0"].keys()


"""
<KeysViewHDF5 ['CoolingRate', 'Coordinates', 'Density', 'DivBcleaningFunctionGradPhi', 'DivBcleaningFunctionPhi', 'DivergenceOfMagneticField', 'ElectronAbundance', 'HeatingRate', 'HydroHeatingRate', 'IMFFormationProperties', 'InternalEnergy', 'MagneticField', 'Masses', 'MetalCoolingRate', 'Metallicity', 'NetHeatingRateQ', 'NeutralHydrogenAbundance', 'ParticleChildIDsNumber', 'ParticleIDGenerationNumber', 'ParticleIDs', 'SmoothingLength', 'StarFormationRate', 'Velocities']>


"""


for uu in range(1):

    ds = yt.load(Dir_path, unit_base=unit_base)
  
    BH_Center = data["PartType5"]["Coordinates"][:]
    Gas_location = data["PartType0"]["Coordinates"][:] - BH_Center
    Gas_mass = 1e10*data["PartType0"]["Masses"][:]
    Gas_Softening = data["PartType0"]["SmoothingLength"][:]

    ions_names=['H I','Mg II','C IV','N V','O VI', 'O VII', 'O VIII', 'Ne VIII']  ## Exploration!
    trident.add_ion_fields(ds, ions=ions_names, ftype="PartType0") #PartType0 is for GIZMO and GADGET 

    Oxygen6_mass = ds.all_data()[('gas', 'O_p6_mass')].in_units('Msun')
    
    NN = 100
    hh = Gas_Softening

    Particles = sphviewer.Particles(Gas_location, Oxygen6_mass, hh)
    Scene = sphviewer.Scene(Particles)
    
    fig = plt.figure(1,figsize=(15,5))
    
    fig.suptitle(r"Precessing kinetic jet with low energy flux: Snapshot #057 with Gas_mass and O6_ion_Fraction", fontsize=17, x=0.5, y=1.6)
    plt.subplots_adjust(top =1.8, bottom=0.2, hspace=0.3, wspace=0.3)
    
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    extendd = [-NN,NN,-NN,NN]

    z = Gas_location[:,2] # MASK!!
    mask_z=np.abs(z)<5  # MASK!!
    Particles1 = sphviewer.Particles(Gas_location[mask_z], Oxygen6_mass[mask_z], hh[mask_z]) # MASK!!

    Scene1 = sphviewer.Scene(Particles1) # MASK!! 
    Scene.update_camera(r='infinity', t=0, p = 0, roll = 0, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent=extendd) 
    
    Render1 = sphviewer.Render(Scene1)
    Render1.set_logscale()
    img1 = Render1.get_image()
    extent1 = Render1.get_extent()
    divider = make_axes_locatable(ax1)

    #ax1.imshow(img1, extent=extent1, origin='lower', cmap=plt.cm.jet, vmax= 3.0, rasterized=True)
    image1 = ax1.imshow(img1, extent=extent1, origin='lower', cmap=plt.cm.jet, rasterized=True)
    #divider = make_axes_locatable(ax1)
    cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    fig.add_axes(cax)
    cb = fig.colorbar(image1, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=15)

    ax1.set_xlabel('$X$(kpc)', size=12)
    ax1.set_ylabel('$Y$(kpc)', size=12)


    x = Gas_location[:,2] # MASK!!
    mask_x=np.abs(x)<5  #5 is the parameter that we can modify 
    Particles2 = sphviewer.Particles(Gas_location[mask_x], Oxygen6_mass[mask_x], hh[mask_x]) 

    Scene2 = sphviewer.Scene(Particles2) # MODIFIED!
    Scene.update_camera(r='infinity', t=-90, p = -90, roll = 0, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent= extendd)
    
    Render2 = sphviewer.Render(Scene2)
    Render2.set_logscale()
    img2 = Render2.get_image()
    extent2 = Render2.get_extent()
    #divider = make_axes_locatable(ax2)
    #ax2.imshow(img2, extent=extent2, origin='lower',cmap=plt.cm.jet, vmax= 3.0, rasterized=True)
    image2 = ax2.imshow(img2, extent=extent2, origin='lower',cmap=plt.cm.jet, rasterized=True)

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

    y = Gas_location[:,2]
    mask_y=np.abs(y)<5  #5 is the parameter that we can modify 
    Particles3 = sphviewer.Particles(Gas_location[mask_y], Oxygen6_mass[mask_y], hh[mask_y]) 

    Scene3 = sphviewer.Scene(Particles3) # MODIFIED!
    Scene.update_camera(r='infinity', t=90, p = 0, roll = -90, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent= extendd)

    Render = sphviewer.Render(Scene3)
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

    #cb = fig.colorbar(img3, cax=cax, orientation="horizontal")
    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    #cb.ax.tick_params(labelsize=15)
    #ax1.set_title(labb_tit[uu])
    #ax3.set_title(labb_tit[uu])
    ax3.set_xlabel('$Z$(kpc)', size=12)
    ax3.set_ylabel('$X$(kpc)', size=12)
    
    plt.savefig('/n/home13/nqutob/AGN_Feedback/hdf5_test/snapshot_test.png', dpi = 400, transparent = True,bbox_inches='tight')
 
    
pdb.set_trace()
