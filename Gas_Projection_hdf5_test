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

"""
Here we present some explanations.

Camera:
    The Camera class is a container that stores the camera parameters. The camera is an object that lives in the space and has spherical coordinates (r,theta,phi), centred around the location (x,y,z). Angles *theta* and *phi* are given in degrees, and enable to rotate the camera along the x-axis and y-axis, respectively. The *roll* angle induces rotations along the line-of-sight, i.e., the z-axis.
"""


import matplotlib.pyplot as plt
import h5py
from sphviewer.tools import QuickView

### Interesting toturial: http://alejandrobll.github.io/py-sphviewer/content/tutorial_projections.html

### Here we study the density profile of the gas particles. We consider two different cases. 1) When we take into account all of the gas particles and monitor their density profile. 2) When we consider the star forming gas and see how do they look like!!!


Dir_path = '/n/holylfs05/LABS/hernquist_lab/AGN_Feedback_Fire/m12_mcvt_m2_t95_3000_tor4/output/snapshot_056.hdf5'

data = h5py.File("/n/holylfs05/LABS/hernquist_lab/AGN_Feedback_Fire/m12_mcvt_m2_t95_3000_tor4/output/snapshot_056.hdf5","r")

#data.keys()

#data["PartType0"].keys()


"""
<KeysViewHDF5 ['CoolingRate', 'Coordinates', 'Density', 'DivBcleaningFunctionGradPhi', 'DivBcleaningFunctionPhi', 'DivergenceOfMagneticField', 'ElectronAbundance', 'HeatingRate', 'HydroHeatingRate', 'IMFFormationProperties', 'InternalEnergy', 'MagneticField', 'Masses', 'MetalCoolingRate', 'Metallicity', 'NetHeatingRateQ', 'NeutralHydrogenAbundance', 'ParticleChildIDsNumber', 'ParticleIDGenerationNumber', 'ParticleIDs', 'SmoothingLength', 'StarFormationRate', 'Velocities']>


"""


for uu in range(1):
    
    BH_Center = data["PartType5"]["Coordinates"][:]
    Gas_location = data["PartType0"]["Coordinates"][:] - BH_Center
    Gas_mass = 1e10*data["PartType0"]["Masses"][:]
    Gas_Softening = data["PartType0"]["SmoothingLength"][:]
    
    NN = 100
    
    #mass = np.ones(len(gas_Loc[:,0]))
    hh = Gas_Softening
    mass = Gas_mass
    Particles = sphviewer.Particles(Gas_location, mass, hh)

    Scene = sphviewer.Scene(Particles)
    
    fig = plt.figure(1,figsize=(15,5))
    
    #set a figure title on top
    #fig.suptitle("Star Forming Gas:" + str(labb_tit[uu]), fontsize=17, x=0.5, y=1.4)
    
    fig.suptitle(r"Precessing kinetic jet with low energy flux: Snapshot #" + str(sorted[oo]+1) + " with Gas_mass and O6_ion_Fraction", fontsize=17, x=0.5, y=2.5)
    plt.subplots_adjust(top =2.5, bottom=0.2, hspace=0.3, wspace=0.3)
    
    #plt.subplots_adjust(top =1.8, bottom=0.2,hspace=0.3, wspace=0.3)
    
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
    
    plt.savefig('/n/home13/nqutob/AGN_Feedback/hdf5_test/snapshot_056.png', dpi = 400, transparent = True,bbox_inches='tight')
 
    
pdb.set_trace()
