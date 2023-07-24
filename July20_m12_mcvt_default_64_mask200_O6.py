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
# you could put a for loop here that replaces 'm12_mcvt_default_64 ' with whaterver directories are in /hernquist_lab/AGN_Feedback_Fire then goes to the 'output' folder in each 
print("Checkpoint 0")
Directory_path = '/n/holylfs05/LABS/hernquist_lab/AGN_Feedback_Fire/m12_mcvt_default_64/output/'

diir = glob(Directory_path + '*.hdf5')

List_diir = np.atleast_1d(diir)
snap_order = []
List_diir_image = []

List_diir = np.atleast_1d(diir)
snap_order = []
for ii in range(len(List_diir)): 
    AA = List_diir[ii]
    snap_ii = int(AA[AA.find('snapshot_')+len('snapshot_'):AA.find('.hdf5')])
    snap_order = np.append(snap_order, [snap_ii])
    List_diir_image = np.append(List_diir_image, List_diir[ii])
    # add txt conversion line later! 
    print("Checkpoint 0B")
# this puts the snapshots in the correct numerical order when they are runa nd saved 
int_snap = snap_order.astype(int) 
sorted = np.sort(int_snap)
List_diir = List_diir_image[np.argsort(int_snap)]
print("Checkpoint 0C")

oo = 0
for fname in List_diir:
    data = h5py.File(fname,"r")
    BH_Center = data["PartType5"]["Coordinates"][:]
    print(BH_Center)
    BH_Center = data["PartType5"]["Coordinates"][:]
    Gas_location = data["PartType0"]["Coordinates"][:] - BH_Center
    Gas_mass = 1e10*data["PartType0"]["Masses"][:]
    Gas_Softening = data["PartType0"]["SmoothingLength"][:]
    print("Checkpoint 1")
    
    # this is where I started adding additional parameters
    ds = yt.load(fname, unit_base=unit_base) 
    #ds.derived_field_list # This displays all the different parameters that I can customize 
    ions_names=['H I','Mg II','C IV','N V','O VI', 'O VII', 'O VIII', 'Ne VIII']  ## Exploration!
    # I think that if we include the number (I, II, III etc) then it only analyzes that singlular ion, but if you include just the atom name (H, C, O, etc.) then it will inlcude all ions from that element 
    trident.add_ion_fields(ds, ions=ions_names, ftype="PartType0") #PartType0 is for GIZMO and GADGET 

    mass = Gas_mass
    Oxygen5_mass = ds.all_data()[('gas', 'O_p5_mass')].in_units('Msun') # should give Msun/kpc^2 
    Oxygen6_mass = ds.all_data()[('gas', 'O_p6_mass')].in_units('Msun')
    Oxygen7_mass = ds.all_data()[('gas', 'O_p7_mass')].in_units('Msun')
    # Magnesium goes here
    print("Checkpoint 2")
    
    NN = 100 # this defines the x, y, and z axis ranges for the plots 
    hh = Gas_Softening #what does gas softening mean?

    Particles = sphviewer.Particles(Gas_location, Oxygen5_mass, hh) # CHANGE PARAMETER HERE!
    Scene = sphviewer.Scene(Particles)
    
    fig = plt.figure(1,figsize=(15,5))
    fig.suptitle(r"No Jet : Snapshot #" + str(sorted[oo]+1) + " with O6 Mass & Mask 200", fontsize=17, x=0.5, y=1.5) #set a figure title on top
    plt.subplots_adjust(top =1.8, bottom=0.2, hspace=0.3, wspace=0.3)
    
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    extendd = [-NN,NN,-NN,NN]
    
    # Begin Masking
    z = Gas_location[:,2] 
    mask_z = np.abs(z)<200 # 5 is the parameter that we can modify 
    Particles1 = sphviewer.Particles(Gas_location[mask_z], Oxygen5_mass[mask_z], hh[mask_z]) # CHANGE PARAMETER HERE! 

    Scene1 = sphviewer.Scene(Particles1) 
    Scene1.update_camera(r='infinity', t=0, p = 0, roll = 0, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent=extendd) 
    
    Render1 = sphviewer.Render(Scene1)
    Render1.set_logscale()
    img1 = Render1.get_image()
    extent1 = Render1.get_extent()
    divider = make_axes_locatable(ax1)

    #ax1.imshow(img1, extent=extent1, origin='lower', cmap=plt.cm.jet, vmax= 3.0, rasterized=True)
    image1 = ax1.imshow(img1, extent=extent1, origin='lower', cmap=plt.cm.jet, rasterized=True, vmin=0, vmax=4)
    cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    fig.add_axes(cax)
    cb = fig.colorbar(image1, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=15)

    ax1.set_xlabel('$X$(kpc)', size=12)
    ax1.set_ylabel('$Y$(kpc)', size=12)

    x = Gas_location[:,0] 
    mask_x=np.abs(x)<200  # 5 is the parameter that we can modify 
    Particles2 = sphviewer.Particles(Gas_location[mask_x], Oxygen5_mass[mask_x], hh[mask_x]) # CHANGE PARAMETER HERE! 

    Scene2 = sphviewer.Scene(Particles2)
    Scene2.update_camera(r='infinity', t=-90, p = -90, roll = 0, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent= extendd)
    
    Render2 = sphviewer.Render(Scene2)
    Render2.set_logscale()
    img2 = Render2.get_image()
    extent2 = Render2.get_extent()
    #divider = make_axes_locatable(ax2)
    #ax2.imshow(img2, extent=extent2, origin='lower',cmap=plt.cm.jet, vmax= 3.0, rasterized=True)
    image2 = ax2.imshow(img2, extent=extent2, origin='lower',cmap=plt.cm.jet, rasterized=True, vmin=0, vmax=4)

    divider = make_axes_locatable(ax2)
    cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    fig.add_axes(cax)
    cb = fig.colorbar(image2, cax=cax, orientation="horizontal")
    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    cb.ax.tick_params(labelsize=15)

    ax2.set_xlabel('$Y$(kpc)', size=12)
    ax2.set_ylabel('$Z$(kpc)', size=12)

    y = Gas_location[:,1]
    mask_y=np.abs(y)<200  #5 is the parameter that we can modify 
    Particles3 = sphviewer.Particles(Gas_location[mask_y], Oxygen5_mass[mask_y], hh[mask_y]) # CHANGE PARAMETER HERE! 

    Scene3 = sphviewer.Scene(Particles3) # MODIFIED!
    Scene3.update_camera(r='infinity', t=90, p = 0, roll = -90, x = 0, y = 0, z = 0, vmin= 6.3, vmax= 7.4, extent= extendd)

    Render3 = sphviewer.Render(Scene3)
    Render3.set_logscale()
    img3 = Render3.get_image()
    extent3 = Render3.get_extent()
    divider = make_axes_locatable(ax3)
    #ax3.imshow(img3, extent=extent3, origin='lower', cmap=plt.cm.jet, vmax= 3.0, rasterized=True)
    image3 = ax3.imshow(img3, extent=extent2, origin='lower',cmap=plt.cm.jet, rasterized=True, vmin=0, vmax=4)

    cax = divider.new_vertical(size="7%", pad=0.7, pack_start=True)
    fig.add_axes(cax)
    cb = fig.colorbar(image3, cax=cax, orientation="horizontal")
    #cb.set_label(label='Temperature ($^{\circ}$C)', size='large', weight='bold')
    cb.ax.tick_params(labelsize=15)

    ax3.set_xlabel('$Z$(kpc)', size=12)
    ax3.set_ylabel('$X$(kpc)', size=12)
    
    Default_dir = '/n/home13/nqutob/AGN_Feedback/ion_snapshots/' + 'm12_mcvt_default_64_mask200_O6'
    #Default_dir_pdf = '/n/home13/nqutob/AGN_Feedback/ion_snapshots/' + 'm12_mcvt_default_64_mask200_pdf'

    try:
        os.mkdir(Default_dir)
    except:
        pass

    #try:
    #    os.mkdir(Default_dir_pdf)
    #except:
    #    pass
    
    plt.savefig(Default_dir + '/Fire' + str(sorted[oo]).zfill(3) + '_m12_mcvt_default_64_mask200_O6' + '.png', dpi = 600, transparent = True, bbox_inches='tight')
    #plt.savefig(Default_dir_pdf + '/Fire' + str(sorted[oo]).zfill(3) + '_m12_mcvt_default_64_mask200_O6' + '.pdf', dpi = 600, transparent = True, bbox_inches='tight')

    plt.close()
    oo += 1
    
    #pdb.set_trace()
