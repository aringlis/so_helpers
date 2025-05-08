from sunpy import map
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
from sunpy.map.maputils import all_coordinates_from_map
import matplotlib.pyplot as plt
import numpy as np

def solo_phi_quicklook_rotate(phi_quicklook_file):

    # do some string manipulation to create a new filename
    name_comp = phi_quicklook_file.split('_')
    newname = name_comp[0] +'_' + name_comp[2] + '_iswa_quicklook_' + name_comp[3] + '_' + name_comp[4] + '_' + name_comp[5].split('.')[0] + '.png'

    # load in the image
    phimap = map.Map(phi_quicklook_file)
    # rotate the image so solar north is up
    phimap2 = phimap.rotate()

    # extract a sub_map from the original image
    rsun = phimap2.meta['rsun_arc']
    top_right = SkyCoord((rsun+200) * u.arcsec, (rsun+200)*u.arcsec,frame = phimap2.coordinate_frame)
    bottom_left = SkyCoord(((rsun*-1)-300) * u.arcsec, ((rsun*-1)-300)* u.arcsec,frame = phimap2.coordinate_frame)

    phisubmap = phimap2.submap(bottom_left, top_right = top_right)

    # create a mask for anything above the limb
    pixel_coords = all_coordinates_from_map(phisubmap)
    solar_center = SkyCoord(0*u.deg, 0*u.deg, frame=phisubmap.coordinate_frame)
    pixel_radii = np.sqrt((pixel_coords.Tx-solar_center.Tx)**2 +
                      (pixel_coords.Ty-solar_center.Ty)**2)

    offlimb_mask = pixel_radii > phisubmap.rsun_obs
    new_cmap = phisubmap.cmap.copy()
    new_cmap.set_bad('black')
    phisubmap.mask = offlimb_mask

    # make the plot
    plt.figure(1,figsize=(10,9))
    plt.subplots_adjust(bottom = 0.1, top = 0.95, left = 0.1, right = 0.95)
    phisubmap.plot(cmap = new_cmap)
    phisubmap.draw_grid()
    cbar = plt.colorbar(shrink=0.9)
    cbar.ax.tick_params(labelsize=12)
    if phisubmap.meta['phidtype'] == 'blos':
        plt.clim([-200,200])
        cbar.set_label('B (Gauss)',size=14)
    plt.tick_params(labelsize=12)
    plt.xlabel('Helioprojective longitude (solar-x)',fontsize=14)
    plt.ylabel('Helioprojective latitude (solar-y)',fontsize=14)

    if phisubmap.meta['phidtype'] == 'blos':
        newtitle = phisubmap.meta['telescop'] + ' LL02 magnetogram: ' + phisubmap.meta['date-obs']
    elif phisubmap.meta['phidtype'] == 'icnt':
        newtitle = phisubmap.meta['telescop'] + ' LL02 continuum intensity: ' + phisubmap.meta['date-obs']
    plt.title(newtitle,fontsize=14)
    
    lon = phimap2.meta['hgln_obs']
    lat = phimap2.meta['hglt_obs']
    dsun = phimap2.meta['dsun_au']
    
    plt.figtext(0.12,0.22,'S/C info:',fontsize=14,color='white')
    plt.figtext(0.12,0.20,'HG lon = ' + str(round(lon,1)) + '$^{\circ}$',fontsize=14,color='white')
    plt.figtext(0.12,0.18,'HG lat = ' + str(round(lat,1)) + '$^{\circ}$',fontsize=14,color='white')
    plt.figtext(0.12,0.16,'$D_{sun}$ = ' + str(round(dsun,2)) + ' AU',fontsize=14,color='white')

    # add ESA copyright label
    plt.figtext(0.48,0.16,'\u00A9 ESA & NASA / Solar Orbiter / PHI team',fontsize=10,color='white')

    #plt.figtext(0.15,0.7,'Quicklook data',fontsize=30,rotation=30,color='red',alpha=0.4)
    #plt.figtext(0.48,0.2,'Quicklook data',fontsize=30,rotation=30,color='red',alpha=0.4)

    # plot the zero longitude line for reference
    stonyhurst_frame = frames.HeliographicStonyhurst(obstime=phisubmap.date)
    constant_lon = SkyCoord(0.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)
    ax = plt.gca()
    ax.plot_coord(constant_lon, color="lightblue")
    
    
    plt.savefig(newname)
    plt.close()
    

