from sunpy import map
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
from sunpy.map.maputils import all_coordinates_from_map
from sunpy.visualization import drawing
import matplotlib.pyplot as plt
import numpy as np

import urllib
import bs4
import requests
import datetime
import os

def solo_phi_quicklook_rotate(phi_quicklook_file):

    # do some string manipulation to create a new filename
    name_comp = phi_quicklook_file.split('_')
    newname = name_comp[0] +'_' + name_comp[2] + '_iswa_lowlatency_' + name_comp[3] + '.png'#'_' + name_comp[4] + '_' + name_comp[5].split('.')[0] + '.png'

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

    # make the plot. Make the entire background black with no frame
    plt.figure(1,figsize=(10,10),frameon=False)
   # plt.subplots_adjust(bottom = 0.1, top = 0.95, left = 0.1, right = 0.95)
    plt.subplots_adjust(bottom = 0.0, top = 1.0, left = 0.0, right = 1.0)
    phisubmap.plot(cmap = new_cmap)
    ax = plt.gca()
    ax.grid(False)
    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[0].set_ticks_visible(False)
    ax.coords[1].set_ticklabel_visible(False)
    ax.coords[1].set_ticks_visible(False)
    phisubmap.draw_grid()
    
   # cbar = plt.colorbar(shrink=0.9)
   # cbar.ax.tick_params(labelsize=12)
    if phisubmap.meta['phidtype'] == 'blos':
        plt.clim([-200,200])
  #      cbar.set_label('B (Gauss)',size=14)
    plt.tick_params(labelsize=12)
   # plt.xlabel('Helioprojective longitude (solar-x)',fontsize=14)
   # plt.ylabel('Helioprojective latitude (solar-y)',fontsize=14)

    # add some informational text to the plot
    lon = phimap2.meta['hgln_obs']
    lat = phimap2.meta['hglt_obs']
    dsun = phimap2.meta['dsun_au']
    
    plt.figtext(0.02,0.92,'S/C info:',fontsize=14,color='white')
    plt.figtext(0.02,0.90,'HG lon = ' + str(round(lon,1)) + '$^{\circ}$',fontsize=14,color='white')
    plt.figtext(0.02,0.88,'HG lat = ' + str(round(lat,1)) + '$^{\circ}$',fontsize=14,color='white')
    plt.figtext(0.02,0.86,'$D_{sun}$ = ' + str(round(dsun,2)) + ' AU',fontsize=14,color='white')

    if phisubmap.meta['phidtype'] == 'blos':
        plt.figtext(0.02,0.05,phisubmap.meta['telescop'] + ' LL02 magnetogram',fontsize=14,color='white')
        plt.figtext(0.02,0.02,phisubmap.meta['date-obs'],fontsize=14,color='white')
    elif phisubmap.meta['phidtype'] == 'icnt':
        plt.figtext(0.02,0.05,phisubmap.meta['telescop'] + ' LL02 continuum intensity',fontsize=14,color='white')
        plt.figtext(0.02,0.025,phisubmap.meta['date-obs'],fontsize=14,color='white')

    # add ESA copyright label
    plt.figtext(0.65,0.03,'\u00A9 ESA & NASA / Solar Orbiter / PHI team',fontsize=10,color='white')

    # plot the zero longitude line for reference - but only at points where it is visible!
    stonyhurst_frame = frames.HeliographicStonyhurst(obstime=phisubmap.date)
    constant_lon = SkyCoord(0.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)
    ax = plt.gca()
    # hide the zero longitude line when it's behind the limb
    visible, hidden = drawing._plot_vertices(constant_lon, ax, phisubmap.coordinate_frame, phisubmap.rsun_meters,color='blue',close_path = False)
    drawing._modify_polygon_visibility(hidden, np.zeros(len(constant_lon)).astype('bool'))
    #ax.plot_coord(constant_lon, color="lightblue")

    # save the plot as a png
    if phisubmap.meta['phidtype'] == 'blos':
        plt.savefig(os.path.join('blos',newname))
    elif phisubmap.meta['phidtype'] == 'icnt':
        plt.savefig(os.path.join('cont',newname))
    else:
        plt.savefig(newname,bbox_inches ='tight')
    plt.close()
    


def find_latest_phi_ll02_files(date = datetime.datetime.today().date()):
    # PHI LL02 files live at the following webpage:
    basepage = 'https://www2.mps.mpg.de/data/outgoing/valori/PHI_LL/'
 
    date_str = date.strftime('%Y-%m-%d')

    # parse the page using BeautifulSoup
    page = requests.get(basepage + date_str)
    soup = bs4.BeautifulSoup(page.text)

    # extract the links to the FITS files
    links = []
    for link in soup.find_all('a', href = True):
        if '.fits' in link['href']:
            links.append(link['href'])

    #download the FITS files
    for link in links:
        full_fname = basepage + date_str + '/' + link
        urllib.request.urlretrieve(full_fname, link)



def solo_phi_quicklook_rotate_gherardo(phi_quicklook_file):
    '''create an example plot for Gherardo.'''
    
    # do some string manipulation to create a new filename
    name_comp = phi_quicklook_file.split('_')
    newname = name_comp[0] +'_' + name_comp[2] + '_iswa_lowlatency_' + name_comp[3] + '_' + name_comp[4] + '_' + name_comp[5].split('.')[0] + '.png'

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
    ax = plt.gca()
    ax.grid(False)
    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[0].set_ticks_visible(False)
    ax.coords[1].set_ticklabel_visible(False)
    ax.coords[1].set_ticks_visible(False)
    phisubmap.draw_grid()
    
    cbar = plt.colorbar(shrink=0.9)
    cbar.ax.tick_params(labelsize=12)
    if phisubmap.meta['phidtype'] == 'blos':
        plt.clim([-200,200])
        cbar.set_label('B (Gauss)',size=14)
    plt.tick_params(labelsize=12)
#    plt.xlabel('Helioprojective longitude (solar-x)',fontsize=14)
 #   plt.ylabel('Helioprojective latitude (solar-y)',fontsize=14)

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

    constant_lon_45 = SkyCoord(45.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)

    constant_lon_90 = SkyCoord(90.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)

    constant_lon_135 = SkyCoord(135.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)

    constant_lon_180 = SkyCoord(180.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)

    constant_lon_225 = SkyCoord(225.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)

    constant_lon_270 = SkyCoord(270.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)

    constant_lon_315 = SkyCoord(315.0*u.deg, np.linspace(-90, 90, 100) * u.deg,
                                frame=stonyhurst_frame)
    
    ax = plt.gca()
    ax.plot_coord(constant_lon, color="blue", label = '0')
   # ax.plot_coord(constant_lon_45, color="tab:blue", label = '45')
   # ax.plot_coord(constant_lon_90, color="tab:orange", label = '90')
    ax.plot_coord(constant_lon_135, color="tab:orange", label = '135')
    ax.plot_coord(constant_lon_180, color="tab:red", label = '180')
    ax.plot_coord(constant_lon_225, color="tab:purple", label = '225')
    #ax.plot_coord(constant_lon_270, color="tab:pink", label = '270')
   # ax.plot_coord(constant_lon_315, color="tab:olive", label = '315')

    plt.legend()
    
    
  #  if phisubmap.meta['phidtype'] == 'blos':
   #     plt.savefig(os.path.join('blos',newname))
   # elif phisubmap.meta['phidtype'] == 'icnt':
    #    plt.savefig(os.path.join('cont',newname))
    #else:
    plt.savefig(newname)
    plt.close()
    
