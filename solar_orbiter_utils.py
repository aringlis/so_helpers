from sunpy.coordinates import frames, get_horizons_coord, get_body_heliographic_stonyhurst
from astropy.coordinates import SkyCoord
import astropy.units as u
import datetime
import matplotlib.pyplot as plt
import numpy as np


def where_is_solar_orbiter(obstime = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M'), parker = True, stereo = True):

    so_location = get_horizons_coord('solar orbiter', obstime)
    so_x = np.cos(np.deg2rad(so_location.lon.value)) * -1.0 * so_location.radius.value
    so_y = np.sin(np.deg2rad(so_location.lon.value)) * -1.0 * so_location.radius.value

    if parker:
        parker_location = get_horizons_coord('parker solar probe', obstime)
        parker_x = np.cos(np.deg2rad(parker_location.lon.value)) * -1.0 * parker_location.radius.value
        parker_y = np.sin(np.deg2rad(parker_location.lon.value)) * -1.0 * parker_location.radius.value

    if stereo:
        stereo_location = get_horizons_coord('stereo-A')
        stereo_x = np.cos(np.deg2rad(stereo_location.lon.value)) * -1.0 * stereo_location.radius.value
        stereo_y = np.sin(np.deg2rad(stereo_location.lon.value)) * -1.0 * stereo_location.radius.value

    #make a simple plot

    plt.figure(1,figsize=(8,8))

    plt.xlim([-1.2,1.2])
    plt.ylim([-1.2,1.2])

    plt.title('Solar Orbiter location at: ' + obstime)
    plt.xlabel('Distance (AU)')
    plt.ylabel('Distance (AU)')
    earth = plt.Circle((-1.0, 0.0),0.03, fill=True, facecolor='slateblue', alpha = 0.7, label='Earth')
    sun = plt.Circle((0.0, 0.0),0.06, fill=True, facecolor='orange', alpha = 0.7, label='Sun')


    solar_orbiter = plt.Circle((so_x, so_y),0.02, fill=True, facecolor='black', alpha = 1.0, label='Solar Orbiter')

    if parker:
        parker =  plt.Circle((parker_x, parker_y),0.02, fill=True, facecolor='purple', alpha = 0.7, label='Parker Solar Probe')
    if stereo:
        stereo = plt.Circle((stereo_x, stereo_y),0.02, fill=True, facecolor='green', alpha = 0.7, label='Stereo-A')

    
    ax = plt.gca()
    ax.add_patch(earth)
    ax.add_patch(sun)
    ax.add_patch(solar_orbiter)
    if parker:
        ax.add_patch(parker)
    if stereo:
        ax.add_patch(stereo)

    plt.axvline(0.0,linestyle='dashed',color='lightgrey')
    plt.axhline(0.0,linestyle='dashed',color='lightgrey')

    radius02 = plt.Circle((0.0, 0.0),0.2, fill=False, facecolor=None, linestyle='dashed',color='lightgrey')
    radius04 = plt.Circle((0.0, 0.0),0.4, fill=False, facecolor=None, linestyle='dashed',color='lightgrey')
    radius06 = plt.Circle((0.0, 0.0),0.6, fill=False, facecolor=None, linestyle='dashed',color='lightgrey')
    radius08 = plt.Circle((0.0, 0.0),0.8, fill=False, facecolor=None, linestyle='dashed',color='lightgrey')
    radius10 = plt.Circle((0.0, 0.0),1.0, fill=False, facecolor=None, linestyle='dashed',color='lightgrey')

    ax.add_patch(radius02)
    ax.add_patch(radius04)
    ax.add_patch(radius06)
    ax.add_patch(radius08)
    ax.add_patch(radius10)

    plt.figtext(0.65,0.25,'Spacecraft location (HG):')
    plt.figtext(0.75,0.2,'Lon: ' + str(np.round(so_location.lon.degree,2)))
    plt.figtext(0.75,0.17,'Lat: ' + str(np.round(so_location.lat.degree,2)))
    plt.figtext(0.75,0.14,'Radius: ' + str(np.round(so_location.radius,2)))

    plt.legend()

    plt.gca().set_aspect('equal')
    plt.show()

    


    
    

    

    

