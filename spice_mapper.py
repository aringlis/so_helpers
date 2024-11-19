import pandas
import matplotlib.pyplot as plt
import numpy as np
import datetime


def spice_mapper(date, catalog_file = 'spice_catalog.csv'):
    '''
    This tool makes a representative plot of all SPICE observations on a given day.
    The SPICE catalog file is needed to use this tool.'''

    if type(date) != datetime.datetime:
        raise ValueError('Input date must be in datetime format, e.g. date = datetime.datetime(2024,3,22)')
    
    
    # read in the SPICE catalog and add two new time columns to enable easy searching
    spicecat = pandas.read_csv(catalog_file)

    start_times = []
    end_times = []
    for i, s in spicecat.iterrows():
        start_times.append(datetime.datetime.fromisoformat(s['DATE-BEG']))
        end_times.append(datetime.datetime.fromisoformat(s['DATE-END']))
    spicecat['DATE-BEG-ISO'] = start_times
    spicecat['DATE-END-ISO'] = end_times

    # now, find and display all the observation on the selected date
    entries = spicecat[(spicecat['DATE-BEG-ISO'] > date) &
                           (spicecat['DATE-END-ISO'] < date + datetime.timedelta(days=1))]

    if len(entries) == 0:

        print(' ')
        print('--------------')
        print('No observations taken by SPICE on ' + date.isoformat())
        print('--------------')
        
        plt.figure(1,figsize=(18,8))
        plt.subplots_adjust(left = 0.06, right=0.5,top=0.95,bottom=0.08)
        plt.title('SPICE observations on ' + date.isoformat())
        ax = plt.gca()
        ax.axis('equal')
        plt.xlim([-3500,3500])
        plt.ylim([-3500,3500])
        plt.figtext(0.12,0.3,'NO SPICE OBSERVATIONS',fontsize=30,fontweight='bold',rotation=-30)
        plt.savefig('spice_observations_on_' + date.strftime('%Y%m%d') + '.png',dpi=300)
        plt.show()
        return
    
    # put the valid entries in order by start time
    entries.sort_values('DATE-BEG-ISO')

    # plot a simple representation of the Sun
    # ---------------------------------------------
    #solar radius in m
    solar_radius = 6.957e8
    #1AU in m
    astronomical_unit = 1.495e11
    distance_from_sun = entries.iloc[0]['DSUN_AU'] * astronomical_unit
    
    #apparent radius of Sun in arcsec
    rsun_apparent = np.rad2deg(np.arctan(solar_radius / distance_from_sun)) * 3600.

    #make a simple plot of the Sun
    plt.figure(1,figsize=(18,8))
    plt.subplots_adjust(left = 0.06, right=0.5,top=0.95,bottom=0.08)
    circle1 = plt.Circle((0, 0), rsun_apparent, fill=True, facecolor='orange',alpha=0.5)

    ax = plt.gca()
    ax.add_patch(circle1)
    plt.set_cmap('tab20')

    #plot each SPICE observation on the map of the Sun
    #---------------------------------

    #don't plot every single file. Instead, plot representations of each unique SPIOBSID
    #i.e. a sequence of sit and stares should have the same SPIOBSID
    unique_spiobsids = entries['SPIOBSID'].unique()

    for id in unique_spiobsids:
        subentries = entries[entries['SPIOBSID'] == id]
        
        if len(subentries) == 1:
            first_entry = subentries.iloc[0]
            lenx = first_entry['NAXIS1'] * first_entry['CDELT1']
            leny = first_entry['NAXIS2'] * first_entry['CDELT2']
            x,y = calculate_corners(first_entry['CRVAL1'],first_entry['CRVAL2'],lenx,leny,first_entry['CROTA'],plot=False)
            description = first_entry['SOOPNAME'] + ' | ' + first_entry['STUDYTYP'] + ' | ' + first_entry['DATE-BEG'] + ' - ' + first_entry['DATE-END'] 
            plt.plot(x,y,label=description)
        elif len(subentries) > 1:
            num_repeats = len(subentries)
            first_entry = subentries.iloc[0]
            last_entry = subentries.iloc[-1]
            lenxfirst = first_entry['NAXIS1'] * first_entry['CDELT1']
            lenyfirst = first_entry['NAXIS2'] * first_entry['CDELT2']
            lenxlast = last_entry['NAXIS1'] * last_entry['CDELT1']
            lenylast = last_entry['NAXIS2'] * last_entry['CDELT2']
            xfirst,yfirst = calculate_corners(first_entry['CRVAL1'],first_entry['CRVAL2'],lenxfirst,lenyfirst,first_entry['CROTA'],plot=False)
            xlast,ylast = calculate_corners(last_entry['CRVAL1'],last_entry['CRVAL2'],lenxlast,lenylast,last_entry['CROTA'],plot=False)
            description = first_entry['SOOPNAME'] + ' | ' + first_entry['STUDYTYP'] + ' (x'+str(num_repeats) + ') | ' + first_entry['DATE-BEG'] + ' - ' + last_entry['DATE-END']
            plt.plot([xfirst[0],xlast[1],xlast[2],xfirst[3],xfirst[4]],[yfirst[0],ylast[1],ylast[2],yfirst[3],yfirst[4]],label=description)
                


  #  plt.xlim([-rsun_apparent,rsun_apparent])
  #  plt.ylim([-rsun_apparent,rsun_apparent])
    
    plt.xlabel('x (arcsec)')
    plt.ylabel('y (arcsec)')
    plt.title('SPICE observations on ' + date.isoformat())

    
    ax = plt.gca()
    
    ax.axis('equal')
    plt.xlim([-3500,3500])
    plt.ylim([-3500,3500])
    plt.grid(color='grey',linestyle='dashed',alpha=0.5)
    ax.legend(bbox_to_anchor=(1.02,0.95),loc='upper left',fontsize = 9,fancybox = True, shadow = True)


    plt.figtext(0.09,0.14,'Solar Orbiter position at: ' + entries.iloc[0]['DATE-BEG'])
    plt.figtext(0.09,0.12,'$D_{sun}$ = ' + str(round(entries.iloc[0]['DSUN_AU'],2)) + ' AU')
    plt.figtext(0.09,0.10,'HG lon = ' + str(round(entries.iloc[0]['HGLN_OBS'],1)) + '$^{\circ}$')

    plt.savefig('spice_observations_on_' + date.strftime('%Y%m%d') + '.png',dpi=300)
    plt.show()



def calculate_corners(crval1,crval2,lenx,leny,crota, plot = False):
    '''Figure out the corners of a FOV box given a rotation angle and other parameters.'''

    mod1 = np.sin(np.deg2rad(crota))
    mod2 = np.cos(np.deg2rad(crota))

    #top right corner
    x_topright = crval1 + ((lenx/2) * mod2) + ((leny/2) * mod1) 
    y_topright = crval2 - ((lenx/2) * mod1) + ((leny/2) * mod2)

    #bottom right corner
    x_bottomright = crval1 + ((lenx/2) * mod2) - ((leny/2) * mod1)
    y_bottomright = crval2 - ((lenx/2) * mod1) - ((leny/2) * mod2)

    #top left corner
    x_topleft = crval1 - ((lenx/2) * mod2) + ((leny/2) * mod1)
    y_topleft = crval2 + ((lenx/2) * mod1) + ((leny/2) * mod2)

    #bottom left corner
    x_bottomleft = crval1 - ((lenx/2) * mod2) - ((leny/2) * mod1)
    y_bottomleft = crval2 + ((lenx/2) * mod1) - ((leny/2) * mod2)

    if plot:
        xx = [crval1 - (lenx/2), crval1 + (lenx/2), crval1 + (lenx/2), crval1 - (lenx/2), crval1 - (lenx/2)]
        yy = [crval2 - (leny/2), crval2 - (leny/2), crval2 + (leny/2), crval2 + (leny/2), crval2 - (leny/2)]

        plt.plot(xx,yy,label = 'unrotated')
        plt.plot([x_bottomleft,x_bottomright,x_topright,x_topleft,x_bottomleft], [y_bottomleft,y_bottomright,y_topright,y_topleft,y_bottomleft],
                     label='rotated by ' + str(crota) + ' deg')
        ax = plt.gca()
        ax.axis('equal')
        plt.legend()
        plt.show()
    
    return [x_bottomleft,x_bottomright,x_topright,x_topleft,x_bottomleft],[y_bottomleft,y_bottomright,y_topright,y_topleft,y_bottomleft]


    

    

    
    


