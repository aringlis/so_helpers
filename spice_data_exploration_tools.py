# import needed libraries
import sunraster
from sunraster.instr.spice import read_spice_l2_fits
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from astropy.io import fits
import sunpy
import sunpy.map
from astropy.coordinates import SkyCoord
import astropy.units as u
import datetime

def make_spice_summary_plots(raster_file = None, sit_and_stare_file = None, eui_file = None):
    '''Make some summary plots of supplied SPICE raster and sit and stare files.
       These plots are meant to provide a quicklook view of the data in the files.'''

    if raster_file:
        raster_result = read_spice_l2_fits(raster_file)
        n_win = len(raster_result.keys())
        print(' ')
        print(' ')
        print('The loaded raster file contains the following spectral windows:')
        print('----------------------')
        for key in raster_result.keys():
            print(key)
        print('----------------------')
        print(' ')
        print(' ')
    else:
        raise ValueError('Need a raster file to use this tool')

    if sit_and_stare_file:
        stare_result = read_spice_l2_fits(sit_and_stare_file)
        print(' ')
        print(' ')
        print('The loaded sit and stare file contains the following spectral windows:')
        print('----------------------')
        for key in stare_result.keys():
            print(key)
        print('----------------------')
        print(' ')
        print(' ')

        # extract the coordinate info for the sit and stare observation
        sit_win = stare_result[list(stare_result.keys())[0]]
        sit_xx = sit_win.axis_world_coords_values('custom:pos.helioprojective.lon')
        sit_yy = sit_win.axis_world_coords_values('custom:pos.helioprojective.lat')
        sit_x0 = sit_xx.custom_pos_helioprojective_lon.to('arcsec')[0]
        sit_x1 = sit_xx.custom_pos_helioprojective_lon.to('arcsec')[-1]
        sit_y0 = sit_yy.custom_pos_helioprojective_lat.to('arcsec')[0]
        sit_y1 = sit_yy.custom_pos_helioprojective_lat.to('arcsec')[-1]
        sit_tstart = sit_win.meta.date_start.strftime('%H:%M:%S')
        

    
    #Fig 1: show the raster image from each spectral window
    #-------------------------------------------------------

    raster_spectral_windows = raster_result.keys()
    #raster_spectral_windows = ['Ne VIII 770 - Peak','Ly Beta 1025 (Merged)']
    plt_rows = 2
    plt_columns = int(np.ceil(n_win/2))
    #plt_columns = int(len(raster_spectral_windows))
    plt.figure(1,figsize = (15,8))
    plt.subplots_adjust(left=0.07,right=0.95,top=0.95,bottom=0.07,hspace=0.3, wspace = 0.5)

    
    
    for i, r in enumerate(raster_spectral_windows):

        plt.subplot(plt_rows,plt_columns,i+1)
        win = raster_result[r]

        if type(win) != sunraster.spectrogram.SpectrogramCube:
            continue
        else:
            
            #need to sum over the spectral dimension
            win2 = np.nansum(win.data,1)

            yy = win.axis_world_coords_values('custom:pos.helioprojective.lat')
            yy_array = yy.custom_pos_helioprojective_lat.to('arcsec')
            win_cmax = np.nanmax(win2[win2 > 0.0])
            win_cmin = np.nanmin(win2[win2 > 0.0])
            tstart = win.meta.date_start.strftime('%H:%M:%S')
            tend = win.meta.date_end.strftime('%H:%M:%S')
            xx = win.axis_world_coords_values('custom:pos.helioprojective.lon')
            xx_array = xx.custom_pos_helioprojective_lon.to('arcsec')

            # make sure we don't have angle wrapping at 0/360 degrees since this will mess up the plotting
            xx_array[xx_array < -180 *u.deg] += 360 *u.deg
            yy_array[yy_array < -180 *u.deg] += 360 *u.deg

            xx_array2 = xx_array.value
            yy_array2 = yy_array.value

            ax = plt.pcolormesh(xx_array.value,yy_array.value, np.squeeze(win2),shading='nearest',
                                    vmin = win_cmin, vmax = win_cmax, cmap = 'plasma')
            ax.axes.set_aspect('equal') 
            plt.colorbar()
            plt.title(win.meta.spectral_window + ':\n  ' + tstart + ' - ' + tend + ' UT',fontsize=8)

            plt.tick_params(labelsize=8)
            plt.xlabel('x (arcsec)',fontsize=8)
            plt.ylabel('y (arcsec)',fontsize=8)

            #overplot the sit and stare location if file supplied
            if sit_and_stare_file:
                plt.plot([sit_x0, sit_x1],[sit_y0, sit_y1],color='lightgrey',label='Sit and stare: ' + sit_tstart)

    plt.savefig('spice_data_exploration_fig1.png',dpi=300)


    #Fig 2: Show an EUI map (if supplied) with the SPICE FOV on top. Also show the SPICE FOV on a generic Sun from Solar Orbiter view
    #---------------------------------------------
        
    fig2 = plt.figure(2,figsize=(12,6))

    if eui_file:
        eui_map = sunpy.map.Map(eui_file)
    
        
        top_right = SkyCoord(xx_array[-1,-1] + (300 *u.arcsec), yy_array[-1,-1] + (300*u.arcsec), frame=eui_map.coordinate_frame)
        bottom_left = SkyCoord(xx_array[0,0] - (300*u.arcsec), yy_array[0,0] - (300*u.arcsec),frame=eui_map.coordinate_frame)
        eui_submap = eui_map.submap(bottom_left, top_right = top_right)
        ax = fig2.add_subplot(1,2,1,projection = eui_submap)
        plt.subplots_adjust(left=0.07,right=0.95,top=0.95,bottom=0.08,hspace=0.3, wspace = 0.2)
        eui_submap.plot(axes = ax)

        
        #plot the SPICE raster FOV on top
        spice_line1 = SkyCoord([xx_array[0,0],xx_array[-1,0]], [yy_array[0,0], yy_array[-1,0]], frame = eui_submap.coordinate_frame)
        spice_line2 = SkyCoord([xx_array[-1,0],xx_array[-1,-1]], [yy_array[-1,0], yy_array[-1,-1]], frame = eui_submap.coordinate_frame)
        spice_line3 = SkyCoord([xx_array[-1,-1],xx_array[0,-1]], [yy_array[-1,-1], yy_array[0,-1]], frame = eui_submap.coordinate_frame)
        spice_line4 = SkyCoord([xx_array[0,-1],xx_array[0,0]], [yy_array[0,-1], yy_array[0,0]], frame = eui_submap.coordinate_frame)
        
        ax.plot_coord(spice_line1,color='white',label='SPICE Raster FOV')
        ax.plot_coord(spice_line2,color='white')
        ax.plot_coord(spice_line3,color='white')
        ax.plot_coord(spice_line4,color='white')

        if sit_and_stare_file:
            spice_slit = SkyCoord([sit_x0, sit_x1],[sit_y0, sit_y1], frame = eui_submap.coordinate_frame)
            ax.plot_coord(spice_slit, color='lightgrey',label='SPICE stare slit loc.')
  
        plt.legend()

    fig2.add_subplot(1,2,2,projection=None)

    win = raster_result[list(raster_result.keys())[0]]
    rsun = win.meta['RSUN_ARC']

    circle1 = plt.Circle((0, 0), rsun, fill=True, facecolor='orange',alpha=0.5,label='SPICE')
    
    ax = plt.gca()
    ax.add_patch(circle1)

    plt.xlim([-rsun,rsun])
    plt.ylim([-rsun,rsun])
    plt.xlabel('x (arcsec)')
    plt.ylabel('y (arcsec)')

    plt.plot([xx_array2[0,0], xx_array2[0,-1], xx_array2[-1,-1], xx_array2[-1,0], xx_array2[0,0]],
                 [yy_array2[0,0], yy_array2[0,-1], yy_array2[-1,-1], yy_array2[-1,0], yy_array2[0,0]])

    if sit_and_stare_file:
        plt.plot([sit_x0, sit_x1],[sit_y0, sit_y1],color='lightgrey')
        plt.text(0.37,0.95,'Stare Time: ' + sit_win.meta['DATE-BEG'], transform=ax.transAxes)

    plt.axis('equal')
    plt.title('Solar Orbiter Sun view')
    ax = plt.gca()
    plt.text(0.74,0.12,'$D_{sun}$ = ' + str(round(win.meta['DSUN_AU'],2)) + ' AU', transform=ax.transAxes)
    plt.text(0.74,0.07,'HG lon = ' + str(round(win.meta['HGLN_OBS'],1)) + '$^{\circ}$', transform=ax.transAxes)
    plt.text(0.35,0.02,'Raster Time: ' + win.meta['DATE-BEG'], transform=ax.transAxes)
    plt.text(0.74,0.90,'$\Delta t_{corr}$ = ' + str(round(win.meta['EAR_TDEL'],1)) + ' s', transform=ax.transAxes)

    plt.savefig('spice_data_exploration_fig2.png',dpi=200)


    
    #Fig 3: Make a timeseries out of each spectral window in the sit and stare file (if supplied)
    #------------------------------------
        
    if sit_and_stare_file:
              
        fig3 = plt.figure(3,figsize=(12,5))
        plt.subplots_adjust(left=0.1,right=0.75,top=0.95,bottom=0.1,hspace=0.3, wspace = 0.2)
        stare_spectral_windows = stare_result.keys()

        colors = mpl.colormaps['tab20'].colors
        tmax_inds = []
        #plt.subplot(plt_rows,2,5)
        for j, s in enumerate(stare_spectral_windows):
            #extract just a single slit exposure (wavelength and y information at some time t=0)
            line = stare_result[s][:,:,:,0]

            line_timeseries = np.nansum(line.data,(1,2))
            tt = line.axis_world_coords_values('time')

            plt.plot(tt.time,line_timeseries / np.max(line_timeseries),label=s, color = colors[j])
            plt.ylabel('Intensity (norm.)',fontsize=8)
            plt.xlabel('Time (s)',fontsize=8)

            # find the time index with maximum amplitude and store for future use
            tmax_inds.append(np.argmax(line_timeseries))
            

        plt.legend(fontsize=8, bbox_to_anchor =(1.04,1.00),loc='upper left')
        plt.title('Line Intensity Timeseries (Sit-and-stare, summed in y and lambda)')
        plt.savefig('spice_data_exploration_fig3.png',dpi=200)

        
        #Fig 4: Show the spectral line profile for each sit and stare spectral window (if supplied), at the t and y corresponding
        # to maximum intensity.
        #---------------------------
        
        fig4 = plt.figure(4,figsize = (15,8))

        n_sit_win = len(stare_spectral_windows)
        plt_rows = 2
        plt_columns = int(np.ceil(n_sit_win/2))
        plt.subplots_adjust(left=0.07,right=0.95,top=0.95,bottom=0.12,hspace=0.4, wspace = 0.5)
        colors = mpl.colormaps['tab20'].colors

        for l, window in enumerate(stare_spectral_windows):
           # plt.subplot(plt_rows,plt_columns,i+1)
            win = stare_result[window]

            #plot the line profile at its max intensity in t and y
            line = win[tmax_inds[l],:,:,0]
            ymax_ind = np.argmax(np.nansum(line.data,0))
            
            #plot all wavelengths in the window at a max y and max t position
            fig4.add_subplot(plt_rows, plt_columns, l+1,projection = line.wcs,slices=(ymax_ind, 'x'))
            ax = plt.gca()
            ax.plot(line.data[:,ymax_ind], linewidth = 2, drawstyle = 'steps-mid',label=window, color = colors[l])
            plt.title(window)

            plt.text(0.7,0.85,'t = ' + str(tmax_inds[l]), transform=ax.transAxes)
            plt.text(0.7,0.77,'y = ' + str(ymax_ind), transform=ax.transAxes)
            #axis formatting
            ax.coords[2].set_ticklabel(exclude_overlapping=True)
            ax.coords[2].set_format_unit(u.angstrom)
            ax.coords[2].set_major_formatter('x.x')
            ax.coords[2].set_minor_frequency(10)
            ax.coords[2].set_ticklabel(rotation=40, pad=30)
            plt.ylabel('Intensity')

        plt.savefig('spice_data_exploration_fig4.png',dpi=200)
        plt.show()


def make_timeseries_from_sit_and_stare_fileset(sit_and_stare_files):

    # read the first file to find out which spectral windows are in the fileset
    sit_and_stare_files.sort()
    first_file = read_spice_l2_fits(sit_and_stare_files[0])
    spectral_windows = first_file.keys()

    print(' ')
    print(' ')
    print('The loaded sit and stare file contains the following spectral windows:')
    print('----------------------')
    for key in spectral_windows:
        print(key)
    print('----------------------')
    print(' ')
    print(' ')


    lightcurves = {}
    lightcurves['time'] = []
    for s in spectral_windows:
        lightcurves[s] = []

    
    # now open all the files in sequence and append them together
    for f in sit_and_stare_files:
        stare_result = read_spice_l2_fits(f)
    
        # extract the coordinate info for the sit and stare observation
        sit_win = stare_result[list(stare_result.keys())[0]]
        sit_xx = sit_win.axis_world_coords_values('custom:pos.helioprojective.lon')
        sit_yy = sit_win.axis_world_coords_values('custom:pos.helioprojective.lat')
        sit_x0 = sit_xx.custom_pos_helioprojective_lon.to('arcsec')[0]
        sit_x1 = sit_xx.custom_pos_helioprojective_lon.to('arcsec')[-1]
        sit_y0 = sit_yy.custom_pos_helioprojective_lat.to('arcsec')[0]
        sit_y1 = sit_yy.custom_pos_helioprojective_lat.to('arcsec')[-1]
        sit_tstart = sit_win.meta.date_start.datetime #strftime('%H:%M:%S')

        
        
        colors = mpl.colormaps['tab20'].colors
        tmax_inds = []

        for j, s in enumerate(spectral_windows):
            #extract just a single slit exposure (wavelength and y information at some time t=0)
            line = stare_result[s][:,:,:,0]

            line_timeseries = np.nansum(line.data,(1,2))
            tt = line.axis_world_coords_values('time').time
            tt2 = []
            for t in tt:
                tt2.append(sit_tstart + datetime.timedelta(seconds = t.value))
                

            lightcurves[s].extend(line_timeseries)
            if j == 0:
                lightcurves['time'].extend(tt2)
                
            
            

    #now make a plot
    num_rows = len(spectral_windows)
    plt.figure(9, figsize=(14,9))
    plt.subplots_adjust(hspace=0, left = 0.08, right = 0.95, bottom = 0.05, top = 0.95)
    for j, s in enumerate(spectral_windows):
        plt.subplot(9,1,j+1)
        plt.plot(lightcurves['time'], lightcurves[s] / np.max(lightcurves[s]),label=s, color = colors[j])
        plt.ylabel('I (norm.)',fontsize=8)
        
        plt.legend(fontsize=8, loc = 'upper right')
        if j == 0:
            plt.title('Line Intensity Timeseries (Sit-and-stare, summed in y and lambda)')
        if j == num_rows:
            plt.xlabel('Time (UT)',fontsize=8)

    # find the time index with maximum amplitude and store for future use
    tmax_inds.append(np.argmax(line_timeseries))
            

   # plt.legend(fontsize=8, bbox_to_anchor =(1.02,1.00),loc='upper left')
   # plt.title('Line Intensity Timeseries (Sit-and-stare, summed in y and lambda)')
    plt.savefig('spice_sit_and_stare_fileset_timeseries.png',dpi=200)
    plt.show()

    return lightcurves

        
        
    
    
