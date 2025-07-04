# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:16:40 2025

@author: sofia
"""

#Phenology by pixel now for the Polar regions

import os #change folders,
import numpy as np # perform calculations and basic math
from scipy.interpolate import interp2d #adjust the spatial resolution of the chlorophyll-a concentration data --> 4km to 10km
from datetime import datetime, timedelta
import matplotlib.pyplot as plt # plot data
import matplotlib.ticker as mticker
from matplotlib.patches import Polygon, Path
import pandas as pd # work with dataframes,tables, spreadsheets, etc.,
import netCDF4 as nc4 # work with netcdf files, the standard file for satellite 2D and 3D data,
import cartopy #work with geographical projections and maps"
from scipy import stats #calculate statistics
import matplotlib as mpl
from scipy import integrate
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.colors as mcolors # to created a colourbar
def start_stop(a, trigger_val):
    """"FINDS INDICES OF START AND END OF BLOOMS"""
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False, np.equal(a, trigger_val), False]
    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])
    # Get the start and end indices with slicing along the shifting ones
    return zip(idx[::2], idx[1::2]-1)

#%% Maps

# Open the variables 
# Calculate yearly average for each pixel
# Add the dif lats
# lat 55to60 
chl_yearlycycle1 = np.empty((len(lat1), len(lon), 242))*np.nan # 242 days because we are doing of a polar area 

for i in range(0, len(lat1)):
    print(i)
    for j in range(0, len(lon)):
        # Extract pixel data
        chl_pixeltemp = chl1[i,j, :]
        # Convert to pandas series
        chl_pixel_series = pd.Series(data=chl_pixeltemp, index=time_date)
        chl_pixel_series_daily = chl_pixel_series.resample('D').mean()
        # Calculate average year
        chl_pixel_series_yearlycycle = chl_pixel_series_daily.groupby([chl_pixel_series_daily.index.month, chl_pixel_series_daily.index.day]).mean()
        # Remove February 29th
        chl_pixel_series_yearlycycle = chl_pixel_series_yearlycycle.drop([(2,29)])
        #Bloom cycle
        chl_pixel_series_yearlycycle_sepdec = chl_pixel_series_yearlycycle.loc[9:12]
        chl_pixel_series_yearlycycle_janapr = chl_pixel_series_yearlycycle.loc[1:4]
        chl_pixel_series_yearlycycle = pd.concat([chl_pixel_series_yearlycycle_sepdec, chl_pixel_series_yearlycycle_janapr])
        # add each pixel cycle to 3D matrix
        chl_yearlycycle[i,j,:] = chl_pixel_series_yearlycycle.ravel()
        
 # pre-allocate
 ROI_bloom_mean = np.empty((len(lat), len(lon)))*np.nan
 ROI_bloom_max = np.empty((len(lat), len(lon)))*np.nan
 ROI_bloom_peak = np.empty((len(lat), len(lon)))*np.nan
 ROI_bloom_amplitude = np.empty((len(lat), len(lon)))*np.nan
 ROI_bloom_start = np.empty((len(lat), len(lon)))*np.nan
 ROI_bloom_end = np.empty((len(lat), len(lon)))*np.nan
 ROI_bloom_duration = np.empty((len(lat), len(lon)))*np.nan
 ROI_bloom_area = np.empty((len(lat), len(lon)))*np.nan
 ROI_bloom_num = np.empty((len(lat), len(lon)))*np.nan
 ROI_total_area = np.empty((len(lat), len(lon)))*np.nan

 # Loop by pixel --> daily and weekly data
 for i in range(0, len(lat)):
     print(i)
     for j in range(0, len(lon)):
         # Extract pixel data
         chl_pixeltemp = chl_yearlycycle[i,j,:]
         # Convert to weekly data
         chl_pixel_series = pd.Series(data=chl_pixeltemp,
                                      index=np.arange(datetime(1998,9,1), datetime(1999,5,1), timedelta(days=1)).astype(datetime))
         chl_pixel_series_weekly = chl_pixel_series.resample('8D').mean()
         # Verificar se a série tem pelo menos 10 valores válidos
         if chl_pixel_series_weekly.notna().sum() < 10:
             continue
         # Look at the year (only if you need)
         #plt.plot(chl_pixel_series_weekly)
         #plt.axhline(np.nanmedian(chl_pixel_series_weekly)*1.05)    
         #plt. close 
         # Find date of chlorophyll-a max (peak of the bloom)
         chl_peak_date = chl_pixel_series_weekly.index[np.nanargmax(chl_pixel_series_weekly)]
         # Calculate threshold (+5% above median)
         median_threshold = np.nanmedian(chl_pixel_series_weekly.values)*1.05
         # Find weeks above threshold
         chl_pixel_abovethresh = chl_pixel_series_weekly >= median_threshold
         mask_threshold = chl_pixel_abovethresh.values
         # Find all blooms (periods of weeks above threshold)
         blooms_indices = pd.DataFrame(start_stop(mask_threshold.ravel(), trigger_val=True))  #identify periods above threshold
         # Exclude blooms with less than 2 weeks (15 days is the min. duration used)
         try:  # the blooms_dur has to be put one space in front 
             blooms_dur = blooms_indices.values[:, 1]-blooms_indices.values[:, 0]+1 # For each bloom, calculates duration
         except: # this removes pixels for which max Chl-a do not coincide with main bloom
             continue      
         blooms_indices = blooms_indices.values[np.argwhere(blooms_dur > 2).ravel()] # Exclude blooms with less than 2 weeks 
         # Find main bloom peak
         max_index = chl_pixel_series_weekly.idxmax()
         # this loop goes through each identifie bloom and extracts info from the main bloom
         for l, item in enumerate(blooms_indices):
             if chl_pixel_series_weekly.index[item[0]] <= max_index < chl_pixel_series_weekly.index[item[1]]: # checks if the bloom includes the yearly maximum. Otherwise, it is not the main bloom of the year and the algorithm moves on
                 b_start_ind = item[0]
                 b_start = chl_pixel_series_weekly.index[b_start_ind].weekofyear #calculates bloom initiation day of the year
                 b_end_ind = item[1]
                 b_end = chl_pixel_series_weekly.index[b_end_ind].weekofyear #calculates bloom termination day of the year
                 if b_end < b_start: #checks if bloom day of the year is laters than the start (may occur for late year blooms)
                     b_duration = 52-b_start+b_end
                 else:
                     b_duration = b_end-b_start+1 # calculates bloom duration                  
         bloom_num = len(blooms_indices) # calculate number of blooms in the year
         bloom_max = chl_pixel_series_weekly.max(skipna=True) # maximum chl-a in the year (bloom peak)
         bloom_mean = chl_pixel_series_weekly.mean(skipna=True) # yearly mean chl-a
         bloom_amplitude = bloom_max - bloom_mean # bloom amplitude
         bloom_peak = max_index.weekofyear # week of bloom peak
         try:  
             bloom_area = integrate.simps(chl_pixel_series_weekly[b_start_ind:b_end_ind+1].dropna()) # area of the bloom (measure of chl-a production)
         except: # this removes pixels for which max Chl-a do not coincide with main bloom
             continue            
         total_area = integrate.simps(chl_pixel_series_weekly.dropna()) # total area during the seasonal cycle (measure of chl-a production of the total seasonal cycle)           
         # attribute to each pre-allocated array
         ROI_bloom_mean[i,j] = bloom_mean
         ROI_bloom_max[i,j] = bloom_max
         ROI_bloom_peak[i,j] = bloom_peak
         ROI_bloom_amplitude[i,j] = bloom_amplitude
         ROI_bloom_start[i,j] = b_start
         ROI_bloom_end[i,j] = b_end
         ROI_bloom_duration[i,j] = b_duration
         ROI_bloom_num[i,j] = bloom_num
         ROI_bloom_area[i,j] = bloom_area
         ROI_total_area[i,j] = total_area                                  
         # delete stuff prior to the next cycle
         del(bloom_mean, bloom_max, bloom_peak, bloom_amplitude, b_start,
             b_end, b_duration, bloom_num, bloom_area, total_area,
             b_start_ind, b_end_ind)       
 
        
 
# Create maps for each metric
# Mean
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_mean, cmap=plt.cm.rainbow,
                    norm = mpl.colors.LogNorm(vmin=0.1,vmax=5), transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map
map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels


cbar = plt.colorbar(f1, ticks=[0.1, 0.3, 0.5, 1, 3, 5]) #add a colorbar
cbar.ax.set_yticklabels([0.1, 0.3, 0.5, 1, 3, 5], fontsize=14)
cbar.set_label('Chl $\\it{a}$ (mg.m$^{-3}$)', fontsize=24) #add a label to the colorbar"

plt.title('Mean Chl-a', fontsize=28)

# bathymetry
# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')



plt.close()

# Max
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_max, cmap=plt.cm.rainbow,
                    norm = mpl.colors.LogNorm(vmin=0.1,vmax=7), transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map

map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))
cbar = plt.colorbar(f1, ticks=[0.1, 0.3, 0.5, 1, 3, 5, 7]) #add a colorbar
cbar.ax.set_yticklabels([0.1, 0.3, 0.5, 1, 3, 5,7], fontsize=14)

cbar.set_label('Chl $\\it{a}$ (mg.m$^{-3}$)', fontsize=24) #add a label to the colorbar"
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
#poly_ROI_A = Polygon(list(ROI_A_verts), facecolor=[1,1,1,0], edgecolor='k', linewidth=3, linestyle='-', zorder=2)
#plt.gca().add_patch(poly_ROI_A)
plt.title('Max Chl-a', fontsize=28)

# bathymetry
# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')

plt.close()

# Amplitude
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_amplitude, cmap=plt.cm.rainbow,
                    norm = mpl.colors.LogNorm(vmin=0.1,vmax=10), transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map


map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))

cbar = plt.colorbar(f1, ticks=[0.1, 0.3, 1, 3, 5, 10]) #add a colorbar
cbar.ax.set_yticklabels([0.1, 0.3, 1, 3, 5, 10], fontsize=14)
cbar.set_label('Chl $\\it{a}$ (mg.m$^{-3}$)', fontsize=24) #add a label to the colorbar"
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
#poly_ROI_A = Polygon(list(ROI_A_verts), facecolor=[1,1,1,0], edgecolor='k', linewidth=3, linestyle='-', zorder=2)
#plt.gca().add_patch(poly_ROI_A)
plt.title('Amplitude', fontsize=28)


# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')


plt.close()

# Change to fit the right scale --> weekly and daily 

ROI_bloom_start_change = np.empty((len(lat), len(lon)))*np.nan
for i in range(0, len(lat)):
    print(i) 
    for y in range(0, len(lon)):
        if ROI_bloom_start[i, y] <= 35: # last week of August
            ROI_bloom_start_change[i, y] = ROI_bloom_start[i, y] + 52 # last week of december 
        else:
            ROI_bloom_start_change[i, y] = ROI_bloom_start[i, y] 
        

# Change to fit the rigt scale -->  monthly 
ROI_bloom_start_change = np.empty((len(lat), len(lon)))*np.nan
for i in range(0, len(lat)):
    print(i) 
    for y in range(0, len(lon)):
        if ROI_bloom_start[i, y] <= 9: # setember 
            ROI_bloom_start_change[i, y] = ROI_bloom_start[i, y] + 12 # december 
        else:
            ROI_bloom_start_change[i, y] = ROI_bloom_start[i, y] 


# Create a custom colormap with specific colors for different ranges --> weekly and daily
# Define your custom color palette
colors = [ 'green', 'LightGreen', 'teal', 'skyblue', 'purple', 'Thistle', 'crimson', 'PeachPuff']
# Define the boundaries for each color
bounds = [36, 40, 44, 48, 52, 56, 60, 64]
# Create a ListedColormap without normalization
cmap = mcolors.ListedColormap(colors)


# Start
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_start_change, cmap=cmap,
                    vmin=36,vmax=68, transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map

map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))
cbar = plt.colorbar(f1) #add a colorbar


#cbar.ax.set_yticks([9, 10, 11, 12, 1+12, 2+12, 3+12, 4+12]) # monthly
cbar.ax.set_yticks([36, 40, 44, 48, 52, 56, 60, 64, 68]) # weekly and daily
cbar.ax.set_yticklabels(['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May']) 
#cbar.ax.set_yticklabels([0.1, 0.3, 1, 3, 5], fontsize=14)
cbar.set_label('Week', fontsize=24) #add a label to the colorbar"
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
#poly_ROI_A = Polygon(list(ROI_A_verts), facecolor=[1,1,1,0], edgecolor='k', linewidth=3, linestyle='-', zorder=2)
#plt.gca().add_patch(poly_ROI_A)
plt.title('Bloom Start', fontsize=28)

# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')

plt.close()

# Change to fit the right scale
ROI_bloom_peak_change = np.empty((len(lat), len(lon)))*np.nan
for i in range(0, len(lat)):
    print(i) 
    for y in range(0, len(lon)):
        if ROI_bloom_peak[i, y] <= 35: # last week of August
            ROI_bloom_peak_change[i, y] = ROI_bloom_peak[i, y] + 52 # last week of december 
        else:
            ROI_bloom_peak_change[i, y] = ROI_bloom_peak[i, y] 

# Change to fit the rigt scale -->  monthly 

ROI_bloom_peak_change = np.empty((len(lat), len(lon)))*np.nan
for i in range(0, len(lat)):
    print(i) 
    for y in range(0, len(lon)):
        if ROI_bloom_peak[i, y] <= 9: # setember 
            ROI_bloom_peak_change[i, y] = ROI_bloom_peak[i, y] + 12 # december 
        else:
            ROI_bloom_peak_change[i, y] = ROI_bloom_peak[i, y] 


# Peak
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_peak_change, cmap=cmap,
                    vmin=36,vmax=68, transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map

map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))
cbar = plt.colorbar(f1) #add a colorbar


cbar.ax.set_yticks([36, 40, 44, 48, 52, 56, 60, 64, 68]) # weekly and daily
#cbar.ax.set_yticks([9, 10, 11, 12, 1+12, 2+12, 3+12, 4+12]) # monthly
cbar.ax.set_yticklabels(['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May'])
#cbar.ax.set_yticklabels([0.1, 0.3, 1, 3, 5], fontsize=14)
cbar.set_label('Week', fontsize=24) #add a label to the colorbar"
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
#poly_ROI_A = Polygon(list(ROI_A_verts), facecolor=[1,1,1,0], edgecolor='k', linewidth=3, linestyle='-', zorder=2)
#plt.gca().add_patch(poly_ROI_A)
plt.title('Bloom Peak', fontsize=28)

# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')

plt.close()

# Change to fit the right scale
ROI_bloom_end_change = np.empty((len(lat), len(lon)))*np.nan
for i in range(0, len(lat)):
    print(i) 
    for y in range(0, len(lon)):
        if ROI_bloom_end[i, y] <= 35: # last week of August
            ROI_bloom_end_change[i, y] = ROI_bloom_end[i, y] + 52 # last week of december 
        else:
            ROI_bloom_end_change[i, y] = ROI_bloom_end[i, y] 
            
# Change to fit the rigt scale -->  monthly 

ROI_bloom_end_change = np.empty((len(lat), len(lon)))*np.nan
for i in range(0, len(lat)):
    print(i) 
    for y in range(0, len(lon)):
        if ROI_bloom_end[i, y] <= 9: # setember 
            ROI_bloom_end_change[i, y] = ROI_bloom_end[i, y] + 12 # december 
        else:
            ROI_bloom_end_change[i, y] = ROI_bloom_end[i, y] 

# End
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_end_change, cmap=cmap,
                    vmin=36,vmax=68, transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map

map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))
cbar = plt.colorbar(f1) #add a colorbar

cbar.ax.set_yticks([36, 40, 44, 48, 52, 56, 60, 64, 68]) # monthly
#cbar.ax.set_yticks([37, 41, 45, 49, 1+52, 6+52, 10+52, 15+52]) # weekly and daily
cbar.ax.set_yticklabels(['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May'])
#cbar.ax.set_yticklabels([0.1, 0.3, 1, 3, 5], fontsize=14)
cbar.set_label('Monthly', fontsize=24) #add a label to the colorbar"
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
#poly_ROI_A = Polygon(list(ROI_A_verts), facecolor=[1,1,1,0], edgecolor='k', linewidth=3, linestyle='-', zorder=2)
#plt.gca().add_patch(poly_ROI_A)
plt.title('Bloom End', fontsize=28)

# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')

plt.close()

# Bloom Duration
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_duration, cmap=plt.cm.rainbow, transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map

map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))
cbar = plt.colorbar(f1) #add a colorbar


#cbar.ax.set_yticks([1, 6, 10, 15, 19, 23, 28, 32, 37, 41, 45, 49])
#cbar.ax.set_yticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
#                         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
#cbar.ax.set_yticklabels([0.1, 0.3, 1, 3, 5], fontsize=14)
cbar.set_label('Nº of Weeks', fontsize=24) #add a label to the colorbar"
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
#poly_ROI_A = Polygon(list(ROI_A_verts), facecolor=[1,1,1,0], edgecolor='k', linewidth=3, linestyle='-', zorder=2)
#plt.gca().add_patch(poly_ROI_A)
plt.title('Bloom Duration', fontsize=28)

# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')


plt.close()

# Create a custom colormap with specific colors for different ranges --> weekly and daily
# Define your custom color palette
colors = [ 'teal', 'purple',  'crimson']
# Define the boundaries for each color
bounds = [1, 2, 3,4]
# Create a ListedColormap without normalization
cmap = mcolors.ListedColormap(colors)

# Bloom Number
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_num, cmap=cmap, vmin=1,vmax=4,transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map

map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))

#cbar = plt.colorbar(f1) #add a colorbar
cbar = plt.colorbar(f1, ticks=[1, 2, 3,4]) #add a colorbar
cbar.ax.set_yticklabels([1, 2, 3,4], fontsize=14)
#cbar.ax.set_yticks([1, 6, 10, 15, 19, 23, 28, 32, 37, 41, 45, 49])
#cbar.ax.set_yticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
#                         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
#cbar.ax.set_yticklabels([0.1, 0.3, 1, 3, 5], fontsize=14)
cbar.set_label('Blooms', fontsize=24) #add a label to the colorbar"
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
#poly_ROI_A = Polygon(list(ROI_A_verts), facecolor=[1,1,1,0], edgecolor='k', linewidth=3, linestyle='-', zorder=2)
#plt.gca().add_patch(poly_ROI_A)
plt.title('Bloom Number', fontsize=28)

# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')


plt.close()

# Bloom Area
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
f1 = map.pcolormesh(lon, lat, ROI_bloom_area, cmap=plt.cm.rainbow,
                    norm = mpl.colors.LogNorm(vmin=1,vmax=14),transform=cartopy.crs.PlateCarree())
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map

map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))
#cbar = plt.colorbar(f1)
cbar = plt.colorbar(f1, ticks=[1, 2, 4, 6, 8, 10, 12, 14]) #add a colorbar
cbar.ax.set_yticklabels([1, 2, 4, 6, 8, 10, 12, 14], fontsize=14)
cbar.set_label('Chl $\\it{a}$ (mg.m$^{-3}$)', fontsize=24) #add a label to the colorbar"
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
#poly_ROI_A = Polygon(list(ROI_A_verts), facecolor=[1,1,1,0], edgecolor='k', linewidth=3, linestyle='-', zorder=2)
#plt.gca().add_patch(poly_ROI_A)
plt.title('Bloom Area', fontsize=28)

# Plot 10m
for geom in lines_10m['geometry']:
    if geom.geom_type == 'LineString':  # Check if it's a LineString
        x, y = geom.coords.xy
        map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')
    elif geom.geom_type == 'MultiLineString':  # Check if it's a MultiLineString
        for ls in geom.geoms:  # Iterate over individual LineString objects
            x, y = ls.coords.xy
            map.plot(x, y, c='black', transform=cartopy.crs.PlateCarree(), alpha=0.8, label='-10')

plt.close()


#%% Phenologya by the cluster
import os #change folders,
import numpy as np # perform calculations and basic math,
import matplotlib.pyplot as plt # plot data
import matplotlib.ticker as mticker
from matplotlib.patches import Polygon, Path
import pandas as pd # work with dataframes,tables, spreadsheets, etc.,
import netCDF4 as nc4 # work with netcdf files, the standard file for satellite 2D and 3D data,
import cartopy #work with geographical projections and maps"
from scipy import stats #calculate statistics
import matplotlib as mpl
from scipy import integrate
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.stats import linregress
def start_stop(a, trigger_val):
    """"FINDS INDICES OF START AND END OF BLOOMS"""
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False, np.equal(a, trigger_val), False]
    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])
    # Get the start and end indices with slicing along the shifting ones
    return zip(idx[::2], idx[1::2]-1)

# open the chl-a variable
# Extracting chl-a cycles from the cluster
# Spacial Mean
mean_cluster1_chla = np.nanmean(cluster1_chla_cycles, axis=0)  # Cluster 1 (meio)
mean_cluster2_chla = np.nanmean(cluster2_chla_cycles, axis=0)  # Cluster 2 (norte)
mean_cluster3_chla = np.nanmean(cluster3_chla_cycles, axis=0)  # Cluster 3 (sul)
del(cluster1_chla_cycles, cluster2_chla_cycles, cluster3_chla_cycles)

# Create DataFrame using pandas +  Average for every week (avoid NaNs)
# Integrate (average spatially within the ROI)
# In the case that we want to do the entire track
chl = np.float32(chl)
chl_ROI_1D = np.nanmean(chl, (0,1)) #averages within the ROI (3D matrix to 1D timeseries) chl_ROI_pixels --> when the area was to be define

chl_ROI_1D_df = pd.DataFrame(chl_ROI_1D, index=time_date)
chl_ROI_1D_df_daily = chl_ROI_1D_df.resample('D').mean() 

# after doing the daily values we can do the weekly (values tha we are doing to use)
chl_ROI_1D_df_weekly = chl_ROI_1D_df_daily.resample('8D').mean() 
chl_ROI_1D_weekly = chl_ROI_1D_df_weekly.values.ravel()

# extract years from dataframe index
time_years = pd.to_datetime(chl_ROI_1D_df_weekly.index).year.values
time_ROI_1D_weekly = chl_ROI_1D_df_weekly.index.values.ravel()

# Start calculating indices!
# Create matrices for bloom phenology indices (change the nº in empty accordingly to the difference between years)
# 2022- 1998 = 24
ROI_bloom_mean = np.empty(24)*np.nan
ROI_bloom_max = np.empty(24)*np.nan
ROI_bloom_peak = np.empty(24)*np.nan
ROI_bloom_amplitude = np.empty(24)*np.nan
ROI_bloom_start = np.empty(24)*np.nan
ROI_bloom_end = np.empty(24)*np.nan
ROI_bloom_duration = np.empty(24)*np.nan
ROI_bloom_area = np.empty(24)*np.nan
ROI_bloom_num = np.empty(24)*np.nan
ROI_total_area = np.empty(24)*np.nan

# Loop through each year and calculate phenology indices
for k, year in enumerate(range(1999, 2023)):
    print(k)
    # consider excluding the first and last year of the dataset
    idx_year_chl = np.argwhere(time_years == year).ravel()
    # create pandas
    chl_pixel = pd.Series(chl_ROI_1D_weekly, index=time_ROI_1D_weekly)
    # Separate for the year
    chl_pixel_year = chl_pixel.iloc[idx_year_chl]
    # Look at the year (only if you need)
    plt.plot(chl_pixel_year)
    plt.axhline(np.nanmedian(chl_pixel_year)*1.05)    
    #plt. close
    # Find date of chlorophyll-a max (peak of the bloom)
    chl_peak_date = chl_pixel_year.index[np.nanargmax(chl_pixel_year)]
    # Setup Bloom window in Polar region (Antartica)
    sep_prevyear_arg = np.argwhere([(chl_pixel.index.month == 9) & (chl_pixel.index.year == year-1)])[0][1]
    apr_year_arg = np.argwhere([(chl_pixel.index.month == 4) & (chl_pixel.index.year == year)])[-1][1]+1
    chl_pixel_bloomwindow = chl_pixel[sep_prevyear_arg:apr_year_arg]
    print(chl_pixel_bloomwindow) # check if period is correct
    # Check new window (only if you need)
    #plt.plot(chl_pixel_bloomwindow)
    #plt.axhline(np.nanmedian(chl_pixel_bloomwindow)*1.05)        
    # Calculate threshold (+5% above median)
    median_threshold = np.nanmedian(chl_pixel_bloomwindow.values)*1.05
    # Find weeks above threshold
    chl_pixel_abovethresh = chl_pixel_bloomwindow >= median_threshold
    mask_threshold = chl_pixel_abovethresh.values
    # Find all blooms (periods of weeks above threshold)
    blooms_indices = pd.DataFrame(start_stop(mask_threshold.ravel(), trigger_val=True))  #identify periods above threshold
    # Exclude blooms with less than 2 weeks (15 days is the min. duration used)
    blooms_dur = blooms_indices.values[:, 1]-blooms_indices.values[:, 0]+1 # For each bloom, calculates duration
    blooms_indices = blooms_indices.values[np.argwhere(blooms_dur > 2).ravel()] # Exclude blooms with less than 2 weeks 
    # Find main bloom peak
    max_index = chl_pixel_bloomwindow.idxmax()
    # this loop goes through each identifie bloom and extracts info from the main bloom
    for l, item in enumerate(blooms_indices):
        if chl_pixel_bloomwindow.index[item[0]] <= max_index <= chl_pixel_bloomwindow.index[item[1]]:
            b_start_ind = item[0]
            b_start = chl_pixel_bloomwindow.index[b_start_ind].weekofyear  # calculates bloom initiation day of the year
            b_end_ind = item[1]
            b_end = chl_pixel_bloomwindow.index[b_end_ind].weekofyear  # calculates bloom termination day of the year
            if b_end < b_start:  # checks if bloom day of the year is later than the start (may occur for late year blooms)
                b_duration = 52 - b_start + b_end
            else:
                b_duration = b_end - b_start + 1  # calculates bloom duration
            # Calculate additional bloom-related parameters
            bloom_num = len(blooms_indices)  # calculate the number of blooms in the year
            bloom_max = chl_pixel_bloomwindow.max(skipna=True)  # maximum chl-a in the year (bloom peak)
            bloom_mean = chl_pixel_bloomwindow.mean(skipna=True)  # yearly mean chl-a
            bloom_amplitude = bloom_max - bloom_mean  # bloom amplitude
            bloom_peak = max_index.weekofyear  # week of bloom peak
            bloom_area = integrate.simps(chl_pixel_bloomwindow[b_start_ind:b_end_ind + 1].dropna())  # area of the bloom (measure of chl-a production)
            total_area = integrate.simps(chl_pixel_bloomwindow.dropna())  # total area during the seasonal cycle (measure of chl-a production of the total seasonal cycle)
            # attribute to each pre-allocated array
            ROI_bloom_mean[k] = bloom_mean
            ROI_bloom_max[k] = bloom_max
            ROI_bloom_peak[k] = bloom_peak
            ROI_bloom_amplitude[k] = bloom_amplitude
            ROI_bloom_start[k] = b_start
            ROI_bloom_end[k] = b_end
            ROI_bloom_duration[k] = b_duration
            ROI_bloom_num[k] = bloom_num
            ROI_bloom_area[k] = bloom_area
            ROI_total_area[k] = total_area
    # delete stuff prior to the next cycle
            del(bloom_mean, bloom_max, bloom_peak, bloom_amplitude, b_start, b_end, b_duration, bloom_num, bloom_area, total_area, b_start_ind, b_end_ind)


# Save data in .csv using pandas
roi_phenology = pd.DataFrame([ROI_bloom_mean,ROI_bloom_max,
                                        ROI_bloom_amplitude,ROI_bloom_peak,
                                        ROI_bloom_start, ROI_bloom_end,
                                        ROI_bloom_duration, ROI_bloom_area,
                                        ROI_total_area,
                                        ROI_bloom_num]).transpose()
roi_phenology_columns = ['Mean','Max','Amplitude','Bloom Peak Date','Bloom Start',
                               'Bloom End','Bloom Duration','Bloom Area',
                               'Total Production','Bloom Number']
roi_phenology.columns = roi_phenology_columns