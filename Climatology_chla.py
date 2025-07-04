# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:12:46 2025

@author: sofia
"""

import os #change folders,
import numpy as np # perform calculations and basic math,
import matplotlib.pyplot as plt # plot data
import matplotlib.ticker as mticker
import cartopy #work with geographical projections and maps"
import matplotlib as mpl
import cartopy.crs as ccrs

# Climatology

# Load the chl-a variavel

## Calculate the mean
chl_mean = np.nanmean(chl,2)

## 1.4 - Map + Chl
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 200, -85, -55]) # set the extent of the map
map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))

# Add names of different land regions --> Adjusted longitude
map.text(171.2746, -80.8320, 'Ross Ice Shelf', transform=cartopy.crs.PlateCarree())
map.text(159.7074, -71.6697, 'Cape Adare', transform=cartopy.crs.PlateCarree())
map.text(143.9247, -74.6761, 'Victoria Land', transform=cartopy.crs.PlateCarree())
map.text(149.3335, -73.3997, 'Terra Nova Bay', transform=cartopy.crs.PlateCarree())
map.text(207.8030, -78.7806, 'Cape Colbeck', transform=cartopy.crs.PlateCarree())  
plt.title('Ross Sea Chla (1998-2022)', fontsize=28)


gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([120,140,160,180,-160,-140,-120]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels


f1 = map.pcolormesh(lon, lat, chl_2001_mean, cmap=plt.cm.rainbow,
                    norm = mpl.colors.LogNorm(vmin=0.1,vmax=5),transform=cartopy.crs.PlateCarree()) # QUESTÃO/CONFIRMAR SE TENHO DEE METER O LOG 10
cbar = plt.colorbar(f1, ticks=[0.1, 0.3, 1, 3, 5]) #add a colorbar
cbar.ax.set_yticklabels([0.1, 0.3, 1, 3, 5], fontsize=14)
cbar.set_label('Chl $\\it{a}$ (mg.m$^{-3}$)', fontsize=24) #add a label to the colorbar"


plt.title('Ross Sea 2001', fontsize=28)
plt.savefig('RossSea.png', format='png', bbox_inches='tight', dpi=300)
# Below, just an example of how to save to a specific directory
# plt.savefig('C:\\Users\\User\\Exercises\\MyFolderForMaps\\mynewchlorophyllmap.png', format='png', bbox_inches='tight', dpi=300)
plt.close()