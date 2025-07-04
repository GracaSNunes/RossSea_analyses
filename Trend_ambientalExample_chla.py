# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 12:08:11 2025

@author: sofia
"""

import os #change folders
import numpy as np # perform calculations and basic math
from scipy.stats import linregress
import matplotlib.pyplot as plt # plot data
import matplotlib.ticker as mticker
import pandas as pd
import cartopy #work with geographical projections and maps
# Change to your working directory where datafiles and scripts are stored

# Open the variable--> salinity as example

# Transforme the time series to focus on the months and years
time_months = np.empty_like(time_date, dtype=int)
for t in range(0, len(time_months)):
    time_months[t] = time_date[t].month
    
time_year = np.empty_like(time_date, dtype=int)
for i in range(0,len(time_year)):
    time_year[i]= time_date[i].year

del(time_date)

#Loop, so we have only the bloom months (sep-apri)
sea_ice_fraction_slope = np.empty((len(lat), len(lon)))**np.nan
for x in np.arange(0, len(lat)):
    print(x)
    for y in np.arange(0, len(lon)):
        for i in np.arange(1999, 2022):
            ## months for year i
            #Sep
            sea_ice_fraction_sep = sea_ice_fraction[x, y, (time_year == i-1) & (time_months == 9)]
            #Oct
            sea_ice_fraction_oct = sea_ice_fraction[x, y, (time_year == i-1) & (time_months == 10)]
            #Nov
            sea_ice_fraction_nov = sea_ice_fraction[x, y, (time_year == i-1) & (time_months ==11)]
            #Dez
            sea_ice_fraction_dez = sea_ice_fraction[x, y, (time_year == i-1) & (time_months == 12)]
            #Jan
            sea_ice_fraction_jan = sea_ice_fraction[x, y, (time_year == i) & (time_months == 1)]
            #Fev
            sea_ice_fraction_fev = sea_ice_fraction[x, y, (time_year == i) & (time_months == 2)]
            #Mar
            sea_ice_fraction_mar = sea_ice_fraction[x, y, (time_year == i) & (time_months == 3)]
            #Apr
            sea_ice_fraction_apr = sea_ice_fraction[x, y, (time_year == i) & (time_months == 4)]
            ## join all months for year i
            sea_ice_fraction_sepapr = np.hstack((sea_ice_fraction_sep, sea_ice_fraction_oct, sea_ice_fraction_nov,
                                             sea_ice_fraction_dez, sea_ice_fraction_jan, sea_ice_fraction_fev, sea_ice_fraction_mar,sea_ice_fraction_apr))
            ## average  
            sea_ice_fraction_sepapr_mean = np.nanmean(sea_ice_fraction_sepapr)
            ## join each year
            if i == 1999:
                sea_ice_fraction_sepapr_allyears = sea_ice_fraction_sepapr_mean
            else:
                sea_ice_fraction_sepapr_allyears = np.hstack((sea_ice_fraction_sepapr_allyears, sea_ice_fraction_sepapr_mean))
        #slope 
        slope, _, _, pvalue, _ = linregress(np.arange(1999, 2022), sea_ice_fraction_sepapr_allyears)
        #if pvalue is less thar 0.05 --> significant (and them all together)
        if pvalue < 0.05:
            sea_ice_fraction_slope[x, y] = slope

# Compute min and max values
sst_min = np.nanmin(sea_ice_fraction_slope)
sst_max = np.nanmax(sea_ice_fraction_slope)


# Map
plt.figure(figsize=(12, 12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat))  # Define projection
map.coastlines(resolution='10m', color='black', linewidth=1)  # Plot coastlines
map.set_extent([160, 196, -85, -57])  # Set the extent of the map
map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black')  # Add gridlines
gl.xlocator = mticker.FixedLocator([140, 160, 180, -160, -140])  # Adjust the longitude
gl.xlabels_top = False  # Remove longitude labels from top
gl.ylabels_right = False  # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'}  # Change size and color - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'}  # Change size and color - Y labels

# Plot data with seismic colormap
f1 = map.pcolormesh(lon, lat, sea_ice_fraction_slope,
                    cmap=plt.cm.seismic, transform=cartopy.crs.PlateCarree(), vmin=-0.02, vmax=0.02)

# Create a ScalarMappable object
scalar_mappable = plt.cm.ScalarMappable(cmap=plt.cm.seismic, norm=plt.Normalize( vmin=-0.02, vmax=0.02))
scalar_mappable.set_array(sea_ice_fraction_slope)

# Add colorbar with label
cbar = plt.colorbar(scalar_mappable)
cbar.set_label('Slope', fontsize=24)  # Add a label to the colorbar

# Adjust color bar font size and tick labels
cbar.ax.tick_params(labelsize=18, labelcolor='k')

plt.title('Salinity Trend (1998-2021)', fontsize=28)

plt.show()

plt.close()

# import os #change folders
import numpy as np # perform calculations and basic math
from scipy.stats import linregress
import matplotlib.pyplot as plt # plot data
import matplotlib.ticker as mticker
import cartopy #work with geographical projections and maps

# Open chl-a

# Transforme the time series to focus on the months and years
time_months = np.empty_like(time_date, dtype=int)
for t in range(0, len(time_months)):
    time_months[t] = time_date[t].month
    
time_year = np.empty_like(time_date, dtype=int)
for i in range(0,len(time_year)):
    time_year[i]= time_date[i].year

del(time_date)

# Loop da trend
chl_slope1 = np.empty((len(lat1), len(lon)))**np.nan
pvalue1 = np.empty((len(lat1), len(lon)))**np.nan
for x in np.arange(0, len(lat1)):
    print(x)
    for y in np.arange(0, len(lon)):
        for i in np.arange(1999, 2023):
            ## months for year i
            #Sep
            chl_sep = chl1[x, y, (time_year == i-1) & (time_months == 9)]
            #Oct
            chl_oct = chl1[x, y, (time_year == i-1) & (time_months == 10)]
            #Nov
            chl_nov = chl1[x, y, (time_year == i-1) & (time_months ==11)]
            #Dez
            chl_dez = chl1[x, y, (time_year == i-1) & (time_months == 12)]
            #Jan
            chl_jan = chl1[x, y, (time_year == i) & (time_months == 1)]
            #Fev
            chl_fev = chl1[x, y, (time_year == i) & (time_months == 2)]
            #Mar
            chl_mar = chl1[x, y, (time_year == i) & (time_months == 3)]
            #Apr
            chl_apr = chl1[x, y, (time_year == i) & (time_months == 4)]
            ## join all months for year i
            chl_sepapr = np.hstack((chl_sep, chl_oct,chl_nov,chl_dez,chl_jan,chl_fev,chl_mar,chl_apr))
            ## average 
            chl_sepapr_mean = np.nanmean(chl_sepapr)
            ## join each year
            if i == 1999:
                chl_sepapr_allyears = chl_sepapr_mean
            else:
                chl_sepapr_allyears = np.hstack((chl_sepapr_allyears, chl_sepapr_mean))
        #slope 
        slope, _, _, pvalue, _ = linregress(np.arange(1999, 2023), chl_sepapr_allyears)
        #if pvalue is less thar 0.05 --> significant (and them all together)
        #if pvalue < 0.05:
        chl_slope1[x, y] = slope
        pvalue1[x,y] = pvalue
        
# Ajustar a Longitude para o contour não ser um problema
lon_adjusted = np.where(lon < 0, lon + 360, lon)
# Cria uma matriz true para os valores que correspondem 
significant_mask= p_value < 0.05 


plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection
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

f1 = map.pcolormesh(lon_adjusted, lat, chl_slope, cmap=plt.cm.seismic,transform=cartopy.crs.PlateCarree(),vmin=-0.10, vmax=0.10) 
cbar = plt.colorbar(f1) #add a colorbar --> Podemos tirar
#cbar.ax.set_yticklabels([0, 50, 70, 150, 200,250], fontsize=14)
cbar.set_label('Slope', fontsize=24) #add a label to the colorbar"



# Delimitando as áreas significativas para o lado leste
# Plotar os contornos
contour = map.contour(lon_adjusted, lat, significant_mask, colors='black', linewidths=0.5, transform=cartopy.crs.PlateCarree())

# Adicionar etiquetas de contorno
clabels = plt.clabel(contour, inline=1, fmt='%1.1f')

# Ajustar as etiquetas de contorno (remover etiqueta específica)
for label in clabels:
    # Exemplo: remover a etiqueta que diz '0.5'
    if label.get_text() == '0.5':
        label.remove()


plt.title('Chl-a Trend (1998-2022)', fontsize=28)