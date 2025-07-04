# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:25:38 2025

@author: sofia
"""

# Climatologia maps 
import os #change folders
import numpy as np # perform calculations and basic math
import matplotlib.pyplot as plt # plot data
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy #work with geographical projections and maps"
import matplotlib.ticker as mticker
from datetime import datetime
import pandas as pd

#%% Wind
# Open the variables

# Transform the time series to months only
time_months = np.empty_like(time_date, dtype=int)
for i in range(0, len(time_months)):
    time_months[i] = np.datetime64(time_date[i], 'M').astype(int) % 12 + 1  # Extract the month as an integer

# Alternatively, you can do this in a more vectorized way:
time_months = (time_date.astype('datetime64[M]').astype(int) % 12) + 1

# Different for the diferetent faixas
u10 = data['u10'] #--> m/s, horizontal
u10 = u10[:, :,  np.logical_or.reduce((time_months == 9, time_months == 10, time_months == 11, time_months == 12, 
                                           time_months == 1, time_months == 2, time_months == 3, time_months == 4))]
v10 = data['v10'] #--> m/s, vertical
v10 = v10[:, :, np.logical_or.reduce((time_months == 9, time_months == 10, time_months == 11, time_months == 12, 
                                           time_months == 1, time_months == 2, time_months == 3, time_months == 4))] 

# Velocity and direction 
# Calcular a velocidade do vento
wind_speed = np.sqrt(u10**2 + v10**2)
wind_speed_mean = np.nanmean(wind_speed, axis=2)
# Calcular a direção do vento em radianos
wind_direction_rad = np.arctan2(v10, u10)
# Converter para graus
# Converter a direção do vento para graus
wind_direction_deg = np.degrees(wind_direction_rad)
# Ajustar a direção para o intervalo [0, 360)
wind_direction_deg = (wind_direction_deg + 360) % 360

wind_direction_deg_mean = np.nanmean(wind_direction_deg, axis=2)

# Map

# Configurar a área de interesse
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection


# Plotar a velocidade do vento
speed_plot = map.pcolormesh(lon, lat, wind_speed_mean, cmap=plt.cm.viridis,vmin=0, vmax=10,
                           transform=ccrs.PlateCarree())
cbar = plt.colorbar(speed_plot)
cbar.set_label('Wind Speed (m/s)', fontsize=24)



# Plotar a direção do vento como vetores
x, y = np.meshgrid(lon, lat)
u = wind_speed_mean * np.cos(np.radians(wind_direction_deg_mean))
v = wind_speed_mean * np.sin(np.radians(wind_direction_deg_mean))

# Slicing the arrays to plot every 2x2 grid points
x_slice = x[::5, ::5]
y_slice = y[::5, ::5]
u_slice = u[::5, ::5]
v_slice = v[::5, ::5]

map.quiver(x_slice, y_slice, u_slice, v_slice, transform=ccrs.PlateCarree(),
           color='black', width=0.003, scale=30, scale_units='inches', pivot='middle')

# Adicionar costas e outros elementos do mapa
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map
map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))

# Adicionar linhas de grade
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
# Título do mapa
plt.title('Wind Speed and Direction (1998-2021)', fontsize=28)

# Mostrar o plot
plt.show()


#%% Currents
# open the varible


# Transforme the time series in months only
time_dates = [pd.Timestamp(date) for date in time_date]  # Ensure all elements are Timestamps
time_months = np.array([date.month for date in time_dates])

# tranform
eastward_sea_water = data['eastward_sea_water_velocity'] #--> h
northward_sea_water = data['northward_sea_water_velocity'] #--> v

eastward_sea_water =  np.squeeze(eastward_sea_water_1, axis=2) #remove deep
eastward_sea_water = eastward_sea_water[:, :,  np.logical_or.reduce((time_months == 9, time_months == 10, time_months == 11, time_months == 12, 
                                           time_months == 1, time_months == 2, time_months == 3, time_months == 4))]


northward_sea_water = np.squeeze(northward_sea_water_1, axis=2)
northward_sea_water = northward_sea_water[:, :, np.logical_or.reduce((time_months == 9, time_months == 10, time_months == 11, time_months == 12, 
                                           time_months == 1, time_months == 2, time_months == 3, time_months == 4))]


# Velocity and direction 
# Calcular a velocidade do vento
seawater_speed = np.sqrt(eastward_sea_water**2 + northward_sea_water**2)

seawater_speed_mean = np.nanmean(seawater_speed, axis=2)

# direction of the wind
seawater_direction_rad = np.arctan2(northward_sea_water, eastward_sea_water)
# Converter para graus
# Converter a direção do vento para graus
seawater_direction_deg = np.degrees(seawater_direction_rad)
# Ajustar a direção para o intervalo [0, 360]
seawater_direction_deg = (seawater_direction_deg + 360) % 360


seawater_direction_deg_mean = np.nanmean(seawater_direction_deg, axis=2)

# Map

# Configurar a área de interesse
plt.figure(figsize=(12,12))
extent = [160, 195, -85, -70]
central_lon = (extent[0] + extent[1]) / 2
central_lat = (extent[2] + extent[3]) / 2
map = plt.axes(projection=cartopy.crs.AzimuthalEquidistant(central_longitude=central_lon, central_latitude=central_lat)) # Define projection


# Plotar a velocidade do vento
speed_plot = map.pcolormesh(lon, lat, seawater_speed_mean, cmap=plt.cm.viridis.reversed(), vmin=0, vmax=0.7,
                           transform=ccrs.PlateCarree())
cbar = plt.colorbar(speed_plot)
cbar.set_label('Sea Water Speed (m/s)', fontsize=24)



# Plotar a direção do vento como vetores
# Plotar a direção da sea water como vetores
x, y = np.meshgrid(lon, lat)
u = seawater_speed_mean * np.cos(np.radians(seawater_direction_deg_mean))
v = seawater_speed_mean * np.sin(np.radians(seawater_direction_deg_mean))

# Slicing the arrays to plot every 2x2 grid points
x_slice = x[::20, ::20]
y_slice = y[::20, ::20]
u_slice = u[::20, ::20]
v_slice = v[::20, ::20]

map.quiver(x_slice, y_slice, u_slice, v_slice, transform=ccrs.PlateCarree(),
           color='black', width=0.004, scale=1, pivot='middle')


# Adicionar costas e outros elementos do mapa
map.coastlines(resolution='10m', color='black', linewidth=1) #plot coastlines
map.set_extent([160, 196, -85, -57]) # set the extent of the map
map.add_feature(cartopy.feature.NaturalEarthFeature(category='physical',
                                                    name='land',
                                                    scale='10m',
                                                    facecolor=cartopy.feature.COLORS['land']))

# Adicionar linhas de grade
gl = map.gridlines(draw_labels=True, alpha=0.5, linestyle='dotted', color='black') # Add gridlines
gl.xlocator= mticker.FixedLocator([140,160,180,-160,-140]) #adjust the longitude
gl.xlabels_top = False # Remove longitude labels from top
gl.ylabels_right = False # Remove latitude labels from right
gl.xlabel_style = {'size': 18, 'color': 'k'} # Change size and colour - X labels
gl.ylabel_style = {'size': 18, 'color': 'k'} # Change size and colour - Y labels
# Título do mapa
plt.title('Sea Water Speed and Direction (1998-2021)', fontsize=28)

# Mostrar o plot
plt.show()

#%% Sea Ice
# open varible

# Mean and %
sea_ice_fractionALL_choosemonth_mean = np.nanmean(sea_ice_fractionALL_choosemonth, axis=2)*100
del(sea_ice_fractionALL_choosemonth)

# If <15% = nan
sea_ice_fractionALL_choosemonth_mean[sea_ice_fractionALL_choosemonth_mean < 15] = np.nan

# Compute min and max values
sst_min = np.nanmin(sea_ice_fractionALL_choosemonth_mean)
sst_max = np.nanmax(sea_ice_fractionALL_choosemonth_mean)


# Map


# Map
# Define custom blue colormap
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

# Plot data with Inferno colormap
f1 = map.pcolormesh(lon, lat, sea_ice_fractionALL_choosemonth_mean, 
                    vmin=0, vmax=100, cmap=plt.cm.Blues,
                    transform=cartopy.crs.PlateCarree())

# Create a ScalarMappable object
scalar_mappable = plt.cm.ScalarMappable(cmap=plt.cm.Blues)
scalar_mappable.set_array(sea_ice_fractionALL_choosemonth_mean) 

# Add colorbar with label
cbar = plt.colorbar(scalar_mappable)
cbar.set_label('Sea-ice concentration (%)', fontsize=24) #add a label to the colorbar



# Adjust color bar font size and tick labels
cbar.ax.tick_params(labelsize=18, labelcolor='k')


plt.title('Sea ice concentration (1998-2021)', fontsize=28)
plt.show()


#%% SST

# Open the varible 

# Mean 
analysed_mean = np.nanmean(analysed_sst, axis=2)
analysed_mean_c = analysed_mean −273.15


# Compute min and max values
sst_min = np.nanmin(analysed_sstALL_choosemonth_mean)
sst_max = np.nanmax(analysed_sstALL_choosemonth_mean)

# Map

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

# Plot data with Inferno colormap
f1 = map.pcolormesh(lon, lat, analysed_mean_c, 
                    vmin=-2, vmax=5, cmap=plt.cm.inferno,
                    transform=cartopy.crs.PlateCarree())

# Create a ScalarMappable object
scalar_mappable = plt.cm.ScalarMappable(cmap=plt.cm.inferno)
scalar_mappable.set_array(analysed_sstALL_choosemonth_mean) 

# Add colorbar with label
cbar = plt.colorbar(scalar_mappable)
cbar.set_label('Sea Surface Temperature (ºC)', fontsize=24) #add a label to the colorbar



# Adjust color bar font size and tick labels
cbar.ax.tick_params(labelsize=18, labelcolor='k')


plt.title('Sea Surface temperature (1998-2021)', fontsize=28)
plt.show()

plt.close()