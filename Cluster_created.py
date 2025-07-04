# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:46:52 2025

@author: sofia
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from scipy import stats
import scipy.cluster.hierarchy as hac

## Load each, average for the full period and delete
## Join all chl climatologies
chl_allcycles = np.vstack((chl_7080_yearlycycle, chl_6570_yearlycycle,
                             chl_6065_yearlycycle, chl_5560_yearlycycle))
# Join all latitudes
lat = np.hstack((lat_7080, lat_6570,lat_6065, lat_5560))
lat = np.flip(lat)
chl_allcycles = np.flip(chl_allcycles, axis=0)
#
# Count Ice-free days for each cycle
icefree_days = np.empty((len(lat), len(lon)))*np.nan
for i in range(0, len(lat)):
    for j in range(0, len(lon)):
        icefree_days[i,j] = np.sum(~np.isnan(chl_allcycles[i,j,:]))
icefree_days[icefree_days == 0] = np.nan

# Load and process chlorophyll data --> Load each, average for the full period

# Load and process coefficients --> Load each, average for the full period

# Convert to 1D but keep track of indexes
chl_mean_1D = np.reshape(chl_climatology, -1) # original size is 105 x 310
coeffs_1D = np.reshape(coeffs, -1)
#chl_max_1D = np.reshape(max_chl, -1)
icefree_1D = np.reshape(icefree_days, -1)
# create pandas
clustering_data = pd.DataFrame([chl_mean_1D, coeffs_1D, icefree_1D]).T
clustering_data.columns = ['Chl Mean', 'Coeff', 'Ice Free days']
# drop NaNs and keep indices of NonNans
clustering_data.replace([np.inf, -np.inf], np.nan, inplace=True) # replace infinites with NaNs
clustering_data_withoutNaNs = clustering_data.dropna(axis=0, how='any')
indices_NonNaNs = clustering_data_withoutNaNs.index.values
#%%
#clustering_data_withoutNaNs.isinfinite().values.any()
#np.isinf(clustering_data_withoutNaNs).values.sum() 

#%% Standardize to between -1 and 1
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
data_scaler = StandardScaler()
scaled_data = data_scaler.fit_transform(clustering_data_withoutNaNs)
#%% Run clustering (Euclidean, Ward, 3 clusters)
from sklearn.cluster import AgglomerativeClustering
hierarchical_cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
labels = hierarchical_cluster.fit_predict(scaled_data)
#%% Convert labels back to shape 
# create full list of labels
labels_fulllist = np.empty_like(coeffs_1D, dtype=int) * np.nan
# assign labels from clustering to full list
labels_fulllist[indices_NonNaNs] = labels
# reshape to map dimensions
labels_2D = np.reshape(labels_fulllist, (105, 310)) #this is our map with the clusters

# use the following line to save clusters

# Plot map with clusters
import matplotlib.colors as mcolors # to created a colourbar
# Define your custom color palette
colors = [ 'seagreen', 'skyblue', 'coral']
# Define the boundaries for each color
bounds = [0.00, 0.49, 1.1, 1.25]
labels = ['Int', 'Oce', 'Coa']  # colorbar text
# Create a ListedColormap without normalization
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(bounds, cmap.N)

import matplotlib.ticker as mticker

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

f1 = map.pcolormesh(lon, lat, labels_2D[:-1, :-1], shading='flat', cmap=cmap, norm=norm, transform=ccrs.PlateCarree())


# mid point in the colourbar
mid_points = [(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds) - 1)]

# Color bar
cbar = plt.colorbar(f1, ticks=mid_points, fraction=0.04, pad=0.1)
cbar.ax.set_yticklabels(labels, fontsize=14)  # Define os textos no centro de cada cor
map.text(171.2746, -80.8320, 'Ross Ice Shelf', transform=cartopy.crs.PlateCarree())
plt.show()