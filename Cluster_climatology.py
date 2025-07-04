# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:06:24 2025

@author: sofia
"""

# Climatologies
import os #change folders
import geopandas as gpd #shape_files
import numpy as np # perform calculations and basic math
import pandas as pd
from shapely.geometry import Point
import matplotlib.pyplot as plt # plot data
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import cartopy
from matplotlib.lines import Line2D

#%% climatology --> 1D cluster

# open the variavel

# open the cluster

# Extracting chl-a cycles
# cluster do meio (cluster 1)
cluster1_chla_cycles = chl[clusters == 0, :]
# cluster do norte (cluster 2)
cluster2_chla_cycles = chl[clusters == 1, :]
# cluster do sul (cluster 3)
cluster3_chla_cycles = chl[clusters == 2, :]

# Spacial Mean
mean_cluster1_chla = np.nanmean(cluster1_chla_cycles, axis=0)  # Cluster 1 (meio)
mean_cluster2_chla = np.nanmean(cluster2_chla_cycles, axis=0)  # Cluster 2 (norte)
mean_cluster3_chla = np.nanmean(cluster3_chla_cycles, axis=0)  # Cluster 3 (sul)
del(cluster1_chla_cycles, cluster2_chla_cycles, cluster3_chla_cycles)

# Create the time series
ts_chl = pd.Series(mean_cluster1_chla, index=pd.to_datetime(time_date))
# No Feb 29
ts_chl_nofeb29 = ts_chl[~((ts_chl.index.month == 2) & (ts_chl.index.day == 29))]

# Group by month and day --> fica do tamanho de 1 ano, fica com um valor "para cada dia do ano"
chl_mean_grouped = ts_chl_nofeb29.groupby([ts_chl_nofeb29.index.month, ts_chl_nofeb29.index.day]).mean()

# Select data from September to December
chl_mean_sept_to_dec = chl_mean_grouped.loc[9:12]
# Select data from January to April
chl_mean_jan_to_april = chl_mean_grouped.loc[1:3]
# Concatenate the two dataframes
chl_mean_sept_to_april = pd.concat([chl_mean_sept_to_dec, chl_mean_jan_to_april])

# Calculate the moving average after grouping
moving_average = chl_mean_sept_to_april.rolling(window=15, min_periods=15, center=True).mean()
print(moving_average.head(10))  # Print some values after the moving average calculation

# Calculate the 10th and 90th percentiles
percentile_10 = ts_chl_nofeb29.groupby([ts_chl_nofeb29.index.month, ts_chl_nofeb29.index.day]).quantile(0.1)
percentile_90 = ts_chl_nofeb29.groupby([ts_chl_nofeb29.index.month, ts_chl_nofeb29.index.day]).quantile(0.9)

# Adjust Data --> min & max
# Select data from September to December
percentile_10_sept_to_dec = percentile_10.loc[9:12]
percentile_90_sept_to_dec = percentile_90.loc[9:12]
# Select data from January to April
percentile_10_jan_to_april = percentile_10.loc[1:3]
percentile_90_jan_to_april = percentile_90.loc[1:3]
# Concatenate the two dataframes
percentile_10_sept_to_april = pd.concat([percentile_10_sept_to_dec, percentile_10_jan_to_april])
percentile_90_sept_to_april = pd.concat([percentile_90_sept_to_dec, percentile_90_jan_to_april])

# Plotting
plt.plot(moving_average.values, color='green', linestyle='-', linewidth=2.5, label='Chla mean')
plt.fill_between(range(len(moving_average)), percentile_10_sept_to_april.values, percentile_90_sept_to_april.values, color='green', alpha=0.3, label='10th to 90th Percentile')

plt.title('Intermediate (1998-2022)', fontsize=26)
plt.xlim(0, len(moving_average) - 1)  # Set x-axis limits to match the data range
plt.xticks(ticks=[0, 30, 61, 91, 121, 152, 180],  # Update the tick positions to include April
           labels=['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar'], fontsize=18)
plt.xlabel('Months', fontsize=20)
plt.ylim(0,np.nanmax(percentile_90_sept_to_april + 0.05))   # Ensure y-axis starts at 0 (before:  ) # 
plt.yticks(fontsize=18)
plt.ylabel('Chl $\\it{a}$ (mg.m$^{-3}$)', fontsize=20)
plt.yticks(fontsize=18)  # Set the fontsize for y-axis tick labels
plt.legend(fontsize=18)
plt.show()
