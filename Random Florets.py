# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:41:26 2025

@author: sofia
"""

# Random Florest Cluster

import os
import datetime
import matplotlib.ticker as ticker
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import pandas as pd
from matplotlib.patches import Polygon
from matplotlib.path import Path
from tqdm import tqdm
from scipy import integrate
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score
from sklearn.base import clone
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import RFECV
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.inspection import PartialDependenceDisplay
import pandas as pd
import matplotlib.patches as mpatches
import shap
import sklearn
from matplotlib import patches
import math
from netCDF4 import Dataset
import seaborn as sns
from sklearn.inspection import permutation_importance
def serial_date_to_string(srl_no): 
    """Converts serial number time to datetime"""
    new_date = datetime.datetime(1981, 1, 1, 0, 0) + datetime.timedelta(seconds=srl_no)
    return new_date
def permutation_importances(rf, X_train, y_train, metric):
    """Calculates permutation-based variable importances"""
    baseline = metric(rf, X_train, y_train)
    imp = []
    for col in X_train.columns:
        save = X_train[col].copy()
        X_train[col] = np.random.permutation(X_train[col])
        m = metric(rf, X_train, y_train)
        X_train[col] = save
        imp.append(baseline - m)
    return np.array(imp)


# Open variables--> Phenology metrics and diferent abiotic variables

#%% Time: Seprember- Abril

## SST
# Assuming time_date_sst_ice is a numpy array of datetime64 objects
time_date_sstice_year = np.empty(len(time_date_sst_ice), dtype=int)  # Adjusted to create an empty array of the correct shape and type
time_date_sstice_month = np.empty(len(time_date_sst_ice), dtype=int)  # Adjusted to create an empty array of the correct shape and type
for i in range(len(time_date_sst_ice)):
    # Convert numpy.datetime64 to pandas.Timestamp and extract year and month
    timestamp = pd.to_datetime(time_date_sst_ice[i])
    time_date_sstice_year[i] = timestamp.year
    time_date_sstice_month[i] = timestamp.month

del(time_date_sst_ice, timestamp, i)
    
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    sst_sep_temp = mean_cluster_analysed_sst[(time_date_sstice_year == i-1) & (time_date_sstice_month == 9)]
    sst_oct_temp = mean_cluster_analysed_sst[(time_date_sstice_year == i-1) & (time_date_sstice_month == 10)]
    sst_nov_temp = mean_cluster_analysed_sst[(time_date_sstice_year == i-1) & (time_date_sstice_month == 11)]
    sst_dez_temp = mean_cluster_analysed_sst[(time_date_sstice_year == i-1) & (time_date_sstice_month == 12)]
    # Data from January to April of the current year (i)
    sst_jan_temp = mean_cluster_analysed_sst[(time_date_sstice_year == i) & (time_date_sstice_month == 1)]
    sst_feb_temp = mean_cluster_analysed_sst[(time_date_sstice_year == i) & (time_date_sstice_month == 2)]
    sst_mar_temp = mean_cluster_analysed_sst[(time_date_sstice_year == i) & (time_date_sstice_month == 3)]
    sst_apr_temp = mean_cluster_analysed_sst[(time_date_sstice_year == i) & (time_date_sstice_month == 4)]
    # Join all months from September of the previous year to April of the current year
    sst_sep_to_apr_temp = np.hstack((sst_sep_temp, sst_oct_temp, sst_nov_temp, sst_dez_temp, sst_jan_temp, sst_feb_temp, sst_mar_temp, sst_apr_temp))
    # Calculate the mean for each individual month
    sst_sep_temp = np.nanmean(sst_sep_temp)
    sst_oct_temp = np.nanmean(sst_oct_temp)
    sst_nov_temp = np.nanmean(sst_nov_temp)
    sst_dez_temp = np.nanmean(sst_dez_temp)
    sst_jan_temp = np.nanmean(sst_jan_temp)
    sst_feb_temp = np.nanmean(sst_feb_temp)
    sst_mar_temp = np.nanmean(sst_mar_temp)
    sst_apr_temp = np.nanmean(sst_apr_temp)
    # Calculate the mean for the entire period from September to April
    sst_sep_to_apr_temp = np.nanmean(sst_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        sst_sep = sst_sep_temp
        sst_oct = sst_oct_temp
        sst_nov = sst_nov_temp
        sst_dez = sst_dez_temp
        sst_jan = sst_jan_temp
        sst_feb = sst_feb_temp
        sst_mar = sst_mar_temp
        sst_apr = sst_apr_temp
        sst_sep_to_apr = sst_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        sst_sep = np.hstack((sst_sep, sst_sep_temp))
        sst_oct = np.hstack((sst_oct, sst_oct_temp))
        sst_nov = np.hstack((sst_nov, sst_nov_temp))
        sst_dez = np.hstack((sst_dez, sst_dez_temp))
        sst_jan = np.hstack((sst_jan, sst_jan_temp))
        sst_feb = np.hstack((sst_feb, sst_feb_temp))
        sst_mar = np.hstack((sst_mar, sst_mar_temp))
        sst_apr = np.hstack((sst_apr, sst_apr_temp))
        sst_sep_to_apr = np.hstack((sst_sep_to_apr, sst_sep_to_apr_temp))

del(sst_sep, sst_oct, sst_nov, sst_dez, sst_jan, sst_feb, sst_mar, sst_apr, 
    sst_sep_temp, sst_oct_temp, sst_nov_temp, sst_dez_temp, sst_jan_temp, sst_feb_temp, sst_mar_temp, sst_apr_temp, sst_sep_to_apr_temp,
    mean_cluster_analysed_sst)

## Sea Ice

for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    ice_sep_temp = mean_cluster_sea_ice_fraction[(time_date_sstice_year == i-1) & (time_date_sstice_month == 9)]
    ice_oct_temp = mean_cluster_sea_ice_fraction[(time_date_sstice_year == i-1) & (time_date_sstice_month == 10)]
    ice_nov_temp = mean_cluster_sea_ice_fraction[(time_date_sstice_year == i-1) & (time_date_sstice_month == 11)]
    ice_dez_temp = mean_cluster_sea_ice_fraction[(time_date_sstice_year == i-1) & (time_date_sstice_month == 12)]
    # Data from January to April of the current year (i)
    ice_jan_temp = mean_cluster_sea_ice_fraction[(time_date_sstice_year == i) & (time_date_sstice_month == 1)]
    ice_feb_temp = mean_cluster_sea_ice_fraction[(time_date_sstice_year == i) & (time_date_sstice_month == 2)]
    ice_mar_temp = mean_cluster_sea_ice_fraction[(time_date_sstice_year == i) & (time_date_sstice_month == 3)]
    ice_apr_temp = mean_cluster_sea_ice_fraction[(time_date_sstice_year == i) & (time_date_sstice_month == 4)]
    # Join all months from September of the previous year to April of the current year
    ice_sep_to_apr_temp = np.hstack((ice_sep_temp, ice_oct_temp, ice_nov_temp, ice_dez_temp, ice_jan_temp, ice_feb_temp, ice_mar_temp, ice_apr_temp))
    # Calculate the mean for each individual month
    ice_sep_temp = np.nanmean(ice_sep_temp)
    ice_oct_temp = np.nanmean(ice_oct_temp)
    ice_nov_temp = np.nanmean(ice_nov_temp)
    ice_dez_temp = np.nanmean(ice_dez_temp)
    ice_jan_temp = np.nanmean(ice_jan_temp)
    ice_feb_temp = np.nanmean(ice_feb_temp)
    ice_mar_temp = np.nanmean(ice_mar_temp)
    ice_apr_temp = np.nanmean(ice_apr_temp)
    # Calculate the mean for the entire period from September to April
    ice_sep_to_apr_temp = np.nanmean(ice_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        ice_sep = ice_sep_temp
        ice_oct = ice_oct_temp
        ice_nov = ice_nov_temp
        ice_dez = ice_dez_temp
        ice_jan = ice_jan_temp
        ice_feb = ice_feb_temp
        ice_mar = ice_mar_temp
        ice_apr = ice_apr_temp
        ice_sep_to_apr = ice_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        ice_sep = np.hstack((ice_sep, ice_sep_temp))
        ice_oct = np.hstack((ice_oct, ice_oct_temp))
        ice_nov = np.hstack((ice_nov, ice_nov_temp))
        ice_dez = np.hstack((ice_dez, ice_dez_temp))
        ice_jan = np.hstack((ice_jan, ice_jan_temp))
        ice_feb = np.hstack((ice_feb, ice_feb_temp))
        ice_mar = np.hstack((ice_mar, ice_mar_temp))
        ice_apr = np.hstack((ice_apr, ice_apr_temp))
        ice_sep_to_apr = np.hstack((ice_sep_to_apr, ice_sep_to_apr_temp))
        
del(ice_sep, ice_oct, ice_nov, ice_dez, ice_jan, ice_feb, ice_mar, ice_apr, ice_sep_to_apr_temp,
    ice_sep_temp, ice_oct_temp, ice_nov_temp, ice_dez_temp, ice_jan_temp, ice_feb_temp, ice_mar_temp, 
    ice_apr_temp)


## Wind (Velocity and Direction) 
# Assuming time_date_wind is a numpy array of datetime64 objects
time_date_wind_year = np.empty(len(time_date_wind), dtype=int)  # Create an empty array for years
time_date_wind_month = np.empty(len(time_date_wind), dtype=int)  # Create an empty array for months
for i in range(len(time_date_wind)):
    # Convert numpy.datetime64 to pandas.Timestamp and extract year and month
    timestamp = pd.to_datetime(time_date_wind[i])
    time_date_wind_year[i] = timestamp.year
    time_date_wind_month[i] = timestamp.month
    
del(time_date_wind, i, timestamp)

# Velocity   
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    windSpeed_sep_temp = mean_cluster_wind_speed[(time_date_wind_year == i-1) & (time_date_wind_month == 9)]
    windSpeed_oct_temp = mean_cluster_wind_speed[(time_date_wind_year == i-1) & (time_date_wind_month == 10)]
    windSpeed_nov_temp = mean_cluster_wind_speed[(time_date_wind_year == i-1) & (time_date_wind_month == 11)]
    windSpeed_dez_temp = mean_cluster_wind_speed[(time_date_wind_year == i-1) & (time_date_wind_month == 12)]
    # Data from January to April of the current year (i)
    windSpeed_jan_temp = mean_cluster_wind_speed[(time_date_wind_year == i) & (time_date_wind_month == 1)]
    windSpeed_feb_temp = mean_cluster_wind_speed[(time_date_wind_year == i) & (time_date_wind_month == 2)]
    windSpeed_mar_temp = mean_cluster_wind_speed[(time_date_wind_year == i) & (time_date_wind_month == 3)]
    windSpeed_apr_temp = mean_cluster_wind_speed[(time_date_wind_year == i) & (time_date_wind_month == 4)]
    # Join all months from September of the previous year to April of the current year
    windSpeed_sep_to_apr_temp = np.hstack((windSpeed_sep_temp, windSpeed_oct_temp, windSpeed_nov_temp, windSpeed_dez_temp, windSpeed_jan_temp, 
                                           windSpeed_feb_temp, windSpeed_mar_temp, windSpeed_apr_temp))
    # Calculate the mean for each individual month
    windSpeed_sep_temp = np.nanmean(windSpeed_sep_temp)
    windSpeed_oct_temp = np.nanmean(windSpeed_oct_temp)
    windSpeed_nov_temp = np.nanmean(windSpeed_nov_temp)
    windSpeed_dez_temp = np.nanmean(windSpeed_dez_temp)
    windSpeed_jan_temp = np.nanmean(windSpeed_jan_temp)
    windSpeed_feb_temp = np.nanmean(windSpeed_feb_temp)
    windSpeed_mar_temp = np.nanmean(windSpeed_mar_temp)
    windSpeed_apr_temp = np.nanmean(windSpeed_apr_temp)
    # Calculate the mean for the entire period from September to April
    windSpeed_sep_to_apr_temp = np.nanmean(windSpeed_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        windSpeed_sep = windSpeed_sep_temp
        windSpeed_oct = windSpeed_oct_temp
        windSpeed_nov = windSpeed_nov_temp
        windSpeed_dez = windSpeed_dez_temp
        windSpeed_jan = windSpeed_jan_temp
        windSpeed_feb = windSpeed_feb_temp
        windSpeed_mar = windSpeed_mar_temp
        windSpeed_apr = windSpeed_apr_temp
        windSpeed_sep_to_apr = windSpeed_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        windSpeed_sep = np.hstack((windSpeed_sep, windSpeed_sep_temp))
        windSpeed_oct = np.hstack((windSpeed_oct, windSpeed_oct_temp))
        windSpeed_nov = np.hstack((windSpeed_nov, windSpeed_nov_temp))
        windSpeed_dez = np.hstack((windSpeed_dez, windSpeed_dez_temp))
        windSpeed_jan = np.hstack((windSpeed_jan, windSpeed_jan_temp))
        windSpeed_feb = np.hstack((windSpeed_feb, windSpeed_feb_temp))
        windSpeed_mar = np.hstack((windSpeed_mar, windSpeed_mar_temp))
        windSpeed_apr = np.hstack((windSpeed_apr, windSpeed_apr_temp))
        windSpeed_sep_to_apr = np.hstack((windSpeed_sep_to_apr, windSpeed_sep_to_apr_temp))

del(mean_cluster_wind_speed, windSpeed_sep_temp, windSpeed_oct_temp, windSpeed_nov_temp, windSpeed_dez_temp,
    windSpeed_jan_temp, windSpeed_feb_temp, windSpeed_mar_temp, windSpeed_apr_temp, windSpeed_sep_to_apr_temp, 
    windSpeed_sep, windSpeed_oct, windSpeed_nov, windSpeed_dez, windSpeed_jan, windSpeed_feb, windSpeed_mar, windSpeed_apr)

# Direction   
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    windDirection_sep_temp = mean_cluster_wind_direction_deg[(time_date_wind_year == i-1) & (time_date_wind_month == 9)]
    windDirection_oct_temp = mean_cluster_wind_direction_deg[(time_date_wind_year == i-1) & (time_date_wind_month == 10)]
    windDirection_nov_temp = mean_cluster_wind_direction_deg[(time_date_wind_year == i-1) & (time_date_wind_month == 11)]
    windDirection_dez_temp = mean_cluster_wind_direction_deg[(time_date_wind_year == i-1) & (time_date_wind_month == 12)]
    # Data from January to April of the current year (i)
    windDirection_jan_temp = mean_cluster_wind_direction_deg[(time_date_wind_year == i) & (time_date_wind_month == 1)]
    windDirection_feb_temp = mean_cluster_wind_direction_deg[(time_date_wind_year == i) & (time_date_wind_month == 2)]
    windDirection_mar_temp = mean_cluster_wind_direction_deg[(time_date_wind_year == i) & (time_date_wind_month == 3)]
    windDirection_apr_temp = mean_cluster_wind_direction_deg[(time_date_wind_year == i) & (time_date_wind_month == 4)]
    # Join all months from September of the previous year to April of the current year
    windDirection_sep_to_apr_temp = np.hstack((windDirection_sep_temp, windDirection_oct_temp, windDirection_nov_temp, windDirection_dez_temp, windDirection_jan_temp, 
                                           windDirection_feb_temp, windDirection_mar_temp, windDirection_apr_temp))
    # Calculate the mean for each individual month
    windDirection_sep_temp = np.nanmean(windDirection_sep_temp)
    windDirection_oct_temp = np.nanmean(windDirection_oct_temp)
    windDirection_nov_temp = np.nanmean(windDirection_nov_temp)
    windDirection_dez_temp = np.nanmean(windDirection_dez_temp)
    windDirection_jan_temp = np.nanmean(windDirection_jan_temp)
    windDirection_feb_temp = np.nanmean(windDirection_feb_temp)
    windDirection_mar_temp = np.nanmean(windDirection_mar_temp)
    windDirection_apr_temp = np.nanmean(windDirection_apr_temp)
    # Calculate the mean for the entire period from September to April
    windDirection_sep_to_apr_temp = np.nanmean(windDirection_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        windDirection_sep = windDirection_sep_temp
        windDirection_oct = windDirection_oct_temp
        windDirection_nov = windDirection_nov_temp
        windDirection_dez = windDirection_dez_temp
        windDirection_jan = windDirection_jan_temp
        windDirection_feb = windDirection_feb_temp
        windDirection_mar = windDirection_mar_temp
        windDirection_apr = windDirection_apr_temp
        windDirection_sep_to_apr = windDirection_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        windDirection_sep = np.hstack((windDirection_sep, windDirection_sep_temp))
        windDirection_oct = np.hstack((windDirection_oct, windDirection_oct_temp))
        windDirection_nov = np.hstack((windDirection_nov, windDirection_nov_temp))
        windDirection_dez = np.hstack((windDirection_dez, windDirection_dez_temp))
        windDirection_jan = np.hstack((windDirection_jan, windDirection_jan_temp))
        windDirection_feb = np.hstack((windDirection_feb, windDirection_feb_temp))
        windDirection_mar = np.hstack((windDirection_mar, windDirection_mar_temp))
        windDirection_apr = np.hstack((windDirection_apr, windDirection_apr_temp))
        windDirection_sep_to_apr = np.hstack((windDirection_sep_to_apr, windDirection_sep_to_apr_temp))

del(mean_cluster_wind_direction_deg, time_date_wind_year, time_date_wind_month, windDirection_sep_temp, windDirection_oct_temp,
    windDirection_nov_temp, windDirection_dez_temp, windDirection_jan_temp, windDirection_feb_temp, windDirection_mar_temp, 
    windDirection_apr_temp, windDirection_sep_to_apr_temp, windDirection_sep, windDirection_oct, windDirection_nov, 
    windDirection_dez, windDirection_jan, windDirection_feb, windDirection_mar, windDirection_apr)

## Current 0m (Velocity and Direction) 
# Assuming time_date_current_0 is a numpy array of datetime64 objects
time_date_current_year = np.empty(len(time_date_current_0), dtype=int)  # Create an empty array for years
time_date_current_month = np.empty(len(time_date_current_0), dtype=int)  # Create an empty array for months
for i in range(len(time_date_current_0)):  # Loop over the correct length
    # Convert numpy.datetime64 to pandas.Timestamp and extract year and month
    timestamp = pd.to_datetime(time_date_current_0[i])
    time_date_current_year[i] = timestamp.year
    time_date_current_month[i] = timestamp.month
    
del(i, timestamp, time_date_current_0)
    
# Velocity
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    current0Speed_sep_temp = mean_cluster_seawater_speed_0[(time_date_current_year == i-1) & (time_date_current_month == 9)]
    current0Speed_oct_temp = mean_cluster_seawater_speed_0[(time_date_current_year == i-1) & (time_date_current_month == 10)]
    current0Speed_nov_temp = mean_cluster_seawater_speed_0[(time_date_current_year == i-1) & (time_date_current_month == 11)]
    current0Speed_dez_temp = mean_cluster_seawater_speed_0[(time_date_current_year == i-1) & (time_date_current_month == 12)]
    # Data from January to April of the current year (i)
    current0Speed_jan_temp = mean_cluster_seawater_speed_0[(time_date_current_year == i) & (time_date_current_month == 1)]
    current0Speed_feb_temp = mean_cluster_seawater_speed_0[(time_date_current_year == i) & (time_date_current_month == 2)]
    current0Speed_mar_temp = mean_cluster_seawater_speed_0[(time_date_current_year == i) & (time_date_current_month == 3)]
    current0Speed_apr_temp = mean_cluster_seawater_speed_0[(time_date_current_year == i) & (time_date_current_month == 4)]
    # Join all months from September of the previous year to April of the current year
    current0Speed_sep_to_apr_temp = np.hstack((current0Speed_sep_temp, current0Speed_oct_temp, current0Speed_nov_temp, current0Speed_dez_temp, current0Speed_jan_temp, 
                                           current0Speed_feb_temp, current0Speed_mar_temp, current0Speed_apr_temp))
    # Calculate the mean for each individual month
    current0Speed_sep_temp = np.nanmean(current0Speed_sep_temp)
    current0Speed_oct_temp = np.nanmean(current0Speed_oct_temp)
    current0Speed_nov_temp = np.nanmean(current0Speed_nov_temp)
    current0Speed_dez_temp = np.nanmean(current0Speed_dez_temp)
    current0Speed_jan_temp = np.nanmean(current0Speed_jan_temp)
    current0Speed_feb_temp = np.nanmean(current0Speed_feb_temp)
    current0Speed_mar_temp = np.nanmean(current0Speed_mar_temp)
    current0Speed_apr_temp = np.nanmean(current0Speed_apr_temp)
    # Calculate the mean for the entire period from September to April
    current0Speed_sep_to_apr_temp = np.nanmean(current0Speed_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        current0Speed_sep = current0Speed_sep_temp
        current0Speed_oct = current0Speed_oct_temp
        current0Speed_nov = current0Speed_nov_temp
        current0Speed_dez = current0Speed_dez_temp
        current0Speed_jan = current0Speed_jan_temp
        current0Speed_feb = current0Speed_feb_temp
        current0Speed_mar = current0Speed_mar_temp
        current0Speed_apr = current0Speed_apr_temp
        current0Speed_sep_to_apr = current0Speed_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        current0Speed_sep = np.hstack((current0Speed_sep, current0Speed_sep_temp))
        current0Speed_oct = np.hstack((current0Speed_oct, current0Speed_oct_temp))
        current0Speed_nov = np.hstack((current0Speed_nov, current0Speed_nov_temp))
        current0Speed_dez = np.hstack((current0Speed_dez, current0Speed_dez_temp))
        current0Speed_jan = np.hstack((current0Speed_jan, current0Speed_jan_temp))
        current0Speed_feb = np.hstack((current0Speed_feb, current0Speed_feb_temp))
        current0Speed_mar = np.hstack((current0Speed_mar, current0Speed_mar_temp))
        current0Speed_apr = np.hstack((current0Speed_apr, current0Speed_apr_temp))
        current0Speed_sep_to_apr = np.hstack((current0Speed_sep_to_apr, current0Speed_sep_to_apr_temp))

del(mean_cluster_seawater_speed_0, current0Speed_sep_temp, current0Speed_oct_temp, current0Speed_nov_temp, 
    current0Speed_dez_temp, current0Speed_jan_temp, current0Speed_feb_temp, current0Speed_mar_temp, current0Speed_apr_temp, 
    current0Speed_sep_to_apr_temp, current0Speed_sep, current0Speed_oct, current0Speed_nov, current0Speed_dez, current0Speed_jan,
    current0Speed_feb, current0Speed_mar, current0Speed_apr)


# Direction
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    current0Direction_sep_temp = mean_cluster_seawater_direction_deg_0[(time_date_current_year == i-1) & (time_date_current_month == 9)]
    current0Direction_oct_temp = mean_cluster_seawater_direction_deg_0[(time_date_current_year == i-1) & (time_date_current_month == 10)]
    current0Direction_nov_temp = mean_cluster_seawater_direction_deg_0[(time_date_current_year == i-1) & (time_date_current_month == 11)]
    current0Direction_dez_temp = mean_cluster_seawater_direction_deg_0[(time_date_current_year == i-1) & (time_date_current_month == 12)]
    # Data from January to April of the current year (i)
    current0Direction_jan_temp = mean_cluster_seawater_direction_deg_0[(time_date_current_year == i) & (time_date_current_month == 1)]
    current0Direction_feb_temp = mean_cluster_seawater_direction_deg_0[(time_date_current_year == i) & (time_date_current_month == 2)]
    current0Direction_mar_temp = mean_cluster_seawater_direction_deg_0[(time_date_current_year == i) & (time_date_current_month == 3)]
    current0Direction_apr_temp = mean_cluster_seawater_direction_deg_0[(time_date_current_year == i) & (time_date_current_month == 4)]
    # Join all months from September of the previous year to April of the current year
    current0Direction_sep_to_apr_temp = np.hstack((current0Direction_sep_temp, current0Direction_oct_temp, current0Direction_nov_temp, current0Direction_dez_temp, current0Direction_jan_temp, 
                                           current0Direction_feb_temp, current0Direction_mar_temp, current0Direction_apr_temp))
    # Calculate the mean for each individual month
    current0Direction_sep_temp = np.nanmean(current0Direction_sep_temp)
    current0Direction_oct_temp = np.nanmean(current0Direction_oct_temp)
    current0Direction_nov_temp = np.nanmean(current0Direction_nov_temp)
    current0Direction_dez_temp = np.nanmean(current0Direction_dez_temp)
    current0Direction_jan_temp = np.nanmean(current0Direction_jan_temp)
    current0Direction_feb_temp = np.nanmean(current0Direction_feb_temp)
    current0Direction_mar_temp = np.nanmean(current0Direction_mar_temp)
    current0Direction_apr_temp = np.nanmean(current0Direction_apr_temp)
    # Calculate the mean for the entire period from September to April
    current0Direction_sep_to_apr_temp = np.nanmean(current0Direction_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        current0Direction_sep = current0Direction_sep_temp
        current0Direction_oct = current0Direction_oct_temp
        current0Direction_nov = current0Direction_nov_temp
        current0Direction_dez = current0Direction_dez_temp
        current0Direction_jan = current0Direction_jan_temp
        current0Direction_feb = current0Direction_feb_temp
        current0Direction_mar = current0Direction_mar_temp
        current0Direction_apr = current0Direction_apr_temp
        current0Direction_sep_to_apr = current0Direction_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        current0Direction_sep = np.hstack((current0Direction_sep, current0Direction_sep_temp))
        current0Direction_oct = np.hstack((current0Direction_oct, current0Direction_oct_temp))
        current0Direction_nov = np.hstack((current0Direction_nov, current0Direction_nov_temp))
        current0Direction_dez = np.hstack((current0Direction_dez, current0Direction_dez_temp))
        current0Direction_jan = np.hstack((current0Direction_jan, current0Direction_jan_temp))
        current0Direction_feb = np.hstack((current0Direction_feb, current0Direction_feb_temp))
        current0Direction_mar = np.hstack((current0Direction_mar, current0Direction_mar_temp))
        current0Direction_apr = np.hstack((current0Direction_apr, current0Direction_apr_temp))
        current0Direction_sep_to_apr = np.hstack((current0Direction_sep_to_apr, current0Direction_sep_to_apr_temp))
     
del(mean_cluster_seawater_direction_deg_0, current0Direction_sep_temp, current0Direction_oct_temp, current0Direction_nov_temp, 
    current0Direction_dez_temp, current0Direction_jan_temp, current0Direction_feb_temp, current0Direction_mar_temp, current0Direction_apr_temp, 
    current0Direction_sep_to_apr_temp, current0Direction_sep, current0Direction_oct, current0Direction_nov, current0Direction_dez,
    current0Direction_jan, current0Direction_feb, current0Direction_mar, current0Direction_apr)

# Mix layer 
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    oceanMixed0_sep_temp = mean_cluster_ocean_mixed_0[(time_date_current_year == i-1) & (time_date_current_month == 9)]
    oceanMixed0_oct_temp = mean_cluster_ocean_mixed_0[(time_date_current_year == i-1) & (time_date_current_month == 10)]
    oceanMixed0_nov_temp = mean_cluster_ocean_mixed_0[(time_date_current_year == i-1) & (time_date_current_month == 11)]
    oceanMixed0_dez_temp = mean_cluster_ocean_mixed_0[(time_date_current_year == i-1) & (time_date_current_month == 12)]
    # Data from January to April of the current year (i)
    oceanMixed0_jan_temp = mean_cluster_ocean_mixed_0[(time_date_current_year == i) & (time_date_current_month == 1)]
    oceanMixed0_feb_temp = mean_cluster_ocean_mixed_0[(time_date_current_year == i) & (time_date_current_month == 2)]
    oceanMixed0_mar_temp = mean_cluster_ocean_mixed_0[(time_date_current_year == i) & (time_date_current_month == 3)]
    oceanMixed0_apr_temp = mean_cluster_ocean_mixed_0[(time_date_current_year == i) & (time_date_current_month == 4)]
    # Join all months from September of the previous year to April of the current year
    oceanMixed0_sep_to_apr_temp = np.hstack((oceanMixed0_sep_temp, oceanMixed0_oct_temp, oceanMixed0_nov_temp, oceanMixed0_dez_temp, 
                                              oceanMixed0_jan_temp, oceanMixed0_feb_temp, oceanMixed0_mar_temp, oceanMixed0_apr_temp))
    # Calculate the mean for each individual month
    oceanMixed0_sep_temp = np.nanmean(oceanMixed0_sep_temp)
    oceanMixed0_oct_temp = np.nanmean(oceanMixed0_oct_temp)
    oceanMixed0_nov_temp = np.nanmean(oceanMixed0_nov_temp)
    oceanMixed0_dez_temp = np.nanmean(oceanMixed0_dez_temp)
    oceanMixed0_jan_temp = np.nanmean(oceanMixed0_jan_temp)
    oceanMixed0_feb_temp = np.nanmean(oceanMixed0_feb_temp)
    oceanMixed0_mar_temp = np.nanmean(oceanMixed0_mar_temp)
    oceanMixed0_apr_temp = np.nanmean(oceanMixed0_apr_temp)
    # Calculate the mean for the entire period from September to April
    oceanMixed0_sep_to_apr_temp = np.nanmean(oceanMixed0_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        oceanMixed0_sep = oceanMixed0_sep_temp
        oceanMixed0_oct = oceanMixed0_oct_temp
        oceanMixed0_nov = oceanMixed0_nov_temp
        oceanMixed0_dez = oceanMixed0_dez_temp
        oceanMixed0_jan = oceanMixed0_jan_temp
        oceanMixed0_feb = oceanMixed0_feb_temp
        oceanMixed0_mar = oceanMixed0_mar_temp
        oceanMixed0_apr = oceanMixed0_apr_temp
        oceanMixed0_sep_to_apr = oceanMixed0_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        oceanMixed0_sep = np.hstack((oceanMixed0_sep, oceanMixed0_sep_temp))
        oceanMixed0_oct = np.hstack((oceanMixed0_oct, oceanMixed0_oct_temp))
        oceanMixed0_nov = np.hstack((oceanMixed0_nov, oceanMixed0_nov_temp))
        oceanMixed0_dez = np.hstack((oceanMixed0_dez, oceanMixed0_dez_temp))
        oceanMixed0_jan = np.hstack((oceanMixed0_jan, oceanMixed0_jan_temp))
        oceanMixed0_feb = np.hstack((oceanMixed0_feb, oceanMixed0_feb_temp))
        oceanMixed0_mar = np.hstack((oceanMixed0_mar, oceanMixed0_mar_temp))
        oceanMixed0_apr = np.hstack((oceanMixed0_apr, oceanMixed0_apr_temp))
        oceanMixed0_sep_to_apr = np.hstack((oceanMixed0_sep_to_apr, oceanMixed0_sep_to_apr_temp))

del(mean_cluster_ocean_mixed_0, time_date_current_year, time_date_current_month, oceanMixed0_sep_temp, oceanMixed0_oct_temp,
    oceanMixed0_nov_temp, oceanMixed0_dez_temp, oceanMixed0_jan_temp, oceanMixed0_feb_temp, oceanMixed0_mar_temp, 
    oceanMixed0_apr_temp, oceanMixed0_sep_to_apr_temp, oceanMixed0_sep, oceanMixed0_oct, oceanMixed0_nov,
    oceanMixed0_dez, oceanMixed0_jan, oceanMixed0_feb, oceanMixed0_apr)

## Current 50m (Velocity and Direction) 
# Assuming time_date_current is a numpy array of datetime64 objects
time_date_current_year = np.empty(len(time_date_current_50), dtype=int)  # Create an empty array for years
time_date_current_month = np.empty(len(time_date_current_50), dtype=int)  # Create an empty array for months
for i in range(len(time_date_current_50)):  # Use the correct length
    # Convert numpy.datetime64 to pandas.Timestamp and extract year and month
    timestamp = pd.to_datetime(time_date_current_50[i])  # Convert to pandas Timestamp
    time_date_current_year[i] = timestamp.year  # Extract year
    time_date_current_month[i] = timestamp.month  # Extract month

del(time_date_current_50, i,timestamp)

# Velocity
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    current50Speed_sep_temp = mean_cluster_seawater_speed_50[(time_date_current_year == i-1) & (time_date_current_month == 9)]
    current50Speed_oct_temp = mean_cluster_seawater_speed_50[(time_date_current_year == i-1) & (time_date_current_month == 10)]
    current50Speed_nov_temp = mean_cluster_seawater_speed_50[(time_date_current_year == i-1) & (time_date_current_month == 11)]
    current50Speed_dez_temp = mean_cluster_seawater_speed_50[(time_date_current_year == i-1) & (time_date_current_month == 12)]
    # Data from January to April of the current year (i)
    current50Speed_jan_temp = mean_cluster_seawater_speed_50[(time_date_current_year == i) & (time_date_current_month == 1)]
    current50Speed_feb_temp = mean_cluster_seawater_speed_50[(time_date_current_year == i) & (time_date_current_month == 2)]
    current50Speed_mar_temp = mean_cluster_seawater_speed_50[(time_date_current_year == i) & (time_date_current_month == 3)]
    current50Speed_apr_temp = mean_cluster_seawater_speed_50[(time_date_current_year == i) & (time_date_current_month == 4)]
    # Join all months from September of the previous year to April of the current year
    current50Speed_sep_to_apr_temp = np.hstack((current50Speed_sep_temp, current50Speed_oct_temp, current50Speed_nov_temp, current50Speed_dez_temp, current50Speed_jan_temp, 
                                           current50Speed_feb_temp, current50Speed_mar_temp, current50Speed_apr_temp))
    # Calculate the mean for each individual month
    current50Speed_sep_temp = np.nanmean(current50Speed_sep_temp)
    current50Speed_oct_temp = np.nanmean(current50Speed_oct_temp)
    current50Speed_nov_temp = np.nanmean(current50Speed_nov_temp)
    current50Speed_dez_temp = np.nanmean(current50Speed_dez_temp)
    current50Speed_jan_temp = np.nanmean(current50Speed_jan_temp)
    current50Speed_feb_temp = np.nanmean(current50Speed_feb_temp)
    current50Speed_mar_temp = np.nanmean(current50Speed_mar_temp)
    current50Speed_apr_temp = np.nanmean(current50Speed_apr_temp)
    # Calculate the mean for the entire period from September to April
    current50Speed_sep_to_apr_temp = np.nanmean(current50Speed_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        current50Speed_sep = current50Speed_sep_temp
        current50Speed_oct = current50Speed_oct_temp
        current50Speed_nov = current50Speed_nov_temp
        current50Speed_dez = current50Speed_dez_temp
        current50Speed_jan = current50Speed_jan_temp
        current50Speed_feb = current50Speed_feb_temp
        current50Speed_mar = current50Speed_mar_temp
        current50Speed_apr = current50Speed_apr_temp
        current50Speed_sep_to_apr = current50Speed_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        current50Speed_sep = np.hstack((current50Speed_sep, current50Speed_sep_temp))
        current50Speed_oct = np.hstack((current50Speed_oct, current50Speed_oct_temp))
        current50Speed_nov = np.hstack((current50Speed_nov, current50Speed_nov_temp))
        current50Speed_dez = np.hstack((current50Speed_dez, current50Speed_dez_temp))
        current50Speed_jan = np.hstack((current50Speed_jan, current50Speed_jan_temp))
        current50Speed_feb = np.hstack((current50Speed_feb, current50Speed_feb_temp))
        current50Speed_mar = np.hstack((current50Speed_mar, current50Speed_mar_temp))
        current50Speed_apr = np.hstack((current50Speed_apr, current50Speed_apr_temp))
        current50Speed_sep_to_apr = np.hstack((current50Speed_sep_to_apr, current50Speed_sep_to_apr_temp))
        
del(mean_cluster_seawater_speed_50, current50Speed_sep_temp, current50Speed_oct_temp, current50Speed_nov_temp, current50Speed_dez_temp,
    current50Speed_jan_temp, current50Speed_feb_temp, current50Speed_mar_temp, current50Speed_apr_temp, current50Speed_sep_to_apr_temp,
    current50Speed_sep, current50Speed_oct, current50Speed_nov, current50Speed_dez, current50Speed_jan, current50Speed_feb,
    current50Speed_mar, current50Speed_apr)

# Direction
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    current50Direction_sep_temp = mean_cluster_seawater_direction_deg_50[(time_date_current_year == i-1) & (time_date_current_month == 9)]
    current50Direction_oct_temp = mean_cluster_seawater_direction_deg_50[(time_date_current_year == i-1) & (time_date_current_month == 10)]
    current50Direction_nov_temp = mean_cluster_seawater_direction_deg_50[(time_date_current_year == i-1) & (time_date_current_month == 11)]
    current50Direction_dez_temp = mean_cluster_seawater_direction_deg_50[(time_date_current_year == i-1) & (time_date_current_month == 12)]
    # Data from January to April of the current year (i)
    current50Direction_jan_temp = mean_cluster_seawater_direction_deg_50[(time_date_current_year == i) & (time_date_current_month == 1)]
    current50Direction_feb_temp = mean_cluster_seawater_direction_deg_50[(time_date_current_year == i) & (time_date_current_month == 2)]
    current50Direction_mar_temp = mean_cluster_seawater_direction_deg_50[(time_date_current_year == i) & (time_date_current_month == 3)]
    current50Direction_apr_temp = mean_cluster_seawater_direction_deg_50[(time_date_current_year == i) & (time_date_current_month == 4)]
    # Join all months from September of the previous year to April of the current year
    current50Direction_sep_to_apr_temp = np.hstack((current50Direction_sep_temp, current50Direction_oct_temp, current50Direction_nov_temp, current50Direction_dez_temp, current50Direction_jan_temp, 
                                           current50Direction_feb_temp, current50Direction_mar_temp, current50Direction_apr_temp))
    # Calculate the mean for each individual month
    current50Direction_sep_temp = np.nanmean(current50Direction_sep_temp)
    current50Direction_oct_temp = np.nanmean(current50Direction_oct_temp)
    current50Direction_nov_temp = np.nanmean(current50Direction_nov_temp)
    current50Direction_dez_temp = np.nanmean(current50Direction_dez_temp)
    current50Direction_jan_temp = np.nanmean(current50Direction_jan_temp)
    current50Direction_feb_temp = np.nanmean(current50Direction_feb_temp)
    current50Direction_mar_temp = np.nanmean(current50Direction_mar_temp)
    current50Direction_apr_temp = np.nanmean(current50Direction_apr_temp)
    # Calculate the mean for the entire period from September to April
    current50Direction_sep_to_apr_temp = np.nanmean(current50Direction_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        current50Direction_sep = current50Direction_sep_temp
        current50Direction_oct = current50Direction_oct_temp
        current50Direction_nov = current50Direction_nov_temp
        current50Direction_dez = current50Direction_dez_temp
        current50Direction_jan = current50Direction_jan_temp
        current50Direction_feb = current50Direction_feb_temp
        current50Direction_mar = current50Direction_mar_temp
        current50Direction_apr = current50Direction_apr_temp
        current50Direction_sep_to_apr = current50Direction_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        current50Direction_sep = np.hstack((current50Direction_sep, current50Direction_sep_temp))
        current50Direction_oct = np.hstack((current50Direction_oct, current50Direction_oct_temp))
        current50Direction_nov = np.hstack((current50Direction_nov, current50Direction_nov_temp))
        current50Direction_dez = np.hstack((current50Direction_dez, current50Direction_dez_temp))
        current50Direction_jan = np.hstack((current50Direction_jan, current50Direction_jan_temp))
        current50Direction_feb = np.hstack((current50Direction_feb, current50Direction_feb_temp))
        current50Direction_mar = np.hstack((current50Direction_mar, current50Direction_mar_temp))
        current50Direction_apr = np.hstack((current50Direction_apr, current50Direction_apr_temp))
        current50Direction_sep_to_apr = np.hstack((current50Direction_sep_to_apr, current50Direction_sep_to_apr_temp))
        
del(mean_cluster_seawater_direction_deg_50, current50Direction_sep_temp, current50Direction_oct_temp, current50Direction_nov_temp, 
    current50Direction_dez_temp, current50Direction_jan_temp, current50Direction_feb_temp, current50Direction_mar_temp, current50Direction_apr_temp,
    current50Direction_sep_to_apr_temp, current50Direction_sep, current50Direction_oct, current50Direction_nov, current50Direction_dez, 
    current50Direction_jan, current50Direction_feb, current50Direction_mar, current50Direction_apr)
        
## Sea ice and Salinity
# Assuming time_date_sst_ice is a numpy array of datetime64 objects
time_date_salinity_year = np.empty(len(time_date_salinity), dtype=int)  # Adjusted to create an empty array of the correct shape and type
time_date_salinity_month = np.empty(len(time_date_salinity), dtype=int)  # Adjusted to create an empty array of the correct shape and type
for i in range(len(time_date_salinity)):
    # Convert numpy.datetime64 to pandas.Timestamp and extract year and month
    timestamp = pd.to_datetime(time_date_salinity[i])
    time_date_salinity_year[i] = timestamp.year
    time_date_salinity_month[i] = timestamp.month

del(time_date_salinity, timestamp, i)

#  Sea Ice Area Fraction

for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    seaf_sep_temp = mean_cluster_sea_ice_area_fraction[(time_date_salinity_year == i-1) & (time_date_salinity_month == 9)]
    seaf_oct_temp = mean_cluster_sea_ice_area_fraction[(time_date_salinity_year == i-1) & (time_date_salinity_month == 10)]
    seaf_nov_temp = mean_cluster_sea_ice_area_fraction[(time_date_salinity_year == i-1) & (time_date_salinity_month == 11)]
    seaf_dez_temp = mean_cluster_sea_ice_area_fraction[(time_date_salinity_year == i-1) & (time_date_salinity_month == 12)]
    # Data from January to April of the current year (i)
    seaf_jan_temp = mean_cluster_sea_ice_area_fraction[(time_date_salinity_year == i) & (time_date_salinity_month == 1)]
    seaf_feb_temp = mean_cluster_sea_ice_area_fraction[(time_date_salinity_year == i) & (time_date_salinity_month == 2)]
    seaf_mar_temp = mean_cluster_sea_ice_area_fraction[(time_date_salinity_year == i) & (time_date_salinity_month == 3)]
    seaf_apr_temp = mean_cluster_sea_ice_area_fraction[(time_date_salinity_year == i) & (time_date_salinity_month == 4)]
    # Join all months from September of the previous year to April of the current year
    seaf_sep_to_apr_temp = np.hstack((seaf_sep_temp, seaf_oct_temp, seaf_nov_temp, seaf_dez_temp, seaf_jan_temp, seaf_feb_temp, seaf_mar_temp, seaf_apr_temp))
    # Calculate the mean for each individual month
    seaf_sep_temp = np.nanmean(seaf_sep_temp)
    seaf_oct_temp = np.nanmean(seaf_oct_temp)
    seaf_nov_temp = np.nanmean(seaf_nov_temp)
    seaf_dez_temp = np.nanmean(seaf_dez_temp)
    seaf_jan_temp = np.nanmean(seaf_jan_temp)
    seaf_feb_temp = np.nanmean(seaf_feb_temp)
    seaf_mar_temp = np.nanmean(seaf_mar_temp)
    seaf_apr_temp = np.nanmean(seaf_apr_temp)
    # Calculate the mean for the entire period from September to April
    seaf_sep_to_apr_temp = np.nanmean(seaf_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        seaf_sep = seaf_sep_temp
        seaf_oct = seaf_oct_temp
        seaf_nov = seaf_nov_temp
        seaf_dez = seaf_dez_temp
        seaf_jan = seaf_jan_temp
        seaf_feb = seaf_feb_temp
        seaf_mar = seaf_mar_temp
        seaf_apr = seaf_apr_temp
        seaf_sep_to_apr = seaf_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        seaf_sep = np.hstack((seaf_sep, seaf_sep_temp))
        seaf_oct = np.hstack((seaf_oct, seaf_oct_temp))
        seaf_nov = np.hstack((seaf_nov, seaf_nov_temp))
        seaf_dez = np.hstack((seaf_dez, seaf_dez_temp))
        seaf_jan = np.hstack((seaf_jan, seaf_jan_temp))
        seaf_feb = np.hstack((seaf_feb, seaf_feb_temp))
        seaf_mar = np.hstack((seaf_mar, seaf_mar_temp))
        seaf_apr = np.hstack((seaf_apr, seaf_apr_temp))
        seaf_sep_to_apr = np.hstack((seaf_sep_to_apr, seaf_sep_to_apr_temp))

del(seaf_sep, seaf_oct, seaf_nov, seaf_dez, seaf_jan, seaf_feb, seaf_mar, seaf_apr, 
    seaf_sep_temp, seaf_oct_temp, seaf_nov_temp, seaf_dez_temp, seaf_jan_temp, seaf_feb_temp, seaf_mar_temp, seaf_apr_temp, seaf_sep_to_apr_temp,
    mean_cluster_sea_ice_area_fraction)

# Sea Ice thickness
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    seat_sep_temp = mean_cluster_sea_ice_thickness[(time_date_salinity_year == i-1) & (time_date_salinity_month == 9)]
    seat_oct_temp = mean_cluster_sea_ice_thickness[(time_date_salinity_year == i-1) & (time_date_salinity_month == 10)]
    seat_nov_temp = mean_cluster_sea_ice_thickness[(time_date_salinity_year == i-1) & (time_date_salinity_month == 11)]
    seat_dez_temp = mean_cluster_sea_ice_thickness[(time_date_salinity_year == i-1) & (time_date_salinity_month == 12)]
    # Data from January to April of the current year (i)
    seat_jan_temp = mean_cluster_sea_ice_thickness[(time_date_salinity_year == i) & (time_date_salinity_month == 1)]
    seat_feb_temp = mean_cluster_sea_ice_thickness[(time_date_salinity_year == i) & (time_date_salinity_month == 2)]
    seat_mar_temp = mean_cluster_sea_ice_thickness[(time_date_salinity_year == i) & (time_date_salinity_month == 3)]
    seat_apr_temp = mean_cluster_sea_ice_thickness[(time_date_salinity_year == i) & (time_date_salinity_month == 4)]
    # Join all months from September of the previous year to April of the current year
    seat_sep_to_apr_temp = np.hstack((seat_sep_temp, seat_oct_temp, seat_nov_temp, seat_dez_temp, seat_jan_temp, seat_feb_temp, seat_mar_temp, seat_apr_temp))
    # Calculate the mean for each individual month
    seat_sep_temp = np.nanmean(seat_sep_temp)
    seat_oct_temp = np.nanmean(seat_oct_temp)
    seat_nov_temp = np.nanmean(seat_nov_temp)
    seat_dez_temp = np.nanmean(seat_dez_temp)
    seat_jan_temp = np.nanmean(seat_jan_temp)
    seat_feb_temp = np.nanmean(seat_feb_temp)
    seat_mar_temp = np.nanmean(seat_mar_temp)
    seat_apr_temp = np.nanmean(seat_apr_temp)
    # Calculate the mean for the entire period from September to April
    seat_sep_to_apr_temp = np.nanmean(seat_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        seat_sep = seat_sep_temp
        seat_oct = seat_oct_temp
        seat_nov = seat_nov_temp
        seat_dez = seat_dez_temp
        seat_jan = seat_jan_temp
        seat_feb = seat_feb_temp
        seat_mar = seat_mar_temp
        seat_apr = seat_apr_temp
        seat_sep_to_apr = seat_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        seat_sep = np.hstack((seat_sep, seat_sep_temp))
        seat_oct = np.hstack((seat_oct, seat_oct_temp))
        seat_nov = np.hstack((seat_nov, seat_nov_temp))
        seat_dez = np.hstack((seat_dez, seat_dez_temp))
        seat_jan = np.hstack((seat_jan, seat_jan_temp))
        seat_feb = np.hstack((seat_feb, seat_feb_temp))
        seat_mar = np.hstack((seat_mar, seat_mar_temp))
        seat_apr = np.hstack((seat_apr, seat_apr_temp))
        seat_sep_to_apr = np.hstack((seat_sep_to_apr, seat_sep_to_apr_temp))

del(seat_sep, seat_oct, seat_nov, seat_dez, seat_jan, seat_feb, seat_mar, seat_apr, 
    seat_sep_temp, seat_oct_temp, seat_nov_temp, seat_dez_temp, seat_jan_temp, seat_feb_temp, seat_mar_temp, seat_apr_temp, seat_sep_to_apr_temp,
    mean_cluster_sea_ice_thickness)


# Sea Water Salinity
for i in np.arange(1999, 2022):
    # Data from September to December of the previous year (i-1)
    sws_sep_temp = mean_cluster_sea_water_salinity[(time_date_salinity_year == i-1) & (time_date_salinity_month == 9)]
    sws_oct_temp = mean_cluster_sea_water_salinity[(time_date_salinity_year == i-1) & (time_date_salinity_month == 10)]
    sws_nov_temp = mean_cluster_sea_water_salinity[(time_date_salinity_year == i-1) & (time_date_salinity_month == 11)]
    sws_dez_temp = mean_cluster_sea_water_salinity[(time_date_salinity_year == i-1) & (time_date_salinity_month == 12)]
    # Data from January to April of the current year (i)
    sws_jan_temp = mean_cluster_sea_water_salinity[(time_date_salinity_year == i) & (time_date_salinity_month == 1)]
    sws_feb_temp = mean_cluster_sea_water_salinity[(time_date_salinity_year == i) & (time_date_salinity_month == 2)]
    sws_mar_temp = mean_cluster_sea_water_salinity[(time_date_salinity_year == i) & (time_date_salinity_month == 3)]
    sws_apr_temp = mean_cluster_sea_water_salinity[(time_date_salinity_year == i) & (time_date_salinity_month == 4)]
    # Join all months from September of the previous year to April of the current year
    sws_sep_to_apr_temp = np.hstack((sws_sep_temp, sws_oct_temp, sws_nov_temp, sws_dez_temp, sws_jan_temp, sws_feb_temp, sws_mar_temp, sws_apr_temp))
    # Calculate the mean for each individual month
    sws_sep_temp = np.nanmean(sws_sep_temp)
    sws_oct_temp = np.nanmean(sws_oct_temp)
    sws_nov_temp = np.nanmean(sws_nov_temp)
    sws_dez_temp = np.nanmean(sws_dez_temp)
    sws_jan_temp = np.nanmean(sws_jan_temp)
    sws_feb_temp = np.nanmean(sws_feb_temp)
    sws_mar_temp = np.nanmean(sws_mar_temp)
    sws_apr_temp = np.nanmean(sws_apr_temp)
    # Calculate the mean for the entire period from September to April
    sws_sep_to_apr_temp = np.nanmean(sws_sep_to_apr_temp)
    if i == 1999:
        # Initialize variables with data from 1998
        sws_sep = sws_sep_temp
        sws_oct = sws_oct_temp
        sws_nov = sws_nov_temp
        sws_dez = sws_dez_temp
        sws_jan = sws_jan_temp
        sws_feb = sws_feb_temp
        sws_mar = sws_mar_temp
        sws_apr = sws_apr_temp
        sws_sep_to_apr = sws_sep_to_apr_temp
    else:
        # Stack the data for each subsequent year
        sws_sep = np.hstack((sws_sep, sws_sep_temp))
        sws_oct = np.hstack((sws_oct, sws_oct_temp))
        sws_nov = np.hstack((sws_nov, sws_nov_temp))
        sws_dez = np.hstack((sws_dez, sws_dez_temp))
        sws_jan = np.hstack((sws_jan, sws_jan_temp))
        sws_feb = np.hstack((sws_feb, sws_feb_temp))
        sws_mar = np.hstack((sws_mar, sws_mar_temp))
        sws_apr = np.hstack((sws_apr, sws_apr_temp))
        sws_sep_to_apr = np.hstack((sws_sep_to_apr, sws_sep_to_apr_temp))

del(sws_sep, sws_oct, sws_nov, sws_dez, sws_jan, sws_feb, sws_mar, sws_apr, 
    sws_sep_temp, sws_oct_temp, sws_nov_temp, sws_dez_temp, sws_jan_temp, sws_feb_temp, sws_mar_temp, sws_apr_temp, sws_sep_to_apr_temp,
    mean_cluster_sea_water_salinity)


# Creates a DataFrame will all your variables, including the response variable
# In this case, the goal is building a model to explain what drives Phenology metrics
dataframe_all = pd.DataFrame(data=pheno_metrics, index= np.arange(1999, 2022), columns=['Chl_Sep_to_Apr'])
# SST Setembro - April
dataframe_all['SST'] = sst_sep_to_apr
# Sea Ice Setembro - April
dataframe_all['Sea Ice concentration'] = ice_sep_to_apr

# Wind Speed Setembro - April
dataframe_all['Wind Speed'] = windSpeed_sep_to_apr
# Wind Direction Setembro - April
dataframe_all['Wind Direction'] = windDirection_sep_to_apr

# Current Speed 0 Setembro - April
dataframe_all['Current Speed 0'] = current0Speed_sep_to_apr
# Current Direction 0 Setembro - April
dataframe_all['Current Direction 0'] = current0Direction_sep_to_apr
# Mixed Layer 0 Setembro - April
dataframe_all['Mixed layer'] = oceanMixed0_sep_to_apr

# Current Speed 50 Setembro - April
#dataframe_all['Current Speed 50'] = current50Speed_sep_to_apr
# Current Direction 50 Setembro - April
#dataframe_all['Current Direction 50'] = current50Direction_sep_to_apr

# Sea Ice and salinity
# Sea Ice fraction
#dataframe_all['Sea Ice fraction'] = seaf_sep_to_apr
# Sea Ice thickness
#dataframe_all['Sea Ice thickness'] = seat_sep_to_apr
# Sea Water Salinity
dataframe_all['Sea Water Salinity'] = sws_sep_to_apr



# Plots a correlation matrix for the previous dataframe (spearman correlation)
# Suggestion: adjust the figure size to the number of variables in the dataframe
corr = dataframe_all.corr(method='spearman')
fig = plt.figure(figsize=(18,18))
ax = fig.add_subplot(111)
cax = ax.matshow(corr, vmin=-1, vmax=1, cmap=plt.cm.seismic)
fig.colorbar(cax)
ticks = np.arange(0,len(dataframe_all.columns),1)
ax.set_xticks(ticks)
plt.xticks(rotation=90)
ax.set_yticks(ticks)
ax.set_xticklabels(dataframe_all.columns, fontsize=12)
ax.set_yticklabels(dataframe_all.columns, fontsize=12)
plt.tight_layout()

# Save


# Checks which variables are co-correlated (R>0.7) and prints them
# Run Correlations --> Remove the ones that ocorre
correlated_features = set()
correlation_matrix = dataframe_all.corr()
for i in range(len(correlation_matrix .columns)):
    for j in range(i):
        if correlation_matrix.iloc[i, j] >= 0.7: # adapt the threshold if needed
            colname = correlation_matrix.columns[i]
            correlated_features.add(colname)
print(correlated_features)

# Remove the ones in the correlated features
dataframe_all.drop(columns=['Sea Water Salinity', 'Current Direction 0'], inplace=True)

# Separate features (predictors) from labels (response variable)
dataframe_all_original = dataframe_all # Creates a copy of the dataframe
## Remove all rows with gaps 
dataframe_all = dataframe_all.dropna()
## You may also opt to fill the gaps with median/mean or any other value
#dataframe_all = dataframe_all.fillna(dataframe_all.median())
X = dataframe_all.drop(['Chl_Sep_to_Apr'],axis=1) # Remove labels 
features = X.columns # Extract names of the features
y = dataframe_all['Chl_Sep_to_Apr'] #Defining labels
## You may also opt to transform your features
#feature_scaler = StandardScaler()
#X = feature_scaler.fit_transform(dataframe_all.drop(['Bloom Start'],axis=1)) # tentar com e sem, mete as td na mesma caixa


# Fit model with best parameters from above or the ones you want
# Use this to test different combinations manually
rfr = RandomForestRegressor(bootstrap = True,
                            criterion = 'squared_error',
                            max_depth = 70,
                            max_features = 'sqrt',
                            min_samples_leaf = 2,
                            min_samples_split = 4,
                            n_estimators = 500,
                            oob_score = True,
                            random_state = 42)
rfr.fit(X, y) # Fit the model
# Print metrics to compare (MAE, MSE, RMSE, R2 and Oob_score)
y_pred = rfr.predict(X)
print('Mean Absolute Error:', metrics.mean_absolute_error(y, y_pred))
print('Mean Squared Error:', metrics.mean_squared_error(y, y_pred))
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y, y_pred)))
print('R2:', metrics.r2_score(y, y_pred))
print('OOB Score:', rfr.oob_score_) #Performance metric that uses the data subsets that were not used in each tree to validade your model. Very useful to know if your model is overfitting. However, if you have low N, it will always be bad and can be +/- ignored.

#Plot Variable Importance - Random forests do not have significant features per se but it offers a ranking of importance (how important each variable is for the model)
# Use default feature importances (only recommended to have a quick look at the importances) 
#f_i = list(zip(features,rfr.feature_importances_))
#f_i.sort(key = lambda x : x[1])
#plt.barh([x[0] for x in f_i],[x[1] for x in f_i])
# Use permutation importance (recommended! More robust. The more permutations you run the better, but takes more time)
scoring = ['r2']
r_multi = permutation_importance(rfr, X, y,
                           n_repeats=100, #here you change the number of permutation, now = 100
                           random_state=42, scoring=scoring)

# Prints list of variables sorted by importance to the model
# Thanks to the permutations, a mean value of importance and standard deviation for each variable is calculated
# This mean and standard deviation is used to exclude variables with importance that is not significant
# Importance is only significant if the mean importance > 2 X standard deviation (based on a alpha=0.05)
# Create an empty DataFrame to store the results
results = []

for metric in r_multi:
    print(f"{metric}")
    r = r_multi[metric]
    for i in r.importances_mean.argsort()[::-1]:
        if r.importances_mean[i] - 2 * r.importances_std[i] > 0:  # Significance check
            feature_name = features[i]
            importance_mean = r.importances_mean[i]
            importance_std = r.importances_std[i]
            
            # Print to console
            print(f"{feature_name:<8}  {importance_mean:.3f}  +/- {importance_std:.3f}")
            
            # Append to results
            results.append([metric, feature_name, importance_mean, importance_std])



#Plots figure with ranking of variable importances + comparison between original values and predicted values from the model
# Bloom Start and End
fig, axs = plt.subplots(1, 2, figsize=(11,4))
## Plot the actual values
axs[1].plot(dataframe_all_original.index, dataframe_all_original['Bloom End'].values, marker='^', c='#36454F', label = 'Observed BEnd', markersize=10)
# Plot the predicted values
axs[1].scatter(X.index, rfr.predict(X), marker='o', label = 'Predicted BEnd', c='#FF7F50', s=75, edgecolor='k')
axs[1].set_xticks(ticks = np.arange(1999, 2022),
           labels = ['99', '00', '01', '02', '03', '04', '05', '06', '07', '08',
                     '09', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                     '19', '20', '21'], fontsize=9)

axs[1].set_yticks(ticks=[40, 44, 48, 53],
                  labels=['Oct', 'Nov', 'Dez', 'Jan'], fontsize=9)
# Adicionar ticks menores entre os ticks principais no eixo y
axs[1].yaxis.set_minor_locator(ticker.MultipleLocator(1))
# Ajustar o comprimento dos ticks menores
axs[1].tick_params(axis='y', which='minor', length=8) 

axs[1].legend()
axs[1].set_xlabel('Years', fontsize=14)
axs[1].set_ylabel('BEnd', fontsize=14)
# Plot Feature Importance
axs[0].barh(variables_names[::-1], variables_importances[::-1], xerr=variables_stds[::-1], facecolor='#FF7F50', edgecolor='k')
#axs[0].axhline(y=10.5, c='#36454F', alpha=0.5)
#axs[0].axhspan(22.5, 27, facecolor='#36454F', alpha=0.2)
#axs[0].set_ylim(5.5, 27)
#axs[0].axvline(x=np.mean(predictor_importances['Importance']), c='#A50021', linestyle='--')
axs[0].set_xlabel('Variable Importance', fontsize=14)
plt.tight_layout()
os.chdir('C:/Users/sofia/OneDrive/Ambiente de Trabalho/FCUL/2º ano (Dissertação)/ROI/Random Forest/7Test_Semana_corrigida/Coastal/Graphic')
graphs_dir = 'randomforests_allvariables_BEnd.png'
plt.savefig(graphs_dir, format='png', bbox_inches='tight', dpi=300)
plt.close()

# Plots figure with ranking of variable importances + comparison between original values and predicted values from the model
# O resto 
fig, axs = plt.subplots(1, 2, figsize=(11,4))
## Plot the actual values
axs[1].plot(dataframe_all_original.index, dataframe_all_original['Chl_Sep_to_Apr'].values, marker='^', c='#36454F', label = 'Observed Chla', markersize=10)
# Plot the predicted values
axs[1].scatter(X.index, rfr.predict(X), marker='o', label = 'Predicted Chla', c='#009F6B', s=75, edgecolor='k')
axs[1].set_xticks(ticks = np.arange(1999, 2022),
           labels = ['99', '00', '01', '02', '03', '04', '05', '06', '07', '08',
                     '09', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                     '19', '20', '21'], fontsize=9)
axs[1].legend()
axs[1].set_xlabel('Years', fontsize=14)
axs[1].set_ylabel('Chla', fontsize=14)
# Plot Feature Importance
axs[0].barh(variables_names[::-1], variables_importances[::-1], xerr=variables_stds[::-1], facecolor='#009F6B', edgecolor='k')
#axs[0].axhline(y=10.5, c='#36454F', alpha=0.5)
#axs[0].axhspan(22.5, 27, facecolor='#36454F', alpha=0.2)
#axs[0].set_ylim(5.5, 27)
#axs[0].axvline(x=np.mean(predictor_importances['Importance']), c='#A50021', linestyle='--')
axs[0].set_xlabel('Variable Importance', fontsize=14)
plt.tight_layout()
os.chdir('C:/Users/sofia/OneDrive/Ambiente de Trabalho/FCUL/2º ano (Dissertação)/ROI/Random Forest/6Test_no50/Intermediate/Graphic')
graphs_dir = 'randomforests_allvariables_Chla.png'
plt.savefig(graphs_dir, format='png', bbox_inches='tight', dpi=300)
plt.close()


# Plot Partial Dependence Plots for the most important variables (See the modelled response of the response variable to changes in a particular feature)
# Bloom End
fig = plt.figure(figsize=(9,9))
PartialDependenceDisplay.from_estimator(rfr,X, ['Current Speed 0'], kind='average',
                        line_kw={"color": "k", "linewidth": 4})


# Adjust the y-axis tick labels to subtract 52
current_ticks = plt.gca().get_yticks()  # Get current y-axis ticks
plt.gca().set_yticklabels([f"{tick - 52:.1f}" for tick in current_ticks])  # Subtract 52 from each tick
# Increase the font size of the tick labels
plt.tick_params(axis='both', which='major', labelsize=14)  # Adjust the font size of tick labels

# Increase the font size of the axis labels
plt.xlabel('Seawater Current Velocity', fontsize=16)  # Set x-axis label font size
plt.ylabel('Partial Dependence', fontsize=16)  # Set y-axis label font size

plt.tight_layout()
os.chdir('C:/Users/sofia/OneDrive/Ambiente de Trabalho/FCUL/2º ano (Dissertação)/ROI/Random Forest/7Test_Semana_corrigida/Coastal/Graphic/BEnd')
graphs_dir = 'parcialplot_BEnd_Sea Ice concentration.png'
plt.savefig(graphs_dir, format='png', bbox_inches='tight', dpi=300)
plt.close()