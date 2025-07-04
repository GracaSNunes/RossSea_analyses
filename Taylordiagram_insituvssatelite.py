# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 11:11:08 2025

@author: sofia
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.grid_finder as gf
import mpl_toolkits.axisartist.floating_axes as fa
import os
import pandas as pd


# Taylor diagram

# Load in situ Chla data

# Load satellite Chla data for different products

# Calculate statistics for Taylor diagram
stdref = np.std(chl_insitu)
stds = [np.std(var) / stdref for var in [chl_sat_CCI, chl_sat_OC4SO, chl_sat_GlobColour, chl_sat_Rrs667Linear]]
cors = [np.corrcoef(chl_insitu, var)[0, 1] for var in [chl_sat_CCI, chl_sat_OC4SO, chl_sat_GlobColour, chl_sat_Rrs667Linear]]

# Define graphic
class TaylorDiagram(object):
    def __init__(self, std_dev, fig=None, rect=111, label=''):
        self.std_dev = std_dev
        tr = PolarAxes.PolarTransform()
        rlocs = np.concatenate((np.arange(10) / 10.0, [0.95, 0.99]))
        tlocs = np.arccos(rlocs)
        gl1 = gf.FixedLocator(tlocs)
        tf1 = gf.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Calculate the maximum standard deviation
        max_std = 6.0 * self.std_dev  # Adjust the scaling factor if needed

        # Adjust the range of the radial axis dynamically
        ghelper = fa.GridHelperCurveLinear(tr, extremes=(0, np.pi/2, 0, max_std),
                                           grid_locator1=gl1, tick_formatter1=tf1)
        if fig is None:
            fig = plt.figure()
        ax = fa.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)
        ax.axis["top"].set_axis_direction("bottom")
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation Coefficient")
        ax.axis["left"].set_axis_direction("bottom")
        ax.axis["left"].label.set_text("Normalized Standard Deviation")
        ax.axis["right"].set_axis_direction("top")
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction("left")
        ax.axis["right"].label.set_axis_direction("top")
        ax.axis["bottom"].set_visible(False)
        ax.grid(True)
        self._ax = ax
        self.ax = ax.get_aux_axes(tr)
        self.ax.plot([0], self.std_dev, 'k*', ls='', ms=10, label=label)

    def add_sample(self, std_dev, corr_coef, *args, **kwargs):
        l, = self.ax.plot(np.arccos(corr_coef), std_dev, *args, **kwargs)
        return l

# Create Taylor Diagram
fig = plt.figure()
# Add in situ Chla data point
dia = TaylorDiagram(1, fig=fig, rect=111, label='In situ Chla')

# Add small random noise to standard deviations and correlation coefficients
stds_noisy = [std + np.random.normal(scale=0.01) for std in stds]
cors_noisy = [corr + np.random.normal(scale=0.01) for corr in cors]

# Add satellite Chla data points with noise
for chl_sat, label in zip([chl_sat_CCI, chl_sat_OC4SO, chl_sat_GlobColour, chl_sat_Rrs667Linear], ['CCI', 'OC4-SO', 'GlobColour', 'Rrs667Linear']):
    std_dev_sat = np.std(chl_sat) + np.random.normal(scale=0.01)
    corr_coef_sat = np.corrcoef(chl_insitu, chl_sat)[0, 1] + np.random.normal(scale=0.01)
    dia.add_sample(std_dev_sat, corr_coef_sat, marker='o', label=label)

# Add legend
dia.ax.legend()

# Show plot
plt.show()
