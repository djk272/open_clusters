#!/usr/bin/env python3

import warnings
import scipy
import pandas as pd
import astropy.units as u
import numpy as np
from astropy.table import hstack
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.table import Table, QTable
from astroquery.ipac.irsa import Irsa
#from isochrones.mist import MISTIsochroneGrid
#import ezpadova
from astroquery.vizier import Vizier
from pathlib import Path
import glob

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# 1 Path to cone search catalog, and cluster member result files

cone_search_path = Path.home() / "open_clusters" / "data" / "out" / "output_catalog_search.ecsv"
cmember_path = Path.home() / "open_clusters" / "data" / "out" / "NGC2682_M67_result_with_membership.csv"

# print(f"Found {len(catalog_matches)} matched objects between the catalogs.")
# if len(catalog_matches) > 0:
#     print("First 5 results:")
#     print(catalog_matches[:5]) # Display selected columns

# the .read() below produces some warnings that we can safely ignore
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)

# 2 import files into astopy table objects
    cone_search_table = QTable.read(cone_search_path)
    cmember_table = QTable.read(cmember_path)

# 3 create objects fo RA and Dec
    cone_search_ra = QTable.read(cone_search_path)["ra"]#*u.degree
    cone_search_dec = QTable.read(cone_search_path)["dec"]#*u.degree

    cmember_ra = QTable.read(cmember_path)["ra"]#*u.degree
    cmember_dec = QTable.read(cmember_path)["dec"]#*u.degree


# 4 We can now create a single SkyCoord object to represent all of the sources from our cone search results
cone_search_coords = SkyCoord(cone_search_ra, cone_search_dec, unit="deg", frame='icrs')
cmember_coords = SkyCoord(cmember_ra, cmember_dec, unit="deg", frame='icrs')

# 5 matching the catalogs
max_sep = 2 * u.arcsec
idx_cmember, sep2d, _ = cmember_coords.match_to_catalog_sky(cone_search_coords)
sep_constraint = sep2d < max_sep
catalog_matches = cone_search_table[idx_cmember[sep_constraint]]


# 6 Look at the distribution of separations (in degrees) for all of the cross-matched sources
output_path = Path.home() / "open_clusters" / "data" / "out"

plt.hist(sep2d[sep_constraint].arcsec, histtype="step", bins=np.logspace(-2, 2.0, 64))
plt.xlabel("separation [arcsec]")
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()

# 7 Save the figure to a file
plt.savefig(output_path / "sep_plot.png", dpi=150)
plt.close()
print("  Saved: sep_plot.png")



# 8 Combine relevant information from both tables, only matched results
# Adds 'column_to_add' from df2 to df1 where 'id' matches

catalog_matches = hstack([cmember_table, catalog_matches])

# 9 Process the results (results are returned as an Astropy Table)
# write the reuslts to a CSV file

catalog_matches.write(str(output_path)+'/output_matches.csv', format='csv', delimiter = ',', overwrite=True)

print(f"Found {len(catalog_matches)} matched objects between the catalogs.")
if len(catalog_matches) > 0:
    print("First 5 results:")
    print(catalog_matches[:5]) # Display selected columns


# Extracting parallax measurements from Gaia
parallax = catalog_matches["parallax"]  # in milliarcseconds
distance_array = 1 / parallax  # distance array from parallax (kiloparsecs)
ra = catalog_matches["ra_1"]  # right ascension
dec = catalog_matches["dec_1"]  # declination

# Sorting data
sorted_indices = np.argsort(distance_array)
ra_sorted = ra[sorted_indices]
dec_sorted = dec[sorted_indices]
distance_sorted = distance_array[sorted_indices]

# Creating a histogram of the parallax measurements and RA/DEC map
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))

# Histogram
ax[0].hist(parallax, bins=100, color="gray")
ax[0].axvline(
    np.nanmedian(parallax),
    color="red",
    label=f"Median Parallax \n {np.nanmedian(parallax):.3f} milliarcsec",
)
ax[0].set_xlim(-5, 5)
ax[0].set_xticks(np.arange(-5, 5, 1))
ax[0].set_xlabel("Parallax (milliarcseconds)")
ax[0].set_ylabel("Frequency")
ax[0].set_title("Histogram of Cluster Parallax Measurements from Gaia")
ax[0].legend()

# RA/DEC map
plot = ax[1].scatter(
    ra_sorted,
    dec_sorted,
    s=2,
    c=distance_sorted,
    cmap="viridis",
    vmin=0,
    vmax=8,
)
cbar = fig.colorbar(plot, ax=ax[1])
cbar.set_label("Distance (kpc)")
ax[1].set_xlabel("RA (degrees)")
ax[1].set_ylabel("DEC (degrees)")
ax[1].set_title("Position of Parallax Measurements Colored by Distance")

plt.savefig(output_path / "distance_check_plot.png", dpi=150)
plt.close()
print("  Saved: distance_check_plot.png")

# Finding the median parallax of the data
cluster_median = np.nanmedian(parallax)

# Calculating distance using stellar parallax formula
distance = 1.0 / cluster_median

print(f"Median Parallax: {cluster_median} milliarcseconds")
print(f"Distance to Cluster: {distance:.4f} kpc")

# Distance to M13
d = distance * 10**3  # parsecs


# Converting absolute magnitudes from isochrones to apparent
def app_mag(M):
    return M + 5 * np.log10(d / 10)


# Importing isochrone .dat file
parsec_data = pd.read_csv(
    Path.home() / "open_clusters" / "data" / "isochrones.dat", sep=r"\s+", comment="#", header=None
)

# Converting the data file to an array
parsec_data = np.array(parsec_data)


# The isochrone data file is broken up into age bins
log_age = [
    9.30103,
    9.47712,
    9.60206,
    9.69897,
    9.77815,
    9.84510,
    9.90309,
    9.95424,
    10.00000,
    10.04139,
    10.07918,
    10.11394,
    10.14613,
]

# Setting up a dictionary
age_data = {}

# Looping through .txt file to separate by age
for i, age in enumerate(log_age, start=1):
    age_data[f"age{i}"] = parsec_data[parsec_data[:, 2] == age]

# Extracting g and r magnitudes, applying distance modulus formula
# Setting up arrays
r_data = []
g_data = []
color_data2 = []

# Sending g and r magnitudes and color to an array
for i in range(1, len(log_age) + 1):
    data = age_data[f"age{i}"]
    g_mag = app_mag(data[:, 29])
    r_mag = app_mag(data[:, 30])
    color = g_mag - r_mag

    r_data.append(r_mag)
    g_data.append(g_mag)
    color_data2.append(color)

# Plotting the isochrones
color = [
    "brown",
    "pink",
    "magenta",
    "purple",
    "navy",
    "darkturquoise",
    "teal",
    "olive",
    "green",
    "gold",
    "orange",
    "lightcoral",
    "red",
]
fig, ax = plt.subplots(figsize=(10, 8))
for i, j, k, r in zip(color_data2, r_data, log_age, color):
    ax.plot(
        i,
        j,
        alpha=0.5,
        label=f"{(10**k) / 1e9:.0f} billion years old",
        color=r,
    )

# 10 With our cross-match done, 
# we can now make Gaia+2MASS color–magnitude diagrams of our candidate NGC 188 members 
# using the information returned by the cross-match:

# combine the color information columns from 2mass and cluster member table.
# bluer - redder color, for the color. "Mag" is for the magnitute "G" optical, "J" infrared

# Jmag = catalog_matches["j_m"]  # note that we use the index array returned above
Gmag = catalog_matches["phot_g_mean_mag"]
Bmag = catalog_matches["bp_rp"]
# H = catalog_matches["h_m"]
# K = catalog_matches["k_m"]
FUVmag = catalog_matches["fuv_mag"]
NUVmag = catalog_matches["nuv_mag"]

#fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# ax = axes[0]
# ax.set_title("CMD of matched M67 sources in GALEX")
# ax.scatter(FUVmag.value - NUVmag.value, FUVmag.value, marker="o", color="k", linewidth=0, alpha=0.5, s=1)
# ax.set_xlabel("$FUV - NUV$")
# ax.set_ylabel("$FUVmag$")
# ax.set_xlim(-2, 4)
# ax.set_ylim(20, 4)  # backwards because magnitudes!

#ax = axes[1]
ax.set_title("CMD of M67 sources from GAIA with Isochrones")
ax.scatter(Bmag, Gmag, marker="x", color="k", linewidth=0.9, alpha=0.5, s=20)
ax.set_xlabel("B-R")
ax.set_ylabel("$G$")
ax.invert_yaxis()
ax.set_xlim(0, 1.5)
ax.set_ylim(14, 12)  # backwards because magnitudes!
# ax.set_xlim(0, 3)
# ax.set_ylim(20, 10)  # backwards because magnitudes!
ax.legend(fontsize=8)

# Adding annotation for main sequence turnoff
circle = plt.Circle(
    (0.80, 12.50),
    0.25,
    color="red",
    fill=False,
    linewidth=2,
    label="Main Sequence Turnoff",
)
ax.add_patch(circle)
ax.annotate(
    "Main Sequence Turnoff",
    xy=(1.0, 12.50),
    xytext=(1.0 + 0.05, 12.50 - 0.2),
    arrowprops=dict(arrowstyle="->", color="red"),
    fontsize=12,
    color="red",
)


fig.tight_layout()

fig.savefig(output_path / "cm_diagram_figure.png", dpi=150)
print("  Saved: cm_diagram_figure.png")