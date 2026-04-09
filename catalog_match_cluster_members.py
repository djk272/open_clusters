#!/usr/bin/env python3

import warnings
import scipy
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.table import Table, QTable
from astroquery.ipac.irsa import Irsa
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
idx_cone_search, sep2d, _ = cmember_coords.match_to_catalog_sky(cone_search_coords)
sep_constraint = sep2d < max_sep
catalog_matches = cone_search_table[idx_cone_search[sep_constraint]]


# 6 Look at the distribution of separations (in degrees) for all of the cross-matched sources
output_path = Path.home() / "open_clusters" / "data" / "out"

#print(str((sep2d_cmembers.deg < 1.0 * u.deg).sum(), len(cone_search_table)))

plt.hist(sep2d[sep_constraint].arcsec, histtype="step", bins=np.logspace(-2, 2.0, 64))
plt.xlabel("separation [arcsec]")
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()

# 7 Save the figure to a file
plt.savefig(output_path / "sep_plot.png", dpi=150)
plt.close()
print("  Saved: sep_plot.png")

# 8 Process the results (results are returned as an Astropy Table)
# write the reuslts to a CSV file

#catalog_matches = cone_search_table[idx_cmembers]
# t = Table()
# t['catalog_matches'] = catalog_matches

catalog_matches.write(str(output_path)+'/output_matches.csv', format='csv', delimiter = ',', overwrite=True)

print(f"Found {len(catalog_matches)} matched objects between the catalogs.")
if len(catalog_matches) > 0:
    print("First 5 results:")
    print(catalog_matches[:5]) # Display selected columns


# With our cross-match done, 
# we can now make Gaia+2MASS color–magnitude diagrams of our candidate NGC 188 members 
# using the information returned by the cross-match:

# Jmag = cmember_table["j_m"][idx_cmembers]  # note that we use the index array returned above
# Gmag = cmember_table["phot_g_mean_mag"]
# Bmag = cmember_table["bp_rp"]

# fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# ax = axes[0]
# ax.scatter(Gmag - Jmag, Gmag, marker="o", color="k", linewidth=0, alpha=0.5)
# ax.set_xlabel("$G - J$")
# ax.set_ylabel("$G$")
# ax.set_xlim(0, 3)
# ax.set_ylim(19, 10)  # backwards because magnitudes!

# ax = axes[1]
# ax.scatter(Bmag - Gmag, Jmag, marker="o", color="k", linewidth=0, alpha=0.5)
# ax.set_xlabel("$G_{BP} - G$")
# ax.set_ylabel("$J$")
# ax.set_xlim(0.2, 1)
# ax.set_ylim(17, 8)  # backwards because magnitudes!

# fig.tight_layout()