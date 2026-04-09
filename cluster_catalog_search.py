#!/usr/bin/env python3

import astropy.units as u
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa import Irsa
from pathlib import Path
import glob
import math

# Path to input files

INPUT_CSV = Path.home() / "open_clusters" / "data" / "out" / "NGC2682_M67_result_with_membership.csv"

# Get a list of all files in the directory
# Get a list of filenames matching a pattern

#cluster_names_list = [file_path.name for file_path in directory_path.glob("*.csv")]

# Load the CSV file into a DataFrame
# the path to one CSV so far to test, open_clusters/data/out/NGC2682_M67_result_with_membership.csv
df = pd.read_csv(INPUT_CSV)

ra_column_name = 'ra'

dec_column_name = 'dec'

# Find the maximum value of ra
ra_max_value = df[ra_column_name].max()

# Find the minimum value of ra
ra_min_value = df[ra_column_name].min()

# Find the maximum value of dec
dec_max_value = df[dec_column_name].max()

# Find the minimum value of dec
dec_min_value = df[dec_column_name].min()


# Find center value for RA

center_ra = (ra_max_value + ra_min_value)/2

# Find center value for Dec

center_dec = (dec_max_value + dec_min_value)/2



# 1. Define the search center, specify coordinates directly
center_coord = SkyCoord(ra=center_ra, dec=center_dec, unit='deg', frame='icrs')
print(f"Searching around coordinates: {center_coord.ra.deg} RA, {center_coord.dec.deg} Dec\n")

# 2. Define the search radius #Find the radius using the the center which of the 2 differences is larger, that will be the radius for the cone
#search. (decmax-center_dec), (ramax-center_ra).
radius = math.sqrt((dec_max_value-center_dec)**2 + (ra_max_value-center_ra)**2) * u.deg # 0.1 degrees radius

# 3. Perform the cone search
# 'catalog' parameter specifies which catalog to search.
# The 'Cone' spatial parameter indicates a cone search.
results_table = Irsa.query_region(
    center_coord,
    radius=radius,
    catalog='fp_psc', # Example catalog (2MASS All-Sky Point Source Catalog (PSC)')
    spatial='Cone'
)

# 4. Process the results (results are returned as an Astropy Table)

output_path = Path.home() / "open_clusters" / "data" / "out"

results_table.write(str(output_path)+'/output_catalog_search.ecsv', format='ecsv', overwrite=True)

print(f"Found {len(results_table)} objects within the search radius.")
if len(results_table) > 0:
    print("First 5 results:")
    print(results_table[:5]['ra', 'dec']) # Display selected columns
