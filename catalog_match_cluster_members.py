#!/usr/bin/env python3

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa import Irsa # Example service, many others available

# 1. Define the search center (e.g., the Andromeda Galaxy M31)
# The coordinates can be specified by name (resolved online) or as RA/Dec values.
#center_coord = SkyCoord.from_name("M31", frame="icrs")
#print(f"Searching around coordinates: {center_coord.ra.deg} RA, {center_coord.dec.deg} Dec\n")

# Alternatively, specify coordinates directly
center_coord = SkyCoord(ra=10.68470833, dec=41.26875, unit='deg', frame='icrs')
print(f"Searching around coordinates: {center_coord.ra.deg} RA, {center_coord.dec.deg} Dec\n")

# 2. Define the search radius
radius = 0.1 * u.deg # 0.1 degrees radius

# 3. Perform the cone search
# 'catalog' parameter specifies which catalog to search.
# The 'Cone' spatial parameter indicates a cone search.
results_table = Irsa.query_region(
    center_coord,
    radius=radius,
    catalog='ztf_objects_dr20', # Example catalog (WISE All-sky source catalog)
    spatial='Cone'
)

# 4. Process the results (results are returned as an Astropy Table)
print(f"Found {len(results_table)} objects within the search radius.")
if len(results_table) > 0:
    print("First 5 results:")
    print(results_table[:5]['ra', 'dec']) # Display selected columns
