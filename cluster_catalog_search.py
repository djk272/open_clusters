#!/usr/bin/env python3

import sys
from pathlib import Path
import glob
import math
import json
import warnings
import matplotlib.pyplot as plt
import scipy
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.table import Table, hstack, QTable
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astroquery.ipac.irsa import Irsa
from astroquery.mast import Mast
from astroquery.vizier import Vizier
#from isochrones.mist import MISTIsochroneGrid
#import ezpadova




print("Starting cluster Catalog Search.\n")
# Path to input files

INITAL_INPUT_CSV = Path.home() / "open_clusters" / "data" / "out" # / "NGC2682_M67_result_with_membership.csv"


# Get a list of all files in the directory
# Get a list of filenames matching a pattern

cluster_names_list = [file.stem for file in INITAL_INPUT_CSV.glob('*.csv')]
# [file_path.name for file_path in INPUT_CSV.glob("*.csv")]

if len(cluster_names_list) < 1:
    print("There are no CSV files in this directory.\n")
    sys.exit()
else:
    print("CSV files found!\n")
#/end of if

# For loop to construct input filepath for each cluster for searching
for cluster in cluster_names_list:
    print(f"Cluster: {cluster}\n")
    INPUT_CSV = f"{INITAL_INPUT_CSV}/{cluster}.csv"
    print(f"The filepath is {INPUT_CSV}\n")


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

    # Cluster members as a table
    cmember_table = QTable.read(INPUT_CSV)

    # RA and Dec values for cluster
    cmember_ra = cmember_table["ra"]#*u.degree
    cmember_dec = cmember_table["dec"]#*u.degree


    # 1. Define the search center, specify coordinates directly
    center_coord = SkyCoord(ra=center_ra, dec=center_dec, unit='deg', frame='icrs')
    print(f"Searching around coordinates: {center_coord.ra.deg} RA, {center_coord.dec.deg} Dec\n")

    # 2. Define the search radius #Find the radius using the the center which of the 2 differences is larger, that will be the radius for the cone
    #search. (decmax-center_dec), (ramax-center_ra).
    radius = math.sqrt((dec_max_value-center_dec)**2 + (ra_max_value-center_ra)**2) * u.deg # 0.1 degrees radius

    # 3. Perform the cone search
    
    print(f"Cone searching 2MASS All-Sky Point Source Catalog (PSC) for, {cluster}, radius: {radius} \n")
    # 'catalog' parameter specifies which catalog to search.
    # The 'Cone' spatial parameter indicates a cone search.
    results_table_1 = Irsa.query_region(
        center_coord,
        radius=radius,
        catalog='fp_psc', # Example catalog (2MASS All-Sky Point Source Catalog (PSC))
        spatial='Cone'
    )

    # 4.1 Process the results (results are returned as an Astropy Table)

    output_path = INITAL_INPUT_CSV / cluster
    # Create directory (and parents if needed)
    # exist_ok=True prevents an error if it already exists
    output_path.mkdir(parents=True, exist_ok=True)

    results_table_1.write(f"{output_path}/{cluster}_fp_psc_output_catalog_search.csv", format='csv', overwrite=True)

    print(f"Found {len(results_table_1)} objects within the fp_psc search radius.\n")
    if len(results_table_1) > 0:
        print("First 5 results:")
        print(f"{results_table_1[:5]['ra', 'dec']}\n") # Display selected columns
    #/end of if

    # 4.2 import files into astopy table objects
    cone_search_table_1 = results_table_1

    # Create objects for RA and Dec
    cone_search_1_ra = cone_search_table_1["ra"]#*u.degree
    cone_search_1_dec = cone_search_table_1["dec"]#*u.degree

    print(f"Cone searching Galex MAST server for, {cluster}, radius: {radius} \n")
    # JSON library only supports basic types like strings, integers, and dictionaries. 
    # Complex objects like Quantity (Ra and Dec have units) must be converted to a serializable format 
    # before saving or sending them as JSON. 
    json_ra = float(center_ra)
    json_dec = float(center_dec)
    json_radius = radius.value

    # Mast_query Sends a request to the MAST server and returns the response(in JSON).
    results_table_2 = Mast.mast_query('Mast.Galex.Catalog', 
                                    ra=json_ra, 
                                    dec=json_dec, 
                                    radius=json_radius
                                    )

    # 4.3 Process the results (results are returned as an Astropy Table)

    output_path = Path.home() / "open_clusters" / "data" / "out" / cluster
    # Create directory (and parents if needed)
    # exist_ok=True prevents an error if it already exists
    output_path.mkdir(parents=True, exist_ok=True)

    results_table_2.write(f"{output_path}/{cluster}_galex_output_catalog_search.csv", format='csv', overwrite=True)

    print(f"Found {len(results_table_2)} objects within the galex search radius.\n")
    if len(results_table_2) > 0:
        print("First 5 results:")
        print(f"{results_table_2[:5]['ra', 'dec']}\n") # Display selected columns
    #/end of if

    # 4.4 import files into astopy table objects
    cone_search_table_2 = results_table_2

    # Create objects for RA and Dec
    cone_search_2_ra = cone_search_table_2["ra"]#*u.degree
    cone_search_2_dec = cone_search_table_2["dec"]#*u.degree

    iter_list = [cone_search_table_1,cone_search_table_2]
    count = 0

    try:
        # the .read() below produces some warnings that we can safely ignore
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)

        for i in iter_list:
            count += 1
            if count == 1:
                catalog = 'fp_psc'
                cone_search_ra = cone_search_1_ra
                cone_search_dec = cone_search_1_dec
                cone_search_table = cone_search_table_1
            else:
                catalog = 'galex'
                cone_search_ra = cone_search_2_ra
                cone_search_dec = cone_search_2_dec
                cone_search_table = cone_search_table_2
            #/end of if

            # 4.5 We can now create a single SkyCoord object to represent all of the sources from our cone search results
            cone_search_coords = SkyCoord(cone_search_ra, cone_search_dec, unit="deg", frame='icrs')
            cmember_coords = SkyCoord(cmember_ra, cmember_dec, unit="deg", frame='icrs')
            print(f"4.5 Found {len(cone_search_coords)} cone search objects.\n   Found {len(cmember_coords)} cluster member objects.\n")

            # 5 matching the catalogs
            max_sep = 2 * u.arcsec
            idx_cmember, sep2d, _ = cmember_coords.match_to_catalog_sky(cone_search_coords)
            sep_constraint = sep2d < max_sep
            catalog_matches = cone_search_table[idx_cmember[sep_constraint]]

            print(f"5. Found {len(catalog_matches)} matched objects between the catalogs.")
            if len(catalog_matches) > 0:
                print("5. First 5 results:")
                print(catalog_matches[:5]) # Display selected columns
            #/end of if

            # 6 Look at the distribution of separations (in degrees) for all of the cross-matched sources

            plt.figure()
            plt.hist(sep2d[sep_constraint].arcsec, histtype="step", bins=np.logspace(-2, 2.0, 64))
            plt.xlabel("separation [arcsec]")
            plt.xscale("log")
            plt.yscale("log")
            plt.tight_layout()

            # 7 Save the figure to a file
            plt.savefig(f"{output_path}/{catalog}_{cluster}_sep_plot.png", dpi=150)
            plt.close()
            print(f"7.  Saved: {catalog}_{cluster}_sep_plot.png\n")



            # 8 Combine relevant information from both tables, only matched results
            # Adds 'column_to_add' from df2 to df1 where 'id' matches

            catalog_matches = hstack([cmember_table, catalog_matches])


            # 9 Process the results (results are returned as an Astropy Table)
            # write the reuslts to a CSV file

            print(f"9. Found {len(catalog_matches)} objects after concating the columns.")
            if len(catalog_matches) > 0:
                print("9. First 5 results:")
                print(catalog_matches[:5]) # Display selected columns
            #/end of if

            catalog_matches.write(f"{output_path}/{catalog}_{cluster}_output_matches.csv", format='csv', delimiter = ',', overwrite=True)


            # 10 Extracting parallax measurements from Gaia
            parallax = catalog_matches["parallax"]  # in milliarcseconds
            distance_array = 1 / parallax  # distance array from parallax (kiloparsecs)
            ra = catalog_matches["ra_1"]  # right ascension
            dec = catalog_matches["dec_1"]  # declination

            # 10.1 Sorting data
            sorted_indices = np.argsort(distance_array)
            ra_sorted = ra[sorted_indices]
            dec_sorted = dec[sorted_indices]
            distance_sorted = distance_array[sorted_indices]

            # 10.2 Creating a histogram of the parallax measurements and RA/DEC map
            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))

            # 10.3 Histogram
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

            # 10.4 RA/DEC map
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

            plt.savefig(f"{output_path}/{catalog}_{cluster}_distance_check_plot.png", dpi=150)
            plt.close()
            print(f"10.4.  Saved: {catalog}_{cluster}_distance_check_plot.png\n")

            # 10.5 Finding the median parallax of the data
            cluster_median = np.nanmedian(parallax)

            # 11 Calculating distance using stellar parallax formula
            distance = 1.0 / cluster_median

            print(f"11. Median Parallax: {cluster_median} milliarcseconds")
            print(f"11. Distance to Cluster: {distance:.4f} kpc")

            # 11.1 Distance to M13
            d = distance * 10**3  # parsecs


            # 11.2 Converting absolute magnitudes from isochrones to apparent
            def app_mag(M):
                return M + 5 * np.log10(d / 10)


            # 11.2 Importing isochrone .dat file
            parsec_data = pd.read_csv(
                Path.home() / "open_clusters" / "data" / "isochrones.dat", sep=r"\s+", comment="#", header=None
            )

            # 11.3 Converting the data file to an array
            parsec_data = np.array(parsec_data)


            # 11.4 The isochrone data file is broken up into age bins
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

            # 11.5 Setting up a dictionary
            age_data = {}

            # 11.6 Looping through .txt file to separate by age
            for i, age in enumerate(log_age, start=1):
                age_data[f"age{i}"] = parsec_data[parsec_data[:, 2] == age]

            # 11.7 Extracting g and r magnitudes, applying distance modulus formula
            # Setting up arrays
            r_data = []
            g_data = []
            color_data2 = []
            #/end of for loop

            # 11.8 Sending g and r magnitudes and color to an array
            for i in range(1, len(log_age) + 1):
                data = age_data[f"age{i}"]
                g_mag = app_mag(data[:, 29])
                r_mag = app_mag(data[:, 30])
                color = g_mag - r_mag

                r_data.append(r_mag)
                g_data.append(g_mag)
                color_data2.append(color)
            #/end of for loop

            # 11.9 Plotting the isochrones
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
            #/end of for loop

            for i, j, k, r in zip(color_data2, r_data, log_age, color):
                ax.plot(
                    i,
                    j,
                    alpha=0.5,
                    label=f"{(10**k) / 1e9:.0f} billion years old",
                    color=r,
                )
            #/end of for loop

            # 12 With our cross-match done, 
            # we can now make Gaia+2MASS color–magnitude diagrams of our candidate NGC 188 members 
            # using the information returned by the cross-match:

            # 12.1 combine the color information columns from 2mass and cluster member table.
            # bluer - redder color, for the color. "Mag" is for the magnitute "G" optical, "J" infrared

            if count == 1:
                Gmag = catalog_matches["phot_g_mean_mag"]
                Bmag = catalog_matches["bp_rp"]
            else:
                # Jmag = catalog_matches["j_m"]  # note that we use the index array returned above
                # H = catalog_matches["h_m"]
                # K = catalog_matches["k_m"]
                FUVmag = catalog_matches["fuv_mag"]
                NUVmag = catalog_matches["nuv_mag"]
            #/end of if

            #fig, axes = plt.subplots(1, 2, figsize=(10, 5))

            if count == 1:
                #ax = axes[1]
                ax.set_title(f"CMD of {cluster} sources from GAIA with Isochrones")
                ax.scatter(Bmag, Gmag, marker="x", color="k", linewidth=0.9, alpha=0.5, s=20)
                ax.set_xlabel("B-R")
                ax.set_ylabel("$G$")
                ax.invert_yaxis()
                # ax.set_xlim(0, 1.5)
                # ax.set_ylim(14, 12)  # backwards because magnitudes!
                ax.set_xlim(0, 3)
                ax.set_ylim(20, 10)  # backwards because magnitudes!
                ax.legend(fontsize=8)

                # # 12.2 Adding annotation for main sequence turnoff
                # circle = plt.Circle(
                #     (0.80, 12.50),
                #     0.25,
                #     color="red",
                #     fill=False,
                #     linewidth=2,
                #     label="Main Sequence Turnoff",
                # )
                # ax.add_patch(circle)
                # ax.annotate(
                #     "Main Sequence Turnoff",
                #     xy=(1.0, 12.50),
                #     xytext=(1.0 + 0.05, 12.50 - 0.2),
                #     arrowprops=dict(arrowstyle="->", color="red"),
                #     fontsize=12,
                #     color="red",
                # )

                fig.tight_layout()

                fig.savefig(f"{output_path}/{catalog}_{cluster}_cm_diagram_figure.png", dpi=150)
                print(f"12.2.  Saved: {catalog}_{cluster}_cm_diagram_figure.png\n \n")
            else:
                # ax = axes[0]
                ax.set_title(f"CMD of matched {cluster} sources in GALEX")
                ax.scatter(FUVmag.value - NUVmag.value, FUVmag.value, marker="o", color="k", linewidth=0, alpha=0.5, s=1)
                ax.set_xlabel("$FUV - NUV$")
                ax.set_ylabel("$FUVmag$")
                ax.set_xlim(-2, 4)
                ax.set_ylim(20, 4)  # backwards because magnitudes!
                ax.legend(fontsize=8)

                fig.tight_layout()

                fig.savefig(f"{output_path}/{catalog}_{cluster}_cm_diagram_figure.png", dpi=150)
                print(f"12.2.  Saved: {catalog}_{cluster}_cm_diagram_figure.png\n \n")
            #/end of if 
    except:
        pass # doing nothing on exception
        #/end of for loop
#/end of for loop

print(f"SUCCESS: Matching for {cluster_names_list} is done!")