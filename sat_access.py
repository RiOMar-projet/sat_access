#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# Authors: Robert Schlegel, Romane Pollet
# Date: January 2026
# Description: This script downloads SPM or chlorophyll-a NetCDF files from a specified URL
#              for given date ranges, optionally concatenates them, and plots or animates the data.
"""

import os
import argparse
import urllib.request
import bz2
import netCDF4
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np

from datetime import datetime, timedelta
from matplotlib.animation import FuncAnimation, PillowWriter

def download_and_plot(dl_var, dl_dates, output_dir, overwrite, concatenate, plot_var, gif_var, bbox):
    """
    Downloads a series of NetCDF files from a specified URL and plots the variable or animates as a GIF.

    Parameters:
    - dl_var: Variable to download (e.g., "SPM", "SPIM", "CHLA").
    - dl_dates: Date range in YYYY-MM-DD YYYY-MM-DD format. Also accepts a single YYYY-MM-DD date.
    - output_dir: Directory to save the downloaded file and plot.
    - overwrite: Whether to overwrite existing files.
    - concatenate: Whether or not to concatenate all of the downloaded files into one.
    - plot_var: Whether or not to plot the downloaded data.
    - anim_var: Whether or not to animate the downloaded data.
    - bbox: Bounding box as [lon_min, lon_max, lat_min, lat_max].
    """

    # === Vérification dates ===
    if len(dl_dates) > 2:
        raise ValueError("Please provide only one or two dates for the date range.")

    start_date = datetime.strptime(dl_dates[0], '%Y-%m-%d')
    if len(dl_dates) == 2:
        end_date = datetime.strptime(dl_dates[1], '%Y-%m-%d') 
    else:
        end_date = start_date

    # === BOUCLE TÉLÉCHARGEMENT JOURNALIER === #

    # Check if concatenated file already exists and contains the requested dates
    concat_filestub = f"{dl_var}_DAILY"
    concat_files = [f for f in os.listdir(output_dir) if concat_filestub in f and f.endswith('.nc')]
    if len(concat_files) == 1:
        concat_filename = concat_files[0]
        output_concat = os.path.join(output_dir, concat_filename)
        parts = concat_filename.split('_')
        file_start = datetime.strptime(parts[2], '%Y%m%d')
        file_end = datetime.strptime(parts[3].split('.')[0], '%Y%m%d')
    elif len(concat_files) > 1:
        print(f"Multiple concatenated files found: {concat_files}. Checking which one contains the requested dates...")
        concat_filename = "none"
        output_concat = "none"
        file_start = datetime.strptime('00010101', '%Y%m%d')
        file_end = datetime.strptime('00010101', '%Y%m%d')
        
        for file in concat_files:
            parts = file.split('_')
            if len(parts) == 4:
                try:
                    file_start_candidate = datetime.strptime(parts[2], '%Y%m%d')
                    file_end_candidate = datetime.strptime(parts[3].split('.')[0], '%Y%m%d')
                    
                    # Check if this file contains the requested dates
                    if start_date >= file_start_candidate and end_date <= file_end_candidate:
                        concat_filename = file
                        output_concat = os.path.join(output_dir, file)
                        file_start = file_start_candidate
                        file_end = file_end_candidate
                        print(f"Found matching file: {concat_filename}")
                        break
                except (ValueError, IndexError):
                    continue
        
        if concat_filename == "none":
            print(f"No suitable concatenated file found among: {concat_files}. Will download new data.")
    else:
        concat_filename = "none"
        output_concat = "none"
        file_start = datetime.strptime('00010101', '%Y%m%d')
        file_end = datetime.strptime('00010101', '%Y%m%d')
        # print("No matching concatenated files found.")
    
    # Check if requested dates are within file date range
    if start_date >= file_start and end_date <= file_end:
        concat_test = True
        print(f"Concatenated file {concat_filename} already contains the requested dates ({dl_dates[0]} to {dl_dates[-1]}). No files need to be downloaded.")
    else :
        concat_test = False
        print(f"Downloading data from {start_date} to {end_date}.")
        current_date = start_date
        while current_date <= end_date:
            dl_date = current_date.strftime('%Y-%m-%d')
            dl_date_flat = dl_date.replace("-", "")

            # Get year and day of year for URL
            url_year_doy = f"{dl_date[:4]}/{current_date.strftime('%j')}"

            # Get general URL based on desired variables
            if dl_var.upper() in ["SPM", "SPIM", "CHLA"]:
                url_base = "ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic"
            else:
                raise ValueError("Variable not yet available. Choose from: SPM, SPIM, CHLA")

            # Get product specifics
            if dl_var.upper() in ["SPM", "SPIM"]:
                file_name = f"{dl_date_flat}-EUR-L4-SPIM-ATL-v01-fv01-OI.nc.bz2"
                url_product = "EUR-L4-SPIM-ATL-v01"
                nc_var_name = "analysed_spim"
                var_label = "SPM [g m-3]"
                nc_file = os.path.join(output_dir, f"{dl_date_flat}-EUR-L4-SPIM-ATL-v01-fv01-OI.nc")
            elif dl_var.upper() == "CHLA":
                file_name = f"{dl_date_flat}-EUR-L4-CHL-ATL-v01-fv01-OI.nc.bz2"
                url_product = "EUR-L4-CHL-ATL-v01"
                nc_var_name = "analysed_chl_a"
                var_label = "chl a [mg m-3]"
                nc_file = os.path.join(output_dir, f"{dl_date_flat}-EUR-L4-CHL-ATL-v01-fv01-OI.nc")
            else:
                raise ValueError("Variable not yet available")

            # Assemble final URL
            url_final = f"{url_base}/{url_product}/{url_year_doy}/{file_name}"
            file_name_full = os.path.join(output_dir, file_name)

            # Download file(s)
            if os.path.exists(file_name_full) and not overwrite:
                print(f"{file_name_full} already exists. Set --overwrite True to force the download.")
            else:
                # Create directory if it does not yet exist
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                try:
                    # Attempt to download the file
                    urllib.request.urlretrieve(url_final, file_name_full)
                    print(f"File downloaded at: {file_name_full}")
                    # Extract the .bz2 file if successful
                    with open(file_name_full, 'rb') as source, open(file_name_full[:-4], 'wb') as dest:
                        dest.write(bz2.decompress(source.read()))
                    try:
                        os.remove(file_name_full)
                        print(f"Fichier compressé supprimé : {file_name_full}")
                    except Exception as e:
                        print(f"Impossible de supprimer le fichier compressé : {e}")
                    
                except Exception as e:
                    print(f"Failed to download or extract {url_final}: {e}")

            # Move to the next day
            current_date += timedelta(days=1)
    
    # === CONCATÉNATION NETCDF "AU FIL DE L'EAU" === #
    if concatenate:
        if not concat_test:
            print("\n=== Début de la concaténation au fil de l'eau ===")

            # Final name for concatenated files
            concat_filename = f"{dl_var}_DAILY_{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}.nc"
            output_concat = os.path.join(
                output_dir, concat_filename
            )

            # Créer le fichier final vide
            nc_out = netCDF4.Dataset(output_concat, 'w', format='NETCDF4')

            first_file = True
            current_date = start_date

            while current_date <= end_date:
                dl_date_flat = current_date.strftime('%Y%m%d')

                if dl_var.upper() in ["SPM", "SPIM"]:
                    nc_name = f"{dl_date_flat}-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"
                else:
                    nc_name = f"{dl_date_flat}-EUR-L4-CHL-ATL-v01-fv01-OI.nc"

                nc_path = os.path.join(output_dir, nc_name)

                if not os.path.exists(nc_path):
                    current_date += timedelta(days=1)
                    continue

                # Ouvrir fichier du jour
                ds = netCDF4.Dataset(nc_path, 'r')

                if first_file:
                    # Copier dimensions
                    for dim_name, dim in ds.dimensions.items():
                        if dim_name == "time":
                            nc_out.createDimension("time", None)  # illimité
                        else:
                            nc_out.createDimension(dim_name, len(dim))

                    # Copier variables
                    for var_name, var in ds.variables.items():
                        out_var = nc_out.createVariable(
                            var_name,
                            var.datatype,
                            var.dimensions,
                            zlib=True
                        )
                        # copier attributs
                        out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})

                    first_file = False

                # Append des données du jour
                time_len_current = nc_out.variables["time"].shape[0]

                for var_name, var in ds.variables.items():
                    if "time" in var.dimensions:
                        nc_out.variables[var_name][time_len_current:time_len_current+1] = var[:]
                    else:
                        if time_len_current == 0:
                            nc_out.variables[var_name][:] = var[:]

                ds.close()

                # Supprimer fichier individuel
                try:
                    os.remove(nc_path)
                    print(f"Supprimé : {nc_path}")
                except:
                    print(f"Impossible de supprimer : {nc_path}")

                current_date += timedelta(days=1)

            nc_out.close()
            print(f"\n=== Fichier final créé : {output_concat} ===")

    # === Plotting code === #
    if plot_var:

        print("Plotting...")

        # Determine which file to open
        if concatenate:
            nc_file = output_concat
        else:
            end_date_flat = end_date.strftime('%Y%m%d')
            if dl_var.upper() in ["SPM", "SPIM"]:
                nc_file = os.path.join(output_dir, f"{end_date_flat}-EUR-L4-SPIM-ATL-v01-fv01-OI.nc")
            else:
                nc_file = os.path.join(output_dir, f"{end_date_flat}-EUR-L4-CHL-ATL-v01-fv01-OI.nc")
        
        if not os.path.exists(nc_file):
            if os.path.exists(output_concat):
                nc_file = output_concat
            else:
                raise FileNotFoundError(f"Cannot find the required NetCDF file: {nc_file}")

        # Get product specifics
        if dl_var.upper() in ["SPM", "SPIM"]:
            nc_var_name = "analysed_spim"
            var_label = "SPM [g m-3]"
        elif dl_var.upper() == "CHLA":
            nc_var_name = "analysed_chl_a"
            var_label = "chl a [mg m-3]"
        else:
            raise ValueError("Variable not yet available")

        # Open the NetCDF file using xarray
        ds = xr.open_dataset(nc_file)

        # Select data for the requested date
        if len(dl_dates) == 2:
            date_to_select = dl_dates[1]
            print("Two dates provided; the last date will be used for plotting.")
        else:
            date_to_select = dl_dates[0]
        
        # date_to_select = np.datetime64(date_to_select)
        # Get date value as a label for plotting
        # date_value = ds.time.values[-1]
        date_time_obj = pd.to_datetime(date_to_select)
        date_label = date_time_obj.strftime('%Y-%m-%d')

        # Filter to a single date
        ds_subset = ds.sel(time=date_to_select, method='nearest')

        # Subset NetCDF file to desired bounding box
        if(bbox is None):
            bbox = [ds_subset.lon.min().item(), ds_subset.lon.max().item(), ds_subset.lat.min().item(), ds_subset.lat.max().item()]
            print(f"No bounding box provided, using full extent: {bbox}")
            
        var_subset = ds_subset[nc_var_name].where(
            (ds.lon >= bbox[0]) & (ds.lon <= bbox[1]) &
            (ds.lat >= bbox[2]) & (ds.lat <= bbox[3]),
            drop=True
        )

        # Create a figure and axes with a Plate Carree projection
        fig = plt.figure(figsize=(10, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Plot the data
        raster = var_subset.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap='viridis',
            add_colorbar=True,
            cbar_kwargs={'label': var_label}
        )

        # Add geospatial features
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.gridlines(draw_labels=True)

        # Add a title
        plt.title(f'Map of {nc_var_name} on {date_label}')

        # Show the plot
        plt.show()

    else:
        print("Plotting skipped.")

    # === ANIMATION === #
    if gif_var:
        if not concatenate:
            print("Set --concatenate True to be able to animate the data.")
        else:
            print("Animating...")

            # Get product specifics
            if dl_var.upper() in ["SPM", "SPIM"]:
                nc_var_name = "analysed_spim"
                var_label = "SPM [g m-3]"
            elif dl_var.upper() == "CHLA":
                nc_var_name = "analysed_chl_a"
                var_label = "chl a [mg m-3]"
            else :
                raise ValueError("Variable not yet available")

            ds = xr.open_dataset(output_concat)
            # var = ds[nc_var_name].sel(lat=slice(41.5, 44), lon=slice(2.5, 7))
            
            # Subset NetCDF file to desired bounding box
            if(bbox is None):
                bbox = [ds.lon.min().item(), ds.lon.max().item(), ds.lat.min().item(), ds.lat.max().item()]
                print(f"No bounding box provided, using full extent: {bbox}")
                
            var_subset = ds[nc_var_name].where(
                (ds.lon >= bbox[0]) & (ds.lon <= bbox[1]) &
                (ds.lat >= bbox[2]) & (ds.lat <= bbox[3]),
                drop=True
            )

            times = var_subset.time.values

            fig = plt.figure(figsize=(10, 6), dpi=120)
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.BORDERS, linestyle=":")
            ax.gridlines(draw_labels=True, ls=':')

            v0 = var_subset.isel(time=0) # Première frame

            # Raster initial
            mesh = ax.pcolormesh(v0.lon, v0.lat, v0.values,
                cmap="viridis", vmin=0, vmax=16,transform=ccrs.PlateCarree())
            cbar = plt.colorbar(mesh, ax=ax, shrink=0.8, pad=0.1)
            cbar.set_label(var_label)

            title = ax.set_title(f"Map of {nc_var_name} on {np.datetime_as_string(times[0], unit='D')}")

            # --- Fonctions Animation ---
            def update(frame):
                t = times[frame]
                data = var_subset.sel(time=t).values
                mesh.set_array(data.ravel()) # mise à jour rapide du raster
                title.set_text(f"Map of {nc_var_name} on {np.datetime_as_string(t, unit='D')}")
                return mesh,

            anim = FuncAnimation(fig, update, frames=len(times), blit=False)

            anim_filename = f"{output_dir}/{dl_var}_{dl_dates[0]}_to_{dl_dates[-1]}.gif"
            anim.save(anim_filename, writer=PillowWriter(fps=5))
            plt.show()
        
    else:
        print("Animating skipped.")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download SPM or chl a NetCDF files and plot them as a map.")
    parser.add_argument("-v", "--variable", type=str, required=True, help="Variable to download (e.g., 'SPM', 'CHLA')")
    parser.add_argument("-d", "--daterange", type=str, nargs='+', required=True, 
                        help="Date range for desired data in YYYY-MM-DD format. Provide only one or two values. Note that if two dates are provided, only the first will be used for plotting.")
    parser.add_argument("-od", "--outputdir", type=str, required=True, help="Directory to save the downloaded file and plot")
    parser.add_argument("-ov", "--overwrite", type=bool, default=False, help="Overwrite existing files. Default = False")
    parser.add_argument("-c", "--concatenate", type=bool, default=False, help="Whether or not to concatenate all downloaded files into one. Necessary to create animation. Default = False")
    parser.add_argument("-p", "--plot", type=bool, default=False, help="Whether or not to plot one day of the downloaded data. Default = False")
    parser.add_argument("-a", "--animate", type=bool, default=False, help="Whether or not to animate all of the downloaded data. Default = False")
    parser.add_argument("-bbox", "--boundingbox", type=float, nargs=4, required=False, help="The bounding box for plotting. Must be given as: lonmin, lonmax, latmin, latmax")

    # ← ici on ajoute les arguments par défaut si aucun n’est passé
    # import sys
    # if len(sys.argv) == 1:
    #     sys.argv += [
    #         "--variable", "SPM",
    #         "--daterange", "1999-01-01", "2024-12-31",
    #         "--outputdir", "C:/Users/your_name/Downloads/output/",
    #         "--plot", "True",
    #         "--overwrite", "False",
    #         "--boundingbox", "2","7","41","44"
    #     ]

    args = parser.parse_args()
    download_and_plot(args.variable, args.daterange, args.outputdir, args.overwrite, args.concatenate, args.plot, args.animate, args.boundingbox)

