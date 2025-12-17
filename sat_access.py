#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import urllib.request
import bz2
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from datetime import datetime, timedelta
import netCDF4
import imageio.v2 as imageio

from matplotlib.animation import FuncAnimation, PillowWriter

def download_and_plot(dl_var, dl_dates, bbox, output_dir, plot_var, overwrite):
    """
    Downloads a NetCDF file from a specified URL and plots a variable as a map.

    Parameters:
    - dl_var: Variable to download (e.g., "SPM", "SPIM", "CHLA").
    - dl_date: Date in YYYY-MM-DD format.
    - bbox: Bounding box as [lon_min, lon_max, lat_min, lat_max].
    - output_dir: Directory to save the downloaded file and plot.
    - overwrite: Whether to overwrite existing files.
    """

    # === Vérification dates ===
    if len(dl_dates) > 2:
        raise ValueError("Please provide only one or two dates for the date range.")

    start_date = datetime.strptime(dl_dates[0], '%Y-%m-%d')
    if len(dl_dates) == 2:
        end_date = datetime.strptime(dl_dates[1], '%Y-%m-%d') 
    else:
        end_date = start_date

    # === BOUCLE TÉLÉCHARGEMENT JOURNALIER ===
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
            raise ValueError("Variable not yet available")

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

        # Download file
        if os.path.exists(file_name_full) and not overwrite:
            print(f"{file_name_full} already exists. Set --overwrite True to force the download.")
        else:
            # Create directory if it does not yet exist
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            try:
                urllib.request.urlretrieve(url_final, file_name_full)
                print(f"File downloaded at: {file_name_full}")
                # Extract the .bz2 file
                with open(file_name_full, 'rb') as source, open(file_name_full[:-4], 'wb') as dest:
                    dest.write(bz2.decompress(source.read()))
                print(f"File extracted at: {file_name_full[:-4]}")
                
                try:
                    os.remove(file_name_full)
                    print(f"Fichier compressé supprimé : {file_name_full}")
                except Exception as e:
                    print(f"Impossible de supprimer le fichier compressé : {e}")
                    
            except Exception as e:
                print(f"Failed to download or extract {url_final}: {e}")

        # Move to the next day
        current_date += timedelta(days=1)
    
    
    # === CONCATÉNATION NETCDF "AU FIL DE L'EAU" ===
    print("\n=== Début de la concaténation au fil de l'eau ===")

    output_final = os.path.join(
        output_dir,
        f"{dl_var}_DAILY_{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}.nc"
    )

    # Créer le fichier final vide
    nc_out = netCDF4.Dataset(output_final, 'w', format='NETCDF4')

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
    print(f"\n=== Fichier final créé : {output_final} ===")

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download SPM or chl a NetCDF files and plot them as a map.")
    parser.add_argument("--variable", type=str, required=True, help="Variable to download (e.g., 'SPM', 'CHLA')")
    parser.add_argument("--daterange", type=str, nargs='+', required=True, 
                        help="Date range for desired data in YYYY-MM-DD format. Provide only one or two values. Note that if two dates are provided, only the first will be used for plotting.")
    parser.add_argument("--outputdir", type=str, required=True, help="Directory to save the downloaded file and plot")
    parser.add_argument("--overwrite", type=bool, default=False, help="Overwrite existing files. Default = False")
    parser.add_argument("--plot", type=bool, default=False, help="Whether or not to plot the downloaded data. Default = False")
    parser.add_argument("--boundingbox", type=float, nargs=4, required=False, help="The bounding box for plotting. Must be given as: lonmin, lonmax, latmin, latmax")

    # ← ici on ajoute les arguments par défaut si aucun n’est passé
    import sys
    if len(sys.argv) == 1:
        sys.argv += [
            "--variable", "SPM",
            "--daterange", "1999-01-01", "1999-01-02",
            "--outputdir", "C:/Users/rpollet/Downloads/output/",
            "--plot", "True",
            "--overwrite", "False",
            "--boundingbox", "2","7","41","44"
        ]

    args = parser.parse_args()
    download_and_plot(args.variable, args.daterange, args.boundingbox, args.outputdir, args.plot, args.overwrite)


# ###### ANIMATION ######
ds = xr.open_dataset("C:/Users/rpollet/Downloads/output/SPM_DAILY_19990101_20241231.nc")
nc_var_name = "analysed_spim"
var = ds[nc_var_name].sel(lat=slice(41.5, 44), lon=slice(2.5, 7))
name = "Inorganic suspended matters (mg/m3)"

times = var.time.values

fig = plt.figure(figsize=(10, 6), dpi=120)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=":")
ax.gridlines(draw_labels=True, ls=':')

v0 = var.isel(time=0) # Première frame

# Raster initial
mesh = ax.pcolormesh(v0.lon, v0.lat, v0.values,
    cmap="viridis", vmin=0, vmax=16,transform=ccrs.PlateCarree())
cbar = plt.colorbar(mesh, ax=ax, shrink=0.8, pad=0.1)
cbar.set_label(name)

title = ax.set_title(f"Map of {nc_var_name} on {np.datetime_as_string(times[0], unit='D')}")

# --- Fonctions Animation ---
def update(frame):
    t = times[frame]
    data = var.sel(time=t).values
    mesh.set_array(data.ravel())   # mise à jour rapide du raster
    title.set_text(f"Map of {nc_var_name} on {np.datetime_as_string(t, unit='D')}")
    return mesh,

anim = FuncAnimation(fig, update, frames=len(times), blit=False)

anim.save("C:/Users/rpollet/Downloads/output/animation_spm_1999-2024.gif",
          writer=PillowWriter(fps=5))
plt.show()
