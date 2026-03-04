import os
import re
import warnings
from datetime import datetime, timedelta
from typing import Optional, List, Tuple, Union

import numpy as np
import pandas as pd
# import urllib
import urllib.request
import xarray as xr
# import requests
import bz2
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
from urllib.parse import quote

# from compression import bz2

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# Constants
SEXTANT_BASE_URL = "ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic"
ODATIS_MR_BASE_URL = "https://tds-odatis.aviso.altimetry.fr/thredds/{access_type}/dataset-l3-ocean-color-odatis-mr-v_1_0.xml/FRANCE"

# Helper functions
def check_date_range(dates: List[str]) -> Tuple[datetime, datetime]:
    """Check and parse date range."""
    if len(dates) > 2:
        raise ValueError("Please provide only one or two dates for the date range.")
    start_date = datetime.strptime(dates[0], "%Y-%m-%d")
    end_date = start_date if len(dates) == 1 else datetime.strptime(dates[1], "%Y-%m-%d")
    return start_date, end_date

def check_time_step(time_step: Optional[str]) -> str:
    """Check and standardize time step."""
    if time_step is None:
        warnings.warn("No time-step chosen, defaulting to daily data.")
        return "day"
    time_step = time_step.lower()
    if time_step in ("day", "daily"):
        return "day"
    elif time_step in ("month", "monthly"):
        return "month"
    elif time_step in ("8d", "8day", "8-day", "8-daily", "weekly"):
        return "8-day"
    else:
        raise ValueError("'dl_time_step' value not recognised, please choose one of: 'day', '8-day', 'month'")

def check_bbox(bbox: Optional[List[float]]) -> Optional[List[float]]:
    """Check and validate bounding box."""
    if bbox is None:
        return None
    if len(bbox) != 4:
        raise ValueError("Please provide a bounding box in the form of a list with 4 values: [lonmin, lonmax, latmin, latmax]")
    if any(l < -180 or l > 180 for l in bbox[:2]):
        raise ValueError("Please ensure that the longitude values in 'dl_bbox' are between -180 and 180.")
    if any(l < -90 or l > 90 for l in bbox[2:]):
        raise ValueError("Please ensure that the latitude values in 'dl_bbox' are between -90 and 90.")
    if bbox[0] >= bbox[1]:
        raise ValueError("Please ensure that the first value in 'dl_bbox' (lonmin) is less than the second value (lonmax).")
    if bbox[2] >= bbox[3]:
        raise ValueError("Please ensure that the third value in 'dl_bbox' (latmin) is less than the fourth value (latmax).")
    return bbox

def infer_product(var: str) -> str:
    """Infer product from variable."""
    var = var.upper()
    if var in ("SPM", "SPIM", "CHL", "CHLA"):
        warnings.warn("No data product chosen, defaulting to 'SEXTANT' for SPM and Chl a data.")
        return "SEXTANT"
    elif var in ("CDOM", "RRS", "NRRS", "T", "TUR", "SST"):
        warnings.warn("No data product chosen, defaulting to 'ODATIS-MR'.")
        return "ODATIS-MR"
    else:
        raise ValueError("Variable not available, please check the value given for 'dl_var'")

def infer_sensor(var: str, product: str) -> str:
    """Infer sensor from variable and product."""
    if product == "ODATIS-MR":
        if var == "SST":
            warnings.warn("No sensor chosen, defaulting to 'MODIS' because 'dl_var = SST'.")
            return "MODIS"
        else:
            warnings.warn("No sensor chosen, defaulting to 'OLCI-A'.")
            return "OLCI-A"
    return None

def infer_correction(sensor: str) -> str:
    """Infer correction from sensor."""
    if sensor == "MODIS":
        warnings.warn("No atmospheric correction chosen, defaulting to 'nirswir' to match 'dl_sensor = MODIS'.")
        return "nirswir"
    elif sensor in ("MERIS", "OLCI-A", "OLCI-B"):
        warnings.warn(f"No atmospheric correction chosen, defaulting to 'polymer' to match 'dl_sensor = {sensor}'.")
        return "polymer"
    return None

def get_sensor_dates(sensor: str) -> Tuple[datetime, datetime]:
    """Get floor and ceiling dates for sensor."""
    if sensor == "MODIS":
        return datetime(2002, 7, 4), datetime(2024, 12, 31)
    elif sensor == "MERIS":
        return datetime(2002, 6, 19), datetime(2012, 4, 8)
    elif sensor == "OLCI-A":
        return datetime(2016, 4, 26), datetime(2024, 12, 31)
    elif sensor == "OLCI-B":
        return datetime(2018, 5, 15), datetime(2024, 12, 31)
    else:
        return datetime(1998, 1, 1), datetime.now() - timedelta(days=7)

def download_nc(
    dl_var: str,
    dl_dates: List[str],
    dl_product: Optional[str] = None,
    dl_correction: Optional[str] = None,
    dl_sensor: Optional[str] = None,
    dl_time_step: Optional[str] = None,
    dl_bbox: Optional[List[float]] = None,
    username: Optional[str] = None,
    password: Optional[str] = None,
    output_dir: str = ".",
    overwrite: bool = False,
) -> None:
    """
    Download NetCDF files for given parameters.
    """
    # Check and parse inputs
    start_date, end_date = check_date_range(dl_dates)
    dl_time_step = check_time_step(dl_time_step)
    dl_bbox = check_bbox(dl_bbox)
    dl_product = infer_product(dl_var) if dl_product is None else dl_product.upper()
    dl_sensor = infer_sensor(dl_var, dl_product) if dl_sensor is None else dl_sensor.upper()
    dl_correction = infer_correction(dl_sensor) if dl_correction is None else dl_correction.lower()

    # Check sensor dates
    floor_date_sensor, ceiling_date_sensor = get_sensor_dates(dl_sensor)
    if start_date < floor_date_sensor or end_date > ceiling_date_sensor:
        raise ValueError(
            f"The chosen date range is outside of the available data range for: {dl_product} - {dl_sensor} - {dl_correction}.\n"
            f"Available data range is from {floor_date_sensor} to {ceiling_date_sensor}. Please adjust your date range accordingly."
        )

    # Iterate over each date in the range
    current_date = start_date
    while current_date <= end_date:
        dl_date = current_date
        dl_date_flat = dl_date.strftime("%Y%m%d")
        dl_year = dl_date.strftime("%Y")
        dl_month = dl_date.strftime("%m")
        dl_day = dl_date.strftime("%d")
        url_year_doy = f"{dl_year}/{dl_date.strftime('%j')}"

        # Build URL and download
        if dl_product == "SEXTANT":
            url_product_stub = "EUR-L4-SPIM-ATL-v01" if dl_var.upper() in ("SPM", "SPIM") else "EUR-L4-CHL-ATL-v01"
            url_product = f"{url_product_stub}/{url_year_doy}"
            file_name = f"{dl_date_flat}-{url_product_stub}-fv01-OI.nc.bz2"
            url_final = f"{SEXTANT_BASE_URL}/{url_product}/{file_name}"
            file_name_full = os.path.join(output_dir, file_name)

        elif dl_product == "ODATIS-MR":
            username_fix = quote(username, safe="")
            password_fix = quote(password, safe="")
            access_type = "dodsC" if dl_bbox else "fileServer"
            url_base = ODATIS_MR_BASE_URL.format(access_type=access_type)
            url_base = url_base.replace("https://", f"https://{username_fix}:{password_fix}@")

            dl_sensor_flat = dl_sensor.lower().replace("-", "")
            dl_sensor_chunk = dl_sensor.replace("CI-", "") if "OLCI" in dl_sensor else dl_sensor[:3]
            dl_correction_flat = dl_correction.lower().replace("-", "")
            dl_correction_chunk = "PO" if dl_correction == "polymer" else "NS"

            if dl_var.upper() in ("CDOM",):
                dl_var_chunk = dl_var.upper()
            elif dl_var.upper() in ("CHL", "CHLA"):
                dl_var_chunk = "CHL-OC5"
            elif dl_var.upper() in ("RRS", "NRRS"):
                dl_var_chunk = "NRRS555" if dl_sensor == "MODIS" else "NRRS560"
            elif dl_var.upper() in ("SPM", "SPIM"):
                dl_var_chunk = "SPM-G"
            elif dl_var.upper() in ("T", "TUR"):
                dl_var_chunk = "T-FNU"
            elif dl_var.upper() in ("SST",):
                dl_var_chunk = "SST-NIGHT"
            else:
                raise ValueError(f"Variable {dl_var} not available for {dl_sensor}.")

            dl_time_step_chunk = "DAY" if dl_time_step == "day" else "8D" if dl_time_step == "8-day" else "MO"
            file_name = f"L3m_{dl_date_flat}__FRANCE_03_{dl_sensor_chunk}_{dl_var_chunk}-{dl_correction_chunk}_{dl_time_step_chunk}_00.nc"
            url_product = f"{dl_correction_flat}/{dl_sensor_flat}/{dl_time_step}/{dl_year}/{dl_month}/{dl_day}"
            url_final = f"{url_base}/{url_product}/{file_name}"
            file_name_full = os.path.join(output_dir, file_name)

        else:
            raise ValueError(f"Product {dl_product} not supported.")

        # Download logic
        if os.path.exists(file_name_full) and not overwrite:
            warnings.warn(f"{file_name_full} already exists. Set 'overwrite=True' to force the download.")
        else:
            os.makedirs(output_dir, exist_ok=True)
            try:
                # print(f"Downloading {file_name} from {url_final}...")
                # Attempt to download the file
                urllib.request.urlretrieve(url_final, file_name_full)
                print(f"File downloaded at: {file_name_full}")
                # Extract the .bz2 file if successful
                with open(file_name_full, 'rb') as source, open(file_name_full[:-4], 'wb') as dest:
                    dest.write(bz2.decompress(source.read()))
                try:
                    os.remove(file_name_full)
                    print(f"Compressed file removed: {file_name_full}")
                except Exception as e:
                    print(f"Failed to remove compressed file: {e}")

            except Exception as e:
                print(f"Failed to download or extract {url_final}: {e}")

        # Move to next time step
        if dl_time_step == "day":
            current_date += timedelta(days=1)
        elif dl_time_step == "8-day":
            current_date += timedelta(days=8)
        elif dl_time_step == "month":
            current_date = (current_date.replace(day=28) + timedelta(days=4)).replace(day=1)

def plot_nc(
    nc_file: str,
    bbox: Optional[List[float]] = None,
    plot_width: Optional[float] = None,
    plot_height: Optional[float] = None,
    output_dir: str = ".",
) -> None:
    """
    Plot a NetCDF file.
    """
    if not os.path.exists(nc_file):
        warnings.warn(f"{nc_file} cannot be found.")
        return

    # Open the NetCDF file
    ds = xr.open_dataset(nc_file)

    # Infer variable name and label
    file_snippets = os.path.basename(nc_file).split("_")
    if len(file_snippets) == 1:
        if "SPIM-ATL" in nc_file:
            nc_var_name = "analysed_spim"
        elif "CHL-ATL" in nc_file:
            nc_var_name = "analysed_chl_a"
        else:
            raise ValueError("File structure cannot be inferred from file name.")
    elif len(file_snippets) == 9:
        nc_var_name = f"{file_snippets[6]}_mean"
    else:
        raise ValueError("File structure cannot be inferred from file name.")

    if "spim" in nc_var_name.lower():
        var_label = "SPM [g m-3]"
    elif "chl" in nc_var_name.lower():
        var_label = "Chl a [mg m-3]"
    elif "cdom" in nc_var_name.lower():
        var_label = "CDOM [m-1]"
    elif "sst" in nc_var_name.lower():
        var_label = "SST [°C]"
    elif "t-fnu" in nc_var_name.lower():
        var_label = "Turbidity [FNU]"
    elif "nrrs560" in nc_var_name.lower():
        var_label = "Rrs at 560 nm [sr-1]"
    elif "nrrs555" in nc_var_name.lower():
        var_label = "Rrs at 555 nm [sr-1]"
    else:
        raise ValueError("Variable label cannot be inferred from variable name.")

    # Extract and prep time for plotting
    time = ds["time"].values
    date_time_obj = pd.to_datetime(time[0])
    # date_label = date_time_obj.strftime('%Y-%m-%d')
    
    # Subset NetCDF file to desired bounding box
    if(bbox is None):
        bbox = [ds.lon.min().item(), ds.lon.max().item(), ds.lat.min().item(), ds.lat.max().item()]
        print(f"No bounding box provided, using full extent: {bbox}")
            
    var_subset = ds[nc_var_name].where(
        (ds.lon >= bbox[0]) & (ds.lon <= bbox[1]) &
        (ds.lat >= bbox[2]) & (ds.lat <= bbox[3]),
        drop=True
    )

    # Create a figure and axes with a Plate Carree projection
    plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Plot the data
    var_subset.plot(
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
    plt.title(f'Map of {nc_var_name} on {date_time_obj.strftime("%Y-%m-%d")}')

    # Show the plot
    plt.show()

# Example usage
if __name__ == "__main__":
    # Example: Download SPM data for a date range
    download_nc(
        dl_var="SPM",
        dl_dates=["2025-09-01"],
        output_dir="../data/SEXTANT",
        overwrite=False,
    )

    # Example: Plot a NetCDF file
    plot_nc(
        nc_file="../data/SEXTANT/20250901-EUR-L4-SPIM-ATL-v01-fv01-OI.nc",
        bbox=[-3.5, -0.5, 44, 48],
        plot_width=5,
        plot_height=9,
        output_dir="../data/SEXTANT",
    )
