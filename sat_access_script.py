import bz2
import os
import re
import shutil
import urllib.request
import warnings
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import List, Optional, Tuple
from urllib.parse import quote

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


warnings.filterwarnings("ignore")

SEXTANT_BASE_URL = "ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic"
ODATIS_MR_BASE_URL = (
    "https://tds-odatis.aviso.altimetry.fr/thredds/{access_type}/"
    "dataset-l3-ocean-color-odatis-mr-v_1_0.xml/FRANCE"
)
ODATIS_BBOX_LIMITS = [-7.8, 10.3, 41.2, 51.5]



def _parse_date(value: str) -> date:
    return datetime.strptime(value, "%Y-%m-%d").date()


def _month_start(value: date) -> date:
    return value.replace(day=1)


def _month_end(value: date) -> date:
    if value.month == 12:
        next_month = value.replace(year=value.year + 1, month=1, day=1)
    else:
        next_month = value.replace(month=value.month + 1, day=1)
    return next_month - timedelta(days=1)


def _next_month(value: date) -> date:
    if value.month == 12:
        return value.replace(year=value.year + 1, month=1, day=1)
    return value.replace(month=value.month + 1, day=1)


def _check_date_range(dl_dates: List[str], dl_time_step: str) -> Tuple[date, date]:
    if len(dl_dates) > 2:
        raise ValueError("Please provide only one or two dates for the date range.")

    if dl_time_step in {"day", "8-day"}:
        start_date = _parse_date(dl_dates[0])
        end_date = _parse_date(dl_dates[1]) if len(dl_dates) == 2 else start_date
    elif dl_time_step == "month":
        start_date = _month_start(_parse_date(dl_dates[0]))
        end_date = _month_start(_parse_date(dl_dates[1])) if len(dl_dates) == 2 else start_date
    else:
        raise ValueError(
            "'dl_time_step' value not recognised, please choose one of: 'day', '8-day', 'month'"
        )

    return start_date, end_date


def _check_time_step(time_step: Optional[str]) -> str:
    if time_step is None:
        warnings.warn("No time-step chosen, defaulting to daily data.")
        return "day"

    value = time_step.lower()
    if value in {"day", "daily"}:
        return "day"
    if value in {"month", "monthly"}:
        return "month"
    if value in {"8d", "8day", "8-day", "8-daily", "weekly"}:
        return "8-day"

    raise ValueError(
        "'dl_time_step' value not recognised, please choose one of: 'day', '8-day', 'month'"
    )


def _check_bbox(bbox: Optional[List[float]]) -> Optional[List[float]]:
    if bbox is None:
        return None
    if len(bbox) != 4:
        raise ValueError(
            "Please provide a bounding box in the form of a list with 4 values: [lonmin, lonmax, latmin, latmax]"
        )

    lon_min, lon_max, lat_min, lat_max = bbox
    if lon_min < -180 or lon_min > 180 or lon_max < -180 or lon_max > 180:
        raise ValueError("Please ensure that the longitude values in 'dl_bbox' are between -180 and 180.")
    if lat_min < -90 or lat_min > 90 or lat_max < -90 or lat_max > 90:
        raise ValueError("Please ensure that the latitude values in 'dl_bbox' are between -90 and 90.")
    if lon_min >= lon_max:
        raise ValueError(
            "Please ensure that the first value in 'dl_bbox' (lonmin) is less than the second value (lonmax)."
        )
    if lat_min >= lat_max:
        raise ValueError(
            "Please ensure that the third value in 'dl_bbox' (latmin) is less than the fourth value (latmax)."
        )
    return bbox


def _infer_product(dl_var: str) -> str:
    value = dl_var.upper()
    if value in {"SPM", "SPIM", "CHL", "CHLA"}:
        warnings.warn("No data product chosen, defaulting to 'SEXTANT' for SPM and Chl a data.")
        return "SEXTANT"
    if value in {"CDOM", "RRS", "NRRS", "T", "TUR", "SST"}:
        warnings.warn("No data product chosen, defaulting to 'ODATIS-MR'.")
        return "ODATIS-MR"
    raise ValueError("Variable not available, please check the value given for 'dl_var'")


def _get_sensor_dates(dl_product: str, dl_sensor: Optional[str], dl_correction: Optional[str]) -> Tuple[date, date]:
    if dl_product == "SEXTANT":
        return date(1998, 1, 1), date.today() - timedelta(days=7)
    if dl_sensor == "MODIS":
        return date(2002, 7, 4), date(2024, 12, 31)
    if dl_sensor == "MERIS":
        return date(2002, 6, 19), date(2012, 4, 8)
    if dl_sensor == "OLCI-A":
        return date(2016, 4, 26), date(2024, 12, 31)
    if dl_sensor == "OLCI-B":
        return date(2018, 5, 15), date(2024, 12, 31)
    raise ValueError(
        f"The chosen date range is outside of the available data range for: {dl_product} - {dl_sensor} - {dl_correction}."
    )


def _validate_product_settings(
    dl_var: str,
    dl_product: str,
    dl_sensor: Optional[str],
    dl_correction: Optional[str],
    dl_bbox: Optional[List[float]],
    username: Optional[str],
    password: Optional[str],
) -> Tuple[str, Optional[str], Optional[str]]:
    dl_var_upper = dl_var.upper()
    dl_product = dl_product.upper()

    if dl_product == "ODATIS-MR":
        if not username or not password:
            print(
                "ODATIS-MR data products require an AVISO+ account.\n"
                "Please provide your username and password via the 'username' and 'password' arguments.\n"
                "One may register for a free account here:\n"
                "https://www.aviso.altimetry.fr/en/data/data-access/registration-form.html"
            )
            return dl_product, None, None

        if dl_bbox is None:
            print("No bounding box provided, data will be downloaded for the full extent of the product.")
        elif (
            dl_bbox[0] < ODATIS_BBOX_LIMITS[0]
            or dl_bbox[1] > ODATIS_BBOX_LIMITS[1]
            or dl_bbox[2] < ODATIS_BBOX_LIMITS[2]
            or dl_bbox[3] > ODATIS_BBOX_LIMITS[3]
        ):
            print(
                "The bounding box provided exceeds the extent of the ODATIS-MR data product:\n"
                "lon: -7.8 to 10.3, lat: 41.2 to 51.5\n"
                "Data will be downloaded up to these limits."
            )
        else:
            print(
                f"Bounding box provided, data will be downloaded from longitude {dl_bbox[0]} to {dl_bbox[1]} "
                f"and latitude {dl_bbox[2]} to {dl_bbox[3]}."
            )

        if dl_sensor is None:
            if dl_var_upper == "SST":
                warnings.warn("No sensor chosen, defaulting to 'MODIS' because 'dl_var = SST'.")
                dl_sensor = "MODIS"
            else:
                warnings.warn("No sensor chosen, defaulting to 'OLCI-A'.")
                dl_sensor = "OLCI-A"
        elif dl_sensor not in {"MODIS", "MERIS", "OLCI-A", "OLCI-B"}:
            raise ValueError("Please set 'dl_sensor' to either 'MODIS', 'MERIS', 'OLCI-A', or 'OLCI-B'")

        if dl_correction is None:
            if dl_sensor == "MODIS":
                warnings.warn(
                    "No atmospheric correction chosen, defaulting to 'nirswir' to match 'dl_sensor = MODIS'."
                )
                dl_correction = "nirswir"
            else:
                warnings.warn(
                    f"No atmospheric correction chosen, defaulting to 'polymer' to match 'dl_sensor = {dl_sensor}'."
                )
                dl_correction = "polymer"
        elif dl_correction not in {"polymer", "nirswir"}:
            raise ValueError("Please set 'dl_correction' to either 'polymer' or 'nirswir'")

    if dl_product == "SEXTANT":
        if dl_bbox is not None:
            print(
                "A bounding box was provided, but SEXTANT products cannot be subset at the source. Downloading the full file."
            )
        if dl_var_upper not in {"SPM", "SPIM", "CHL", "CHLA"}:
            raise ValueError(
                "SEXTANT data product only contains SPM and Chl a data. Please adjust `dl_var` accordingly."
            )
        return dl_product, dl_sensor, dl_correction

    if dl_sensor == "MODIS":
        allowed_vars = {"CDOM", "CHL", "CHLA", "RRS", "NRRS", "SPM", "SPIM", "SST", "T", "TUR"}
        if dl_var_upper not in allowed_vars:
            raise ValueError(
                "ODATIS-MR MODIS data product does not contain the requested variable. Please adjust your variable choice accordingly."
            )
    elif dl_sensor in {"MERIS", "OLCI-A", "OLCI-B"}:
        allowed_vars = {"CDOM", "CHL", "CHLA", "RRS", "NRRS", "SPM", "SPIM", "T", "TUR"}
        if dl_var_upper not in allowed_vars:
            raise ValueError(
                f"ODATIS-MR: {dl_sensor} data product does not contain the requested variable. Please adjust your variable choice accordingly."
            )

    if dl_sensor == "MODIS" and dl_correction != "nirswir":
        print(
            "ODATIS-MR : MODIS data product only uses the 'nirswir' atmospheric correction. 'dl_correction' adjusted accordingly."
        )
        dl_correction = "nirswir"
    elif dl_sensor in {"MERIS", "OLCI-A", "OLCI-B"} and dl_correction != "polymer":
        print(
            f"ODATIS-MR : {dl_sensor} data product only uses the 'polymer' atmospheric correction. 'dl_correction' adjusted accordingly."
        )
        dl_correction = "polymer"

    return dl_product, dl_sensor, dl_correction


def _fetch_8day_start_date(
    current_date: date,
    dl_year: str,
    dl_month: str,
    dl_correction_flat: str,
    dl_sensor_flat: str,
) -> date:
    url = (
        "https://tds-odatis.aviso.altimetry.fr/thredds/catalog/"
        f"dataset-l3-ocean-color-odatis-mr-v_1_0.xml/FRANCE/{dl_correction_flat}/"
        f"{dl_sensor_flat}/8-day/{dl_year}/{dl_month}/catalog.html"
    )
    with urllib.request.urlopen(url) as response:
        html = response.read().decode("utf-8", errors="ignore")

    day_matches = sorted(set(re.findall(r'href="(\d{2})/catalog\.html"', html)))
    url_8day_dates = [date(int(dl_year), int(dl_month), int(day)) for day in day_matches]
    if len(url_8day_dates) < 3 or len(url_8day_dates) > 5:
        raise ValueError("Something has gone wrong with the 8-day HTML scraping.")

    return min(url_8day_dates, key=lambda value: abs((value - current_date).days))


def _build_odatis_var_chunk(dl_var: str, dl_sensor: str) -> str:
    value = dl_var.upper()
    if value == "CDOM":
        return "CDOM"
    if value in {"CHL", "CHLA"}:
        return "CHL-OC5"
    if value in {"RRS", "NRRS"}:
        return "NRRS555" if dl_sensor == "MODIS" else "NRRS560"
    if value in {"SPM", "SPIM"}:
        return "SPM-G"
    if value in {"T", "TUR"}:
        return "T-FNU"
    if value == "SST":
        return "SST-NIGHT"
    raise ValueError("Please check the value given for 'dl_var'")


def _copy_variable_subset(
    source: Dataset, 
    target: Dataset, 
    var_name: str, 
    lon_idx: np.ndarray, 
    lat_idx: np.ndarray
) -> None:
    src_var = source.variables[var_name]
    fill_value = getattr(src_var, "_FillValue", None)
    target_var = target.createVariable(var_name, src_var.datatype, src_var.dimensions, fill_value=fill_value)

    for attr_name in src_var.ncattrs():
        if attr_name == "_FillValue":
            continue
        target_var.setncattr(attr_name, src_var.getncattr(attr_name))

    slices = []
    for dim_name in src_var.dimensions:
        if dim_name == "lon":
            slices.append(slice(lon_idx[0], lon_idx[-1] + 1))
        elif dim_name == "lat":
            slices.append(slice(lat_idx[0], lat_idx[-1] + 1))
        else:
            slices.append(slice(None))
    target_var[:] = src_var[tuple(slices)]


def _subset_odatis_file(url: str, output_path: str, bbox: List[float]) -> None:
    with Dataset(url) as source:
        lon = np.asarray(source.variables["lon"][:])
        lat = np.asarray(source.variables["lat"][:])

        lon_idx = np.where((lon >= bbox[0]) & (lon <= bbox[1]))[0]
        lat_idx = np.where((lat >= bbox[2]) & (lat <= bbox[3]))[0]
        if lon_idx.size == 0 or lat_idx.size == 0:
            raise ValueError("The requested bounding box does not overlap the available data grid.")

        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with Dataset(output_path, "w", format="NETCDF4") as target:
            for dim_name, dim in source.dimensions.items():
                if dim_name == "lon":
                    target.createDimension("lon", len(lon_idx))
                elif dim_name == "lat":
                    target.createDimension("lat", len(lat_idx))
                else:
                    target.createDimension(dim_name, None if dim.isunlimited() else len(dim))

            for attr_name in source.ncattrs():
                target.setncattr(attr_name, source.getncattr(attr_name))

            target.setncattr("northernmost_latitude", float(lat[lat_idx].max()))
            target.setncattr("southernmost_latitude", float(lat[lat_idx].min()))
            target.setncattr("easternmost_longitude", float(lon[lon_idx].max()))
            target.setncattr("westernmost_longitude", float(lon[lon_idx].min()))
            target.setncattr(
                "spatial_subset",
                f"Spatially subset on {date.today().isoformat()} by "
                "https://github.com/RiOMar-projet/sat_access/blob/main/sat_access_script.py",
            )

            for var_name in source.variables:
                _copy_variable_subset(source, target, var_name, lon_idx, lat_idx)


def _build_auth_opener(username: Optional[str], password: Optional[str]) -> Optional[urllib.request.OpenerDirector]:
    if not username or not password:
        return None

    password_manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    password_manager.add_password(None, ODATIS_MR_BASE_URL.format(access_type="fileServer"), username, password)
    password_manager.add_password(None, ODATIS_MR_BASE_URL.format(access_type="dodsC"), username, password)
    auth_handler = urllib.request.HTTPBasicAuthHandler(password_manager)
    return urllib.request.build_opener(auth_handler)


def _download_file(
    url: str,
    output_path: str,
    username: Optional[str] = None,
    password: Optional[str] = None,
) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    opener = _build_auth_opener(username, password)
    if opener is None:
        urllib.request.urlretrieve(url, output_path)
    else:
        with opener.open(url) as response, open(output_path, "wb") as target:
            shutil.copyfileobj(response, target)
    print(f"File downloaded at: {output_path}")


def _decompress_bz2(path: str) -> None:
    output_path = path[:-4]
    with bz2.open(path, "rb") as src, open(output_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    print(f"File unzipped at: {output_path}")
    os.remove(path)
    print(f"File removed: {path}")


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
    dl_time_step = _check_time_step(dl_time_step)
    start_date, end_date = _check_date_range(dl_dates, dl_time_step)
    dl_bbox = _check_bbox(dl_bbox)
    dl_product = _infer_product(dl_var) if dl_product is None else dl_product.upper()
    dl_sensor = dl_sensor.upper() if dl_sensor is not None else None
    dl_correction = dl_correction.lower() if dl_correction is not None else None
    dl_product, dl_sensor, dl_correction = _validate_product_settings(
        dl_var=dl_var,
        dl_product=dl_product,
        dl_sensor=dl_sensor,
        dl_correction=dl_correction,
        dl_bbox=dl_bbox,
        username=username,
        password=password,
    )
    # if dl_product == "ODATIS-MR" and (dl_sensor is None or dl_correction is None):
    #     raise ValueError(
    #         "If downloading ODATIS-MR data, please provide a value for 'dl_sensor'."
    #     )

    floor_date_sensor, ceiling_date_sensor = _get_sensor_dates(dl_product, dl_sensor, dl_correction)
    if start_date < floor_date_sensor or end_date > ceiling_date_sensor:
        raise ValueError(
            f"The chosen date range is outside of the available data range for: {dl_product} - {dl_sensor} - {dl_correction}.\n"
            f"Available data range is from {floor_date_sensor} to {ceiling_date_sensor}. Please adjust your date range accordingly."
        )

    if dl_product == "SEXTANT" and dl_time_step != "day":
        print("SEXTANT data product only contains daily data. 'dl_time_step' adjusted accordingly.")
        dl_time_step = "day"
        start_date, end_date = _check_date_range(dl_dates, dl_time_step)

    current_date = start_date
    while current_date <= end_date:
        dl_date = current_date
        dl_date_flat = dl_date.strftime("%Y%m%d")
        dl_year = dl_date.strftime("%Y")
        dl_month = dl_date.strftime("%m")
        dl_day = dl_date.strftime("%d")
        url_year_doy = f"{dl_year}/{dl_date.strftime('%j')}"

        if dl_product == "SEXTANT":
            url_base = SEXTANT_BASE_URL
            if dl_var.upper() in {"SPM", "SPIM"}:
                url_product_stub = "EUR-L4-SPIM-ATL-v01"
            elif dl_var.upper() in {"CHL", "CHLA"}:
                url_product_stub = "EUR-L4-CHL-ATL-v01"
            else:
                raise ValueError("Variable not available")

            url_product = f"{url_product_stub}/{url_year_doy}"
            file_name = f"{dl_date_flat}-{url_product_stub}-fv01-OI.nc.bz2"
        elif dl_product == "ODATIS-MR":
            username_fix = quote(username or "", safe="")
            password_fix = quote(password or "", safe="")
            access_type = "dodsC" if dl_bbox is not None else "fileServer"
            url_base = ODATIS_MR_BASE_URL.format(access_type=access_type)

            dl_sensor_flat = dl_sensor.lower().replace("-", "")
            dl_sensor_chunk = dl_sensor.replace("CI-", "") if "OLCI" in dl_sensor else dl_sensor[:3]
            dl_correction_flat = dl_correction.lower().replace("-", "")
            dl_correction_chunk = "PO" if dl_correction == "polymer" else "NS"
            dl_var_chunk = _build_odatis_var_chunk(dl_var, dl_sensor)

            if dl_time_step == "day":
                dl_time_step_chunk = "DAY"
            elif dl_time_step == "8-day":
                dl_time_step_chunk = "8D"
                dl_date = _fetch_8day_start_date(
                    current_date=current_date,
                    dl_year=dl_year,
                    dl_month=dl_month,
                    dl_correction_flat=dl_correction_flat,
                    dl_sensor_flat=dl_sensor_flat,
                )
                dl_day = dl_date.strftime("%d")
                dl_date_flat = f"{dl_date.strftime('%Y%m%d')}-{(dl_date + timedelta(days=7)).strftime('%Y%m%d')}"
                current_date = dl_date

                if (dl_date + timedelta(days=7)).year > dl_date.year:
                    dl_date_flat = f"{dl_date.strftime('%Y%m%d')}-{date(dl_date.year, 12, 31).strftime('%Y%m%d')}"
                if dl_date + timedelta(days=7) > ceiling_date_sensor:
                    dl_date_flat = f"{dl_date.strftime('%Y%m%d')}-{ceiling_date_sensor.strftime('%Y%m%d')}"
            elif dl_time_step == "month":
                dl_time_step_chunk = "MO"
                dl_date_flat = (
                    f"{_month_start(dl_date).strftime('%Y%m%d')}-{_month_end(dl_date).strftime('%Y%m%d')}"
                )
            else:
                raise ValueError("Please check the value given for 'dl_time_step'")

            url_product = "/".join(
                [dl_correction_flat, dl_sensor_flat, dl_time_step, dl_year, dl_month, dl_day]
            )
            file_name = (
                f"L3m_{dl_date_flat}__FRANCE_03_{dl_sensor_chunk}_{dl_var_chunk}-{dl_correction_chunk}_"
                f"{dl_time_step_chunk}_00.nc"
            )
        else:
            raise ValueError("Please check the value used for 'dl_product'")

        url_final = f"{url_base}/{url_product}/{file_name}"
        file_path = os.path.join(output_dir, file_name)

        if os.path.exists(file_path.replace(".bz2", "")) and not overwrite:
            print(f"{file_path.replace(".bz2", "")} already exists. Set 'overwrite=True' to force the download.")
        else:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"File removed: {file_path}")

            try:
                if dl_bbox is not None and dl_product == "ODATIS-MR":
                    auth_url = url_final.replace("https://", f"https://{username_fix}:{password_fix}@", 1)
                    _subset_odatis_file(auth_url, file_path, dl_bbox)
                    print(f"Subsetting complete. The results are saved as: {file_path}")
                else:
                    _download_file(url_final, file_path, username=username, password=password)
                    if file_path.endswith(".bz2") and os.path.exists(file_path):
                        _decompress_bz2(file_path)
            except Exception as exc:
                message = str(exc)
                if "Server denied you to change to the given directory" in message:
                    message = "not found on server."
                print(f"{file_name}: {message}")

        if dl_time_step == "day":
            current_date += timedelta(days=1)
        elif dl_time_step == "8-day":
            current_date += timedelta(days=8)
        elif dl_time_step == "month":
            current_date = _next_month(current_date)


def _infer_nc_var_name(nc_file: str) -> str:
    file_snippets = os.path.basename(nc_file).split("_")
    if len(file_snippets) == 1:
        if "SPIM-ATL" in nc_file:
            return "analysed_spim"
        if "CHL-ATL" in nc_file:
            return "analysed_chl_a"
    elif len(file_snippets) == 9:
        return f"{file_snippets[6]}_mean"

    raise ValueError(
        "File structure cannot be inferred from file name. Please ensure the file was downloaded via the function found in this script."
    )


def _infer_var_label(nc_var_name: str) -> str:
    lower_name = nc_var_name.lower()
    if re.search(r"spim|spm", lower_name):
        return "SPM [g m-3]"
    if "chl" in lower_name:
        return "Chl a [mg m-3]"
    if "cdom" in lower_name:
        return "CDOM [m-1]"
    if "sst" in lower_name:
        return "SST [°C]"
    if "t-fnu" in lower_name:
        return "Turbidity [FNU]"
    if "nrrs560" in lower_name:
        return "Rrs at 560 nm [sr-1]"
    if "nrrs555" in lower_name:
        return "Rrs at 555 nm [sr-1]"
    raise ValueError("Variable label cannot be inferred from variable name.")


def _extract_plot_metadata(nc_data: Dataset, nc_file: str, nc_var_name: str) -> Tuple[str, date]:
    if "SPIM-ATL" in nc_file or "CHL-ATL" in nc_file:
        time_value = np.asarray(nc_data.variables["time"][:]).ravel()[0]
        plot_date = (datetime(1998, 1, 1) + timedelta(seconds=float(time_value))).date()
        return f"Map of {nc_var_name} on {plot_date}", plot_date

    if "L3m_" in nc_file:
        start_date = datetime.strptime(nc_data.getncattr("period_start_day"), "%Y%m%d").date()
        end_date = datetime.strptime(nc_data.getncattr("period_end_day"), "%Y%m%d").date()
        product_type = nc_data.getncattr("product_type")
        if product_type == "daily":
            return f"Map of {nc_var_name} on {start_date}", start_date
        return f"Map of {nc_var_name} from {start_date} to {end_date}", start_date

    raise ValueError("Date value cannot be inferred from file structure.")


def _subset_grid(
    lon: np.ndarray,
    lat: np.ndarray,
    var: np.ma.MaskedArray,
    bbox: List[float],
) -> Tuple[np.ndarray, np.ndarray, np.ma.MaskedArray]:
    lon_mask = (lon >= bbox[0]) & (lon <= bbox[1])
    lat_mask = (lat > bbox[2]) & (lat <= bbox[3])
    if not lon_mask.any() or not lat_mask.any():
        raise ValueError("The requested plot bounding box does not overlap the data grid.")

    lon_sub = lon[lon_mask]
    lat_sub = lat[lat_mask]
    var_sub = np.ma.asarray(var)[np.ix_(lat_mask, lon_mask)]
    return lon_sub, lat_sub, var_sub


def plot_nc(
    nc_file: str,
    bbox: Optional[List[float]] = None,
    plot_width: Optional[float] = None,
    plot_height: Optional[float] = None,
    output_dir: str = ".",
) -> None:
    if not os.path.exists(nc_file):
        print(f"{nc_file} cannot be found.")
        return

    nc_var_name = _infer_nc_var_name(nc_file)
    var_label = _infer_var_label(nc_var_name)

    with Dataset(nc_file) as nc_data:
        lon = np.asarray(nc_data.variables["lon"][:])
        lat = np.asarray(nc_data.variables["lat"][:])
        var = np.ma.asarray(nc_data.variables[nc_var_name][:]).squeeze()
        plot_title, plot_date = _extract_plot_metadata(nc_data, nc_file, nc_var_name)

    if bbox is None:
        bbox = [float(np.nanmin(lon)), float(np.nanmax(lon)), float(np.nanmin(lat)), float(np.nanmax(lat))]
        print(f"No bounding box provided, using full extent: {', '.join(f'{value:.4f}' for value in bbox)}")

    lon_sub, lat_sub, var_sub = _subset_grid(lon, lat, var, bbox)
    lon_grid, lat_grid = np.meshgrid(lon_sub, lat_sub)
    var_masked = np.ma.masked_invalid(np.ma.asarray(var_sub))
    var_masked = np.ma.masked_less(var_masked, 0)

    plot_width = 6 if plot_width is None else plot_width
    plot_height = 5 if plot_height is None else plot_height
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    plot_name = os.path.join(output_dir, f"{nc_var_name}_{plot_date}.png")

    fig, ax = plt.subplots(
        figsize=(plot_width, plot_height),
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    ax.set_extent(bbox, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="0.85")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.4)

    mesh = ax.pcolormesh(
        lon_grid,
        lat_grid,
        var_masked,
        shading="auto",
        cmap="viridis",
        transform=ccrs.PlateCarree(),
    )
    colorbar = fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9)
    colorbar.set_label(var_label)

    gridlines = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.5", alpha=0.5)
    gridlines.top_labels = False
    gridlines.right_labels = False

    ax.set_title(plot_title, fontsize=10 if " on " in plot_title else 8)
    fig.text(0.5, 0.02, f"Source: {nc_file}", ha="center", va="bottom", fontsize=6)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(plot_name, dpi=300)
    plt.close(fig)
    print(f"Image saved at: {plot_name}")

