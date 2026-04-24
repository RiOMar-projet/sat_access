# 🛰️ sat_access

> Simple `R` and `Python` tools for downloading and visualizing **SEXTANT** and **ODATIS-MR** satellite products for the RiOMar project.

## ✨ Overview

This repository contains lightweight scripts designed to make satellite data access straightforward for RiOMar users.

It supports:

- **SEXTANT** products
- **ODATIS-MR** products
- Both **`R`** and **`Python`** workflows
- Script-based use inside existing analysis pipelines
- Command-line access for direct downloads and plotting

> [!IMPORTANT]
> Access to **ODATIS-MR** data requires an **AVISO+** account.  
> You can create one for free here: [AVISO+ registration](https://www.aviso.altimetry.fr/en/data/data-access/registration-form.html)

## 🚀 Two Ways To Use This Repository

### 1. Script workflow

Use:

- `sat_access_script.R`
- `sat_access_script.py`

These files are intended to be imported into a normal workflow so you can reuse the download and plotting functions directly in your own scripts or notebooks.

Benefits:

- ✅ Little to no setup
- ✅ Easy to integrate into existing code
- ✅ Good choice for most users

Example notebook:

- [sat_access_script_examples.ipynb](https://github.com/RiOMar-projet/sat_access/blob/main/sat_access_script_examples.ipynb)

### 2. Command-line workflow

Use:

- `sat_access.R`
- `sat_access.py`

These scripts allow you to download **SEXTANT** NetCDF files for **Chl a** and/or **SPM**, and visualize a single day of data directly from the command line.

> [!NOTE]
> Running command-line requests in `R` may require more local setup than the script-based workflow.  
> For most `R` users, `sat_access_script.R` is the better starting point.

CLI usage notebook:

- [sat_access_cli_usage.ipynb](https://github.com/RiOMar-projet/sat_access/blob/main/sat_access_cli_usage.ipynb)

## 🌟 Features

- 📦 No installation required
- ⬇️ Download a script and use it directly
- 🌍 Download **SEXTANT** and **ODATIS-MR** NetCDF products
- 🗺️ Extract and visualize selected dates, bounding boxes, and variables
- 🐍 Available in `Python`
- 📘 Available in `R`

## ✅ Validation Notes

An initial validation of the satellite products available through these scripts has been performed for the RiOMar study regions.

Main takeaways:

- **ODATIS-MR MODIS** generally gives the closest matchups to *in situ* measurements from the **REPHY** and **SOMLIT** networks for both **SPM** and **Chl a**
- **SEXTANT** is often a close second, and sometimes the best option
- **MODIS** data have a **300 m** resolution
- **SEXTANT** data have a **1 km** resolution
- **SEXTANT** is **L4**, meaning gap-free fields
- **MODIS** is **L3**, meaning gaps may occur because of cloud cover and related issues

In practice:

- If you need finer spatial resolution, **MODIS** may be preferable
- If you want easier day-to-day handling with fewer gaps, **SEXTANT** is often the more practical choice
- **SEXTANT SPM** also compares well with *in situ* turbidity and may be used as a replacement where appropriate

## 🗂️ Project Updates

### 2026-03-30

- Created Jupyter notebooks with usage examples for both the `sat_access` and `sat_access_script` workflows in `R` and `Python`

### 2026-03-26

- Expanded `sat_access_script.py` to reach feature parity with `sat_access_script.R`

### 2026-03-03

- Created `sat_access_script.py`, initially with limited functionality
- Improved error trapping so a missing file on the server does not stop the full download process

### 2026-02-26

- Added username/password support for `ODATIS-MR` downloads

### 2026-02-12

- Added source-side subsetting for `ODATIS-MR` files via bounding box arguments

### 2026-01-22

- Added support for downloading all 8-day and monthly `ODATIS-MR` products via `sat_access_script.R`
- Added plotting support for all daily, 8-day, and monthly data files downloadable via `sat_access_script.R`
- Created this **Updates** section and backfilled project history

### 2026-01-21

- Added support for downloading all daily `ODATIS-MR` products via `sat_access_script.R`

### 2025-01-19

- Updated the `README` to document the addition of `sat_access_script.R`

### 2025-01-16

- Created `sat_access_script.R` based on user feedback as a more traditional non-CLI `R` workflow

### 2025-01-15

- Improved logic for recognizing already-downloaded files in `sat_access.R` and `sat_access.py`
- Improved logic for unzipping files in `sat_access.R`

### 2025-12-17

- Extended `sat_access.py` to support concatenating daily files and creating GIFs
- Updated documentation to match new arguments in `sat_access.py`
- Improved logic gates for handling various requests

### 2025-12-16

- Added a **TO DO** section

### 2025-11-26

- Started the `sat_access` project

## 📄 License

This project is provided as-is under the [MIT License](https://github.com/RiOMar-projet/sat_access/blob/main/LICENSE).  
No liability is assumed and no warranty is provided.

## 📬 Support

For questions or issues, contact Robert Schlegel:

- [robert.schlegel@imev-mer.fr](mailto:robert.schlegel@imev-mer.fr)
