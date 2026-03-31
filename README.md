# Overview

This repository contains scripts designed to facilitate access to the **SEXTANT** and **ODATIS-MR** satellite products for members of the RiOMar project. **NB**: Access to **ODATIS-MR** data requires an AVISO+ account. This can be created for free [here](https://www.aviso.altimetry.fr/en/data/data-access/registration-form.html).

The `sat_access_script.R` & `sat_access_script.py` scripts have been created so that they can be used to import the necessary acces/plotting functions into whatever standard workflow the user would be developing that would need access to satellite data. The scripts may simply be downloaded individually as the user needs and run on their local computer. They are designed to require little to no setup in order to download and visualise the data files. The Jupyter notebook giving examples for their usage may be found [here](https://github.com/RiOMar-projet/sat_access/blob/main/sat_access_script_examples.ipynb).

There are also the `sat_access.py` & `sat_access.R` scripts, which allow the user to download **SEXTANT** NetCDF files for Chl _a_ and/or SPM from an FTP server and visualize one day of data on a map directly via command line interface (CLI) calls. Note however that running command line requests for `R` can require a non-neglible amount of setup on one's local computer, so it is advised to use `sat_access_script.R`. The Jupyter notebook that explains how to use these two CLI focussed scripts may be found [here](https://github.com/RiOMar-projet/sat_access/blob/main/sat_access_cli_usage.ipynb).

# Features

-   No installation required
-   Simply download a script and source it or call it from the command-line
-   Downloads a range of NetCDF files for **SEXTANT** and **ODATIS-MR**
-   Extracts and visualizes specified dates, bounding boxes, and variables as a map
-   Available for both `Python` and `R`

# Validation

Note that an initial validation of all of the satellite products available via these scripts has been performed for the regions of study for RiOMar. Generally, for both SPM and CHl _a_, the data from the **ODATIS-MR** _MODIS_ sensor provide the closest matchups to _in situ_ measurements made for both the REPHY and SOMLIT networks. A close second (sometimes first) are the data from **SEXTANT** (one merged 'sensor'). Note that the data for MODIS have a 300 m resolution, while **SEXTANT** is 1 km. **SEXTANT** however is L4 (i.e. no gaps in the data), while MODIS is L3 (i.e. gaps due to cloud cover etc.). Meaning that, unless one has need for a high resolution analysis, the **SEXTANT** data will generally be easier to work with, while still providing an accurate estimate for SPM and Chl _a_ in the coastal zone. Also note that **SEXTANT** SPM also compares well against _in situ_ measured turbidity so may be used as a replacement if necessary.

# Updates

## 2026-03-30

-   Created Jupyter notebooks to provide examples for `R` and `Python` usage for the `sat_access` and `sat_access_script` workflows.

## 2026-03-26

-   Expanded `sat_access_script.py` to have the same functionality as `sat_access_script.R`.

## 2026-03-03

-   Created `sat_access_script.py`, which is similar to `sat_access_script.R` but with limited functionality.

## 2026-03-03

-   Improved error trapping so a missing file on the server doesn't stop the entire download process.

## 2026-02-26

-   Added username/password functionality for `ODATIS-MR` downloads.

## 2026-02-12

-   Added ability to subset `ODATIS-MR` files at the source via bounding box argument.

## 2026-01-22

-   Added the ability to download all of the 8-day and monthly products within `ODATIS-MR` via `sat_access_script.R`.
-   Added the ability to plot all of the daily, 8-day, and monthly data files that can be downloaded via `sat_access_script.R`.
-   Created the __Updates__ section and populated it with the history of the project to date.

## 2026-01-21

-   Added the ability to download all of the daily products within `ODATIS-MR` via `sat_access_script.R`.

## 2025-01-19

-   Updated __README__ documentation to cover the addition of `sat_access_script.R` to the project.

## 2025-01-16

- Based on user feedback, created the `sat_access_script.R`, which is run like a 'traditional' R script and does not use the command line.

## 2025-01-15

-   Improved logic gates around the recognition of already downloaded files in `sat_access.R` and `sat_access.py`.

## 2025-01-15

-   Improved logic gates around the unzipping of files in `sat_access.R`.

## 2025-12-17

-   Extended functionality of `sat_access.py` to be able to concatenate daily files and create GIFs.
-   Updated documentation to match new arguments for `sat_access.py`.
-   Improved logic gates for handling various requests.

## 2025-12-16

-   Added __TO DO__ section.

## 2025-11-26

-   Started `sat_access` project.

# License

This script is provided as-is and is free to use via an [MIT License](https://github.com/RiOMar-projet/sat_access/blob/main/LICENSE).
No liability is assumed and no warranty is provided.

# Support

For questions or issues, please contact Robert Schlegel at: [robert.schlegel\@imev-mer.fr](mailto:robert.schlegel@imev-mer.fr){.email}

