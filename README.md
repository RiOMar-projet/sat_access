# Overview

The `sat_access.py` & `sat_access.R` scripts allow you to download a series NetCDF files for chl a and/or SPM from an FTP server and visualize one day of data on a map.
They are designed specifically to work for the suite of ODATIS products put into production for/during the RiOMar project.
Note that while these scripts have been written so that they can be called from the command line, 
they can also be opened in your IDE of choice (e.g. `VS Code`, `RStudio`, etc.) and the parameters set directly in the code.
Indeed, command line based operations for the `R` language have proven to be more complicated than it may be worth. So a third script,
`sat_access_script.R` has been created which is designed to run like a 'traditional' `R` script.
Meaning it can be downloaded from this repository, and then opened from an IDE (e.g. `RStudio`) and the functions can be run from there.

# Features

-   No installation required, simply download a script and call it from the command-line
-   Uses **`argparse`** for easy command-line argument handling
-   Downloads a range of NetCDF file from an FTP server
-   Extracts and visualizes specified dates, bounding boxes, and variables as a map
-   Available for both `Python` and `R`

# How It Works

1.  FTP Download:

-   The script connects to the specified FTP server based on the requested variable and date range and downloads the NetCDF file to your local machine.

2.  Visualisation:

-   It loads the NetCDF file and clips the data to the specified longitude, latitude, and the chosen variable before plotting it as a map using **`matplotlib`** (`Python`) or **`ggplot2`** (`R`).

3.  Saving plots:

-   The `R` script creates and saves the plot directly to the specified `outputdir` folder, while the `Python` script displays the plot in an interactive window.

# Requirements

## Python

-   Python: 3.x
-   Packages:
    -   netCDF4
    -   xarray
    -   matplotlib
    -   cartopy

To check if you have Python installed, open a terminal and type:

``` bash
python --version
```

Or if you would like to use an existing virtual environment (recommended):

``` bash
conda activate your_env_name
```

If your machine cannot find an active version of Python, but you should have one, see the troubleshooting below.

Otherwise, once you have an active Python environment available in your terminal, install the required packages from the terminal with:

``` bash
pip install xarray matplotlib cartopy
```

Note that if you have not used `cartopy` before, the first time you run the `sat_access.py` script it will download a few shape files.
This may take a few minutes on a slow internet connection.
During which time the map menu will appear to be hanging.
You will know it is working if you see activity in the console as it downloads the necessary shapefiles.

## R

**NB:** The following requirements are only for `sat_access.R`, these requirements can be ignored if one chooses to use `sat_access_script.R` instead.

-   R: 4.x
-   Packages:
    -   argparse
    -   ncdf4
    -   curl
    -   lubridate
    -   reshape2
    -   ggplot2

To check if you have R installed, open a terminal and type:

(Linux/MacOS)

``` bash
R --version
```

(Windows - PowerShell)

``` powershell
Rscript --version
```

If your machine cannot find an active version of R, but you should have one, see the troubleshooting section below.

(Linux/MacOS/Windows) Otherwise, install the required packages from the terminal/PowerShell with:

``` bash
Rscript -e "install.packages(c('argparse', 'ncdf4', 'curl', 'ggplot2', 'reshape2'), repos='https://cloud.r-project.org/')"
```

### Linux/MacOS

Finally, make the script executable (run from the same location as the script):

``` bash
chmod +x sat_access.R
```

# Usage

## Command-Line Arguments

| Argument | Description | Required | Example Value |
|------------------|------------------|------------------|------------------|
| `--variable` | The data layer to download | Yes | SPM |
| `--daterange` | Date range for desired data in YYYY-MM-DD format | Yes | 2025-10-15 |
| `--outputdir` | Local folder to save the NetCDF file | Yes | downloads |
| `--overwrite` | Overwrite existing files; Default = False | No | True |
| `--concatenate` | Whether or not to concatenate the downloaded data into a single file; Default = False | No | True |
| `--plot` | Whether or not to plot the downloaded data; Default = False | No | True |
| `--animate` | Whether or not to animated the downloaded data; Default = False | No | True |
| `--boundingbox` | Bounding box as: lon_min lon_max lat_min lat_max | No | 4 6 42 44 |

__NB:__ The `--concatenate` and `--animate` arguments are currently only available in the `Python` version.

## Example Usage

The following examples assume the user is in a terminal at the same location as the script.
Note that the argument `--outputdir .` specifies the current directory as the output folder.

### Python

Download a single chl a file:

``` bash
python sat_access.py --variable chla --date 2025-10-15 --outputdir .
```

Download multiple SPM files and plot one:

``` bash
python sat_access.py --variable spm --date 2025-10-15 2025-10-17 --outputdir . --plot True --boundingbox 4 6 42 44
```

Download multiple SPM files, concatenate them into a single file, and create an animation:

``` bash
python sat_access.py --variable spm --date 2025-10-01 2025-10-31 --outputdir . --concatenate True --animate True --boundingbox 4 6 42 44
```

### R

**NB:** The following command line examples only work for `sat_access.R`, not `sat_access_script.R`.

#### Linux/MacOS

Download a single chl a file and plot it:

``` bash
./sat_access.R --variable chla --date 2025-09-15 --outputdir . --plot TRUE --boundingbox 4 6 42 44
```

Download multiple SPM files but plot none:

``` bash
./sat_access.R --variable spm --date 2025-11-15 2025-11-17 --outputdir . --plot FALSE
```

**NB:** The `./` before `sat_access.R` is necessary for bash to understand that this R script is meant to be run as an executable.

#### Windows

Download a single chl a file and plot it:

``` powershell
Rscript sat_access.R --variable chla --date 2025-09-15 --outputdir . --plot TRUE --boundingbox 4 6 42 44
```

# Troubleshooting

-   FTP Connection Issues: Check your internet connection and firewall settings. Or wait a few seconds and try again.
-   Missing Variables: If not all required variables have been imputed in the terminal, an automated error message should be generated to help you along.
-   Module Errors: Ensure all required `Python` or `R` packages are installed.
-   `Python` or `R` Not Found: See additional steps below.

## Windows

### Python

#### Error: Python not found

If this returns nothing:

```         
python --version
```

But you are certain Python is installed (e.g. via anaconda), you may need to add Python to your PATH.

1.  Open System Environment Variables:

-   Press Win + S, type "Environment Variables", and select "Edit the system environment variables".
-   Click "Environment Variables".

2.  Edit the PATH Variable:

-   Under "User variables" or "System variables", find the Path variable and click "Edit".
-   Add the following paths (adjust for your installation):
-   `C:\Users\YourUsername\anaconda3`
-   `C:\Users\YourUsername\anaconda3\Scripts`
-   `C:\Users\YourUsername\anaconda3\Library\bin`
-   Click OK to save.

3.  Restart PowerShell:

-   Close and reopen PowerShell for the changes to take effect.

Alternatively, if you already have a virtual environment installed on your system it would be preferable to activate it before running the script.

To check the virtual environments installed on your system (with anaconda/miniconda):

```         
conda env list
```

To activate the environment of choice:

```         
conda activate your_env_name
```

#### Error: 'ScipyArrayWrapper' object has no attribute 'oindex'

This error is related to the NetCDF file not being loaded correctly by `xarray`.
To address this issue one should update both `xarray` and `netCDF4` to the latest versions.

```         
pip install --upgrade xarray netCDF4
```

### R

**NB:** If you are encountering issues trying to get `sat_access.R` to work, it is strongly advised to use `sat_access_script.R` instead.
This is because there can be a very long list of potential issues getting the command line functionality to work with `R`. 
Particularly on computers running a Windows based operating systems.

#### Error: R not found

If this returns nothing:

```         
Rscript --version
```

But you are certain R is installed, you may need to add R to your PATH.
This is done the same as for Python above, but adding the path to your R installation, e.g.: `C:\Program Files\R\R-4.x.x\bin\` **NB:** Replace `x.x` with your installed version number.
Restart the PowerShell and try again.

If this still doesn't work, check that the command `R` hasn't already been mapped to something else by default: In a PowerShell terminal type:

```         
Get-Alias -Name R
```

If this returns an Alias that isn't `R.exe` you will need to add a new alias for R: (**NB:** Replace `x.x` with your installed version number.)

```         
Add-Content -Path $PROFILE -Value "`nNew-Alias -Name R -Value 'C:\Program Files\R\R-4.x.x\bin\R.exe'"
```

Restart PowerShell and heck that it worked:

```         
Get-Alias -Name Rcode
```

You should see that the Alias for `Rcode` is now mapped to your `R` installation (e.g. `R.exe`).
If yes, check the `R` version:

```         
Rcode --version
```

Once `R` is recognised by your system, don't forget to install the necessary packages

```         
Rscript -e "install.packages(c('argparse', 'ncdf4', 'curl', 'ggplot2', 'reshape2'), repos='https://cloud.r-project.org/')"
```

**NB:** Regardless of what new alias you may have created for `R`, the installation of packages is still done from the PowerShell with the command `Rscript`.
Or, you can also simply open `R`/`RStudio` and install the packages from there.
Whatever is easiest for you.

#### Error: Couldn't find sufficient Python binary

This error is caused by the system not understanding how to find `Python` while calling scripts from `R`.
To address this issue open the `sat_access.R` script in your IDE of choice (e.g. `RStudio`) and uncomment the following line of code and change the file pathway to where you have `Python` installed:

```         
options(python_cmd = "C:/Path/To/Your/Python/python.exe")
```

# TO DO

-   Add knowledge of which algorithms produce best results by region
-   Allow function to automatically choose the best processing algorithm for the given request
-   Add `ODATIS-MR` downloading functionality for 8-day and monthly data to all scripts
-   Add `ODATIS-MR` daily downloading functionality to `sat_access.py` and `sat_access.R`
-   Add `ODATIS-MR` plotting functionality to `sat_access.py` and `sat_access.R`

# Updates

## 2026-01-22

-   Added the ability to plot all of the daily data files that can be downloaded via `sat_access_script.R`.
-   Created the __Updates__ section and populated it with the history of the project to date.

## 2026-01-21

-   Added the ability to download all of the different products within `ODATIS-MR` to `sat_access_script.R`.

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

