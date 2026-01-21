# sat_access_script.R

# Title: Satellite Data Access and Visualization Script
# Author: Robert Schlegel
# Date: January 2026
# Description: This script contains functions that can be used to
#              download SPM or chlorophyll-a NetCDF files from a specified URL
#              for given date ranges and plot the data.


# Libraries ---------------------------------------------------------------

# Check for missing libraries and install them if necessary
if (!all(c("ncdf4", "curl", "reshape2", "ggplot2") %in% installed.packages())) {
  install.packages(c("ncdf4", "curl", "lubridate", "reshape2", "ggplot2"), repos = "https://cloud.r-project.org/")
}

# Activate libraries
library(ncdf4)    # For reading NetCDF files
suppressPackageStartupMessages(library(curl)) # For FTP download
library(lubridate) # For working with dates
library(reshape2) # For data reshaping
library(ggplot2)  # For visualization


# The download function ---------------------------------------------------

# testing...
# dl_var = "CDOM"
# dl_dates = c("2022-09-01", "2022-09-05")
# dl_product = "ODATIS-MR"
# dl_sensor = "OLCI-A"
# dl_correction = "polymer"
# dl_time_step = "day"
# output_dir = "~/Downloads"
# overwrite = FALSE 

download_nc <- function(dl_var, dl_dates, 
                        dl_product = NULL, dl_sensor = NULL, 
                        dl_correction = NULL, dl_time_step = NULL,
                        output_dir, overwrite) {
  
  # Check date range
  if(length(dl_dates) > 2){
    stop("Please provide only one or two dates for the date range.")
  }
  
  # Check requested time-step
  if(is.null(dl_time_step)){
    message("No time-step chosen, defaulting to daily data.")
    dl_time_step <- "day"
  } else if(dl_time_step %in% c("day", "daily")){
    dl_time_step <- "day"
  } else if(dl_time_step %in% c("month", "monthly")){
    dl_time_step <- "month"
  } else if(dl_time_step %in% c("8d", "8D", "8day", "8Day", "8-day", "8-daily", "weekly")){
    dl_time_step <- "8-day"
  } else {
    stop("'dl_time_step value not recognised, please choose one of : 'day', '8-day', 'month'")
  }
  
  # Set start and end dates for download
  # TODO: Correct this based on the time-step
  # Need special setup for 8-day and monthly requests
  if(dl_time_step == "day"){
    start_date <- as.Date(dl_dates[1])
    if(length(dl_dates) == 2){
      end_date <- as.Date(dl_dates[2])
    } else{
      end_date <- start_date
    }
  } else if(dl_time_step == "8-day"){
    stop("8-day time-steps not yet implemented")
  } else if(dl_time_step == "month"){
    stop("month time-steps not yet implemented")
  } else {
    stop("'dl_time_step value not recognised, please choose one of : 'day', '8-day', 'month'")
  }
  
  # Check chosen variable against chosen data product
  if(is.null(dl_product)){
    if(toupper(dl_var) %in% c("SPM", "SPIM", "CHL", "CHLA")){
      message("No data product chosen, defaulting to 'SEXTANT' for SPM and Chl a data.")
      dl_product <- "SEXTANT"
    } else if(toupper(dl_var) %in% c("CDOM", "RRS", "NRRS", "T", "TUR", "SST")){
      message("No data product chosen, defaulting 'ODATIS-MR'.")
      dl_product <- "ODATIS-MR"
    } else {
      stop("Variable not available, please check the value given for 'dl_var'")
    }
  }
  
  # Get sensors and corrections for ODATIS-MR data products
  if(dl_product == "ODATIS-MR"){
    
    # Check chosen sensor against data product
    if(is.null(dl_sensor)){
      if(dl_var == "SST"){
        message("No sensor chosen, defaulting to 'MODIS' because 'dl_var = SST'.")
        dl_sensor <- "MODIS"
      } else {
        message("No sensor chosen, defaulting to 'OLCI-A'.")
        dl_sensor <- "OLCI-A"
      }
    } else if(!dl_sensor %in% c("MODIS", "MERIS", "OLCI-A", "OLCI-B")){
      stop("Please set 'dl_sensor' to either 'MODIS', 'MERIS', 'OLCI-A', or 'OLCI-B'")
    }
    
    # Check chosen correction against data product
    if(is.null(dl_correction)){
      if(dl_sensor == "MODIS"){
        message("No atmospheric correction chosen, defaulting to 'nirswir' to match 'dl_sensor = MODIS'.")
        dl_correction <- "nirswir"
      } else if(dl_sensor %in% c("MERIS", "OLCI-A", "OLCI-B")){
        message(paste0("No atmospheric correction chosen, defaulting to 'polymer' to match 'dl_sensor = ",dl_sensor,"'.."))
        dl_correction <- "polymer"
      }
    } else if(!dl_correction %in% c("polymer", "nirswir")){
      stop("Please set 'dl_correction' to either 'polymer' or 'nirswir'")
    }
  }
  
  # Check chosen dates against chosen data product, sensor, and correction
  # TODO: This could be optimized
  if(toupper(dl_product) == "SEXTANT"){
    if(start_date < as.Date("1998-01-01") | end_date > as.Date(Sys.Date())-7){
      stop("SEXTANT data are only available from 1998-01-01 to roughly 1 week before the present date. Please adjust your date range accordingly.")
    }
  } else if(dl_product == "ODATIS-MR"){
    if(dl_sensor == "MODIS"){
      # TODO: Increase data ranges when the data are uploaded through to 2026
      if(start_date < as.Date("2002-07-04") | end_date > as.Date("2023-12-31")){
        stop(paste0("ODATIS-MR MODIS data are only available from 2002-07-04 to 2023-12-31. Please adjust your date range accordingly."))
      }
    } else if(dl_sensor == "MERIS"){
      if(start_date < as.Date("2002-06-19") | end_date > as.Date("2012-04-08")){
        stop(paste0("ODATIS-MR MERIS data are only available from 2013-07-01 to . Please adjust your date range accordingly."))
      }
    } else if(dl_sensor == "OLCI-A"){
      if(start_date < as.Date("2016-04-26") | end_date > as.Date("2023-12-31")){
        stop(paste0("ODATIS-MR OLCI-A data are only available from 2016-04-26 to 2023-12-31. Please adjust your date range accordingly."))
      }
    } else if(dl_sensor == "OLCI-B"){
      if(start_date < as.Date("2018-05-15") | end_date > as.Date("2023-12-31")){
        stop(paste0("ODATIS-MR OLCI-B data area only available from 2018-05-15 to 2023-12-31. Please adjust your date range accordingly."))
      }
    }
  }
  
  # TODO: 
  # Create a more intelligent logic gate that chooses the best product, sensor, and correction based on the variable and download dates chosen
  # It must create messages for the user letting them know what has been chosen
  # This could also be made to favour higher resolution data
  # The messages could briefly explain why the choice was made
  # This would address the TODO item of choosing the 'best' chl-a etc. data given the bbox etc.
  
  # Iterate over each date in the range
  # TODO: This will need to be reworked if the user is requesting 8-day or monthly data
  current_date <- start_date
  while(current_date <= end_date){
    
    # Prep date strings
    dl_date <- current_date
    dl_date_flat <- gsub("-", "", dl_date)
    dl_year <- substr(dl_date, start = 1, stop = 4)
    dl_month <- substr(dl_date, start = 6, stop = 7)
    dl_day <- substr(dl_date, start = 9, stop = 10)
    url_year_doy <- paste0(substr(dl_date, start = 1, stop = 4),"/",strftime(as.Date(dl_date), format = "%j"))
    
    # Get URL specifics for SEXTANT
    if(toupper(dl_product) == "SEXTANT"){
      
      # Base URL
      url_base <- "ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic"
      
      # Prep file stub string
      if(toupper(dl_var) %in% c("SPM", "SPIM")){
        url_product_stub <- "EUR-L4-SPIM-ATL-v01"
      } else if(toupper(dl_var) %in% c("CHL", "CHLA")){
        url_product_stub <- "EUR-L4-CHL-ATL-v01"
      } else {
        stop("Variable not available")
      }
      
      # Variable specifics
      url_product <- paste0(url_product_stub,"/",url_year_doy)
      file_name <- paste0(dl_date_flat,"-",url_product_stub,"-fv01-OI.nc.bz2")
      nc_file <- file.path(output_dir, paste0(dl_date_flat,"-",url_product_stub,"-fv01-OI.nc"))

    # Get URL specifics for ODATIS-MR
    } else if(toupper(dl_product) == "ODATIS-MR"){
      
      # Base URL
      url_base <- "https://tds%40odatis-ocean.fr:odatis@tds-odatis.aviso.altimetry.fr/thredds/fileServer/dataset-l3-ocean-color-odatis-mr-v_1_0.xml/FRANCE"
      
      # Prep sensor strings
      dl_sensor_flat <- tolower(gsub("-", "", dl_sensor))
      if(dl_sensor == "OLCI-A"){
        dl_sensor_chunk <- "OLA"
      } else if(dl_sensor == "OLCI-B"){
        dl_sensor_chunk <- "OLB"
      } else if(dl_sensor == "MERIS"){
        dl_sensor_chunk <- "MER"
      } else if(dl_sensor == "MODIS"){
        dl_sensor_chunk <- "MOD"
      } else {
        stop("Please check the value given for 'dl_sensor'")
      }
      
      # Prep correction strings
      dl_correction_flat <- tolower(gsub("-", "", dl_correction))
      if(dl_correction == "polymer"){
        dl_correction_chunk <- "PO"
      } else if(dl_correction == "nirswir"){
        dl_correction_chunk <- "NS"
      } else {
        stop("Please check the value given for 'dl_correction'")
      }
      
      # Prep time-step string
      if(dl_time_step == "day"){
        dl_time_step_chunk <- "DAY"
      } else if(dl_time_step == "8-day"){
        # TODO: Need to think about how to manage the 8 day jumps
        # Would be ideal to be able to quickly access the HTML structure within the chosen years+months of data
        dl_time_step_chunk <- "8D"
      } else if(dl_time_step == "month"){
        dl_time_step_chunk <- "MO"
        dl_date_flat <- paste0(format(floor_date(ymd(dl_date), "month"), "%Y%m%d"), "-", format(ceiling_date(ymd(dl_date), "month") - days(1), "%Y%m%d"))
      } else {
        stop("Please check the value given for 'dl_time_step'")
      }
      
      # Prep the variable string
      if(dl_var %in% c("CDOM")){
        dl_var_chunk <- dl_var
      } else if(dl_var %in% c("CHL", "CHLA")){
        dl_var_chunk <- "CHL-OC5"
      } else if(dl_var %in% c("RRS", "NRRS")){
        if(dl_sensor == "MODIS"){
          dl_var_chunk <- "NRRS560"
        } else {
          dl_var_chunk <- "NRRS555"
        }
      } else if(dl_var %in% c("SPM", "SPIM")){
        dl_var_chunk <- "SPM-G"
      } else if(dl_var %in% c("T", "TUR")){
        dl_var_chunk <- "T-FNU"
      } else if(dl_var %in% c("SST")){
        dl_var_chunk <- "SST-NIGHT"
      } else {
        stop("Please check the value given for 'dl_sensor'")
      }
      
      # Product URL
      url_product <- paste(dl_correction_flat, dl_sensor_flat, dl_time_step, dl_year, dl_month, dl_day, sep = "/")
      file_name <- paste0("L3m_",dl_date_flat,"__FRANCE_03_",dl_sensor_chunk,"_",dl_var_chunk,"-",dl_correction_chunk,"_",toupper(dl_time_step),"_00.nc")
      nc_file <- file.path(output_dir, file_name)
      
    } else {
      stop("Please check the value used for 'dl_product'")
    }

    # Assemble final URL
    url_final <- paste(url_base, url_product, file_name, sep = "/")
    file_name_full <- file.path(output_dir, file_name)
    
    # Fetch file
    if(file.exists(nc_file) & !overwrite){
      
      message(paste0(nc_file," already exists. Set 'overwrite = TRUE' to force the download."))
      
    } else {
      
      # Ensure output directory exists
      if(!dir.exists(output_dir)){
        dir.create(output_dir, recursive = TRUE)
      }
      
      # Download
      tryCatch({
        curl::curl_download(url_final, destfile = file_name_full)
        message(paste0("File downloaded at: ",file_name_full))
      }, error = function(e) {
        if(grepl("Server denied you to change to the given directory", e$message[1])){
          e$message <- " not found on server."
        }
        message(paste0(file_name,e$message[1]))
      })
      
      # Unzip if necessary
      if(file.exists(file_name_full)){
        if(grepl("bz2", file_name_full)){
          result <- system(paste("bunzip2 -k -f", file_name_full), intern = TRUE, ignore.stderr = FALSE)
          if(length(result) != 0){
            message("Failed to unzip the file: ", file_name_full)
          } else {
            message("File unzipped at: ", gsub(".bz2","",file_name_full))
            tryCatch({
              file.remove(file_name_full)
              message(paste("File removed :", file_name_full))
            }, error = function(e) {
              message(paste("Impossible to remove file :", e$message))
            })
          }
        }
      }
    }
    
    # Move to the next day
    # TODO: Change this to match the requested time-step
    current_date <- current_date + 1
  }
}
#

# The plotting function ---------------------------------------------------

plot_nc <- function(nc_file, bbox = NULL, plot_width = NULL, plot_height = NULL, output_dir) {
  
  # Check for file
  if(!file.exists(nc_file)){
    return(message(paste0(nc_file," cannot be found.")))
  }
  
  # Get product specifics
  if(grepl("SPIM", nc_file)){
    nc_var_name <- "analysed_spim"
    var_label <- "SPM [g m-3]"
  } else if(grepl("CHL", nc_file)){
    nc_var_name <- "analysed_chl_a"
    var_label <- "chl a [mg m-3]"
  } else {
    stop("Variable not yet available")
  }
  
  # Open the NetCDF file
  nc_data <- nc_open(nc_file)
  
  # Extract longitude, latitude, and the specified variable
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat")
  time <- ncvar_get(nc_data, "time")
  var <- ncvar_get(nc_data, nc_var_name)
  
  # Close the NetCDF file
  nc_close(nc_data)
  
  # Reshape data for ggplot
  var_df <- melt(var)
  names(var_df) <- c("lon_idx", "lat_idx", "value")
  var_df$lon <- lon[var_df$lon_idx]
  var_df$lat <- lat[var_df$lat_idx]
  
  if(is.null(bbox)){
    bbox <- c(
      min(lon, na.rm = TRUE),
      max(lon, na.rm = TRUE),
      min(lat, na.rm = TRUE),
      max(lat, na.rm = TRUE)
    )
    message(paste("No bounding box provided, using full extent:", paste(round(bbox, 4), collapse = ", ")))
  }
  
  # Filter df to bounding box
  var_df_sub <- var_df[var_df$lon >= bbox[1] & var_df$lon <= bbox[2],]
  var_df_sub <- var_df_sub[var_df_sub$lat > bbox[3] & var_df_sub$lat <= bbox[4], ]
  
  # Get date for plot label
  plot_date <- as.Date(as.POSIXct(time, origin = "1998-01-01"))
  
  # Set plot name
  plot_name <- file.path(output_dir, paste0(nc_var_name,"_",plot_date,".png"))
  
  # Plot using ggplot2
  p <- ggplot(var_df_sub, aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(var_label) +
    labs(title = paste("Map of", nc_var_name, "on", plot_date),
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Determine plot dimensions, save, and exit
  if(is.null(plot_width)) plot_width <- 6
  if(is.null(plot_height)) plot_height <- 5
  ggsave(filename = plot_name, plot = p, width = plot_width, height = plot_height)
  message(paste0("Image saved at: ",plot_name))
}


# Examples ----------------------------------------------------------------

## Downloading ------------------------------------------------------------

# Download a few days of SPM data
download_nc(
  dl_var = "SPM",
  dl_dates = c("2025-09-01", "2025-09-05"),
  output_dir = "~/Downloads", # Change as desired/required
  overwrite = FALSE # Change to TRUE to force downloads
)

# Download one day of Chl a data
download_nc(
  dl_var = "CHL",
  dl_dates = c("2025-12-25"),
  output_dir = "~/Downloads", # Change as desired/required
  overwrite = FALSE # Change to TRUE to force downloads
)

# Download a few days day of CDOM data
download_nc(
  dl_var = "CDOM",
  dl_dates = c("2022-09-01", "2022-09-05"),
  output_dir = "~/Downloads",
  overwrite = FALSE
)

# Download one day day of SST data
download_nc(
  dl_var = "SST",
  dl_dates = "2014-01-01",
  output_dir = "~/Downloads",
  overwrite = FALSE
)


## Plotting ---------------------------------------------------------------

# Plot an SPM NetCDF file
plot_nc(
  nc_file = "~/Downloads/20250905-EUR-L4-SPIM-ATL-v01-fv01-OI.nc", 
  bbox =  c(-5, 5, 35, 45), 
  plot_width = 6, 
  plot_height = 5, 
  output_dir = "~/Downloads"
)

# Plot a Chl a NetCDF file
plot_nc(
  nc_file = "~/Downloads/20251225-EUR-L4-CHL-ATL-v01-fv01-OI.nc", 
  bbox =  c(-2.5, -0.5, 44, 48), 
  plot_width = 6, 
  plot_height = 9, 
  output_dir = "~/Downloads"
)

