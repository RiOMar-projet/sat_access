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
  install.packages(c("ncdf4", "curl", "reshape2", "ggplot2"), repos = "https://cloud.r-project.org/")
}

# Activate libraries
library(ncdf4)    # For reading NetCDF files
suppressPackageStartupMessages(library(curl)) # For FTP download
library(reshape2) # For data reshaping
library(ggplot2)  # For visualization


# The download function ---------------------------------------------------

download_nc <- function(dl_var, dl_dates, output_dir, overwrite) {
  
  # Check date range
  if(length(dl_dates) > 2){
    stop("Please provide only one or two dates for the date range.")
  }
  
  # Set start and end dates for download
  start_date <- as.Date(dl_dates[1])
  if(length(dl_dates) == 2){
    end_date <- as.Date(dl_dates[2])
  } else{
    end_date <- start_date
  }
  
  # Iterate over each date in the range
  current_date <- start_date
  while(current_date <= end_date){
    
    # Prep date strings
    dl_date <- current_date
    dl_date_flat <- gsub("-", "", dl_date)
    url_year_doy <- paste0(substr(dl_date, start = 1, stop = 4),"/",strftime(as.Date(dl_date), format = "%j"))
    
    # Get general URL based on desired variables
    if(toupper(dl_var) %in% c("SPM", "SPIM", "CHL", "CHLA")){
      url_base <- "ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic"
    } else {
      stop("Variable not yet available")
    }
    
    # Get product specifics
    if(toupper(dl_var) %in% c("SPM", "SPIM")){
      file_name <- paste0(dl_date_flat,"-EUR-L4-SPIM-ATL-v01-fv01-OI.nc.bz2")
      url_product <- "EUR-L4-SPIM-ATL-v01"
      nc_file <- file.path(output_dir, paste0(dl_date_flat,"-EUR-L4-SPIM-ATL-v01-fv01-OI.nc"))
      nc_var_name <- "analysed_spim"
      var_label <- "SPM [g m-3]"
    } else if(toupper(dl_var) %in% c("CHL", "CHLA")){
      file_name <- paste0(dl_date_flat,"-EUR-L4-CHL-ATL-v01-fv01-OI.nc.bz2")
      url_product <- "EUR-L4-CHL-ATL-v01"
      nc_file <- file.path(output_dir, paste0(dl_date_flat,"-EUR-L4-CHL-ATL-v01-fv01-OI.nc"))
      nc_var_name <- "analysed_chl_a"
      var_label <- "chl a [mg m-3]"
    } else {
      stop("Variable not yet available")
    }
    
    # Assemble final URL
    url_final <- paste(url_base, url_product, url_year_doy, file_name, sep = "/")
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
      
      # Unzip
      if(file.exists(file_name_full)){
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
    
    # Move to the next day
    current_date <- current_date + 1
  }
}


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

