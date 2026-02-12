# sat_access_script.R

# Title: Satellite Data Access and Visualization Script
# Author: Robert Schlegel
# Date: January 2026
# Description: This script contains functions that can be used to
#              download SPM or chlorophyll-a NetCDF files from a specified URL
#              for given date ranges and plot the data.
#              The products available via this script are SEXTANT and all files
#              available on the Acri-ST ftp server for ODATIS products.

# Useful documentation:
# ODATIS-MR user guide (en français) :
# https://www.aviso.altimetry.fr/fileadmin/documents/data/tools/hdbk_ODATIS_MR.pdf
# ODATIS-MR URL interface :
# https://tds-odatis.aviso.altimetry.fr/thredds/catalog/dataset-l3-ocean-color-odatis-mr-v_1_0.xml/catalog.html

# rschlegel
# MuDMxqHn6p9p

# Libraries ---------------------------------------------------------------

# Check for missing libraries and install them if necessary
if (!all(c("ncdf4", "curl", "lubridate", "reshape2", "maps", "ggplot2") %in% installed.packages())) {
  install.packages(c("ncdfCF", "curl", "lubridate", "reshape2", "maps", "ggplot2"), repos = "https://cloud.r-project.org/")
}

# Activate libraries
library(ncdf4)     # For accessing/reading NetCDF files
suppressPackageStartupMessages(library(curl)) # For FTP downloading
library(lubridate) # For working with dates
library(reshape2)  # For data reshaping
library(ggplot2)   # For visualization


# The download function ---------------------------------------------------

download_nc <- function(dl_var, dl_dates, 
                        dl_product = NULL, dl_sensor = NULL, 
                        dl_correction = NULL, dl_time_step = NULL, dl_bbox = NULL,
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
    stop("'dl_time_step' value not recognised, please choose one of : 'day', '8-day', 'month'")
  }
  
  # Set start and end dates for download
  # NB: 8-day files are initially handled as daily timesteps
  if(dl_time_step %in% c("day", "8-day")){
    start_date <- as.Date(dl_dates[1])
    if(length(dl_dates) == 2){
      end_date <- as.Date(dl_dates[2])
    } else{
      end_date <- start_date
    }
  } else if(dl_time_step == "month"){
    # NB: Always start with the first day of the year via floor_date()
    # Later in the script the full range of days for the chosen month are accounted for
    start_date <- floor_date(ymd(dl_dates[1]), "month")
    if(length(dl_dates) == 2){
      end_date <- floor_date(ymd(dl_dates[2]), "month")
    } else{
      end_date <- start_date
    }
  } else {
    stop("'dl_time_step value not recognised, please choose one of : 'day', '8-day', 'month'")
  }
  
  # Check if dl_bbox was given and return message accordingly
  if(is.null(dl_bbox)){
    # Intentionally blank
  } else if(length(dl_bbox) != 4){
    stop("Please provide a bounding box in the form of a vector with 4 values: c(lonmin, lonmax, latmin, latmax)")
  } else if(any(dl_bbox < -180) | any(dl_bbox > 180)){
    stop("Please ensure that the longitude values in 'dl_bbox' are between -180 and 180.")
  } else if(any(dl_bbox < -90) | any(dl_bbox > 90)){
    stop("Please ensure that the latitude values in 'dl_bbox' are between -90 and 90.")
  } else if(dl_bbox[1] >= dl_bbox[2]){ 
    stop("Please ensure that the first value in 'dl_bbox' (lonmin) is less than the second value (lonmax).") 
  } else if(dl_bbox[3] >= dl_bbox[4]){ 
    stop("Please ensure that the third value in 'dl_bbox' (latmin) is less than the fourth value (latmax).") 
  }
  
  # Check that start_date and end_date objects are proper date format
  if(!inherits(start_date, "Date") | !inherits(end_date, "Date")){
    stop("Please ensure that 'dl_dates' are provided in 'YYYY-MM-DD' format.")
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
  # TODO: Add a logic gate that checks the extent given versus what is available in ODATIS-MR
  # c(-7.8, 10.3, 41.2, 51.5) # The ODATIS-MR extent
  if(dl_product == "ODATIS-MR"){
    
    if(is.null(dl_bbox)){
      message("No bounding box provided, data will be downloaded for the full extent of the product.")
    } else {
      message(paste("Bounding box provided, data will be downloaded from longitude ",
                    dl_bbox[1]," to ",dl_bbox[2]," and latitude ",dl_bbox[3]," to ",dl_bbox[4],".", sep = ""))
    }
    
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
        message(paste0("No atmospheric correction chosen, defaulting to 'polymer' to match 'dl_sensor = ",dl_sensor,"'."))
        dl_correction <- "polymer"
      }
    } else if(!dl_correction %in% c("polymer", "nirswir")){
      stop("Please set 'dl_correction' to either 'polymer' or 'nirswir'")
    }
  }
  
  # Check that the chosen product, sensor, correction combos are possible
  if(dl_product == "SEXTANT"){
    if(!is.null(dl_bbox)){
      message("A bounding box was provided, but SEXTANT products cannot be subset at the source. Downloading the full file.")
    }
    if(!toupper(dl_var) %in% c("SPM", "SPIM", "CHL", "CHLA")){
      stop("SEXTANT data product only contains SPM and Chl a data. Please adjust `dl_var` accordingly.")
    }
  } else if(dl_product == "ODATIS-MR"){
    if(dl_sensor == "MODIS"){
      if(!toupper(dl_var) %in% c("CDOM", "CHL", "CHLA", "RRS", "NRRS", "SPM", "SPIM", "SST", "T", "TUR")){
        stop("ODATIS-MR MODIS data product does not contain the requested variable. Please adjust your variable choice accordingly.")
      }
    } else if(dl_sensor %in% c("MERIS", "OLCI-A", "OLCI-B")){
      if(!toupper(dl_var) %in% c("CDOM", "CHL", "CHLA", "RRS", "NRRS", "SPM", "SPIM", "T", "TUR")){
        stop(paste0("ODATIS-MR : ",dl_sensor," data product does not contain the requested variable. Please adjust your variable choice accordingly."))
      }
    }
  }
  
  # Check that the chosen correction against chosen data product and sensor and correct if necessary
  if(dl_product == "ODATIS-MR"){
    if(dl_sensor == "MODIS" & dl_correction != "nirswir"){
      message("ODATIS-MR : MODIS data product only uses the 'nirswir' atmospheric correction. 'dl_correction' adjusted accordingly.")
      dl_correction <- "nirswir"
    } else if(dl_sensor %in% c("MERIS", "OLCI-A", "OLCI-B") & dl_correction != "polymer"){
      message(paste0("ODATIS-MR : ",dl_sensor," data product only uses the 'polymer' atmospheric correction. 'dl_correction' adjusted accordingly."))
      dl_correction <- "polymer"
    }
  }
  
  # Create floor and ceiling date objects based on which sensor has been selected
  # TODO: Increase data ranges when ODATIS-MR data are uploaded through to 2026
  if(dl_product == "SEXTANT"){
    floor_date_sensor <- as.Date("1998-01-01")
    ceiling_date_sensor <- as.Date(Sys.Date())-7
  } else if(dl_sensor == "MODIS"){
    floor_date_sensor <- as.Date("2002-07-04")
    ceiling_date_sensor <- as.Date("2023-12-31")
  } else if(dl_sensor == "MERIS"){
    floor_date_sensor <- as.Date("2002-06-19")
    ceiling_date_sensor <- as.Date("2012-04-08")
  } else if(dl_sensor == "OLCI-A"){
    floor_date_sensor <- as.Date("2016-04-26")
    ceiling_date_sensor <- as.Date("2023-12-31")
  } else if(dl_sensor == "OLCI-B"){
    floor_date_sensor <- as.Date("2018-05-15")
    ceiling_date_sensor <- as.Date("2023-12-31")
  }
  
  # Check chosen dates against chosen data product, sensor, and correction
  if(start_date < floor_date_sensor | end_date > ceiling_date_sensor){
    stop(paste0("The chosen date range is outside of the available data range for : ",
                  dl_product, " - ", dl_sensor, " - ", dl_correction, ". \n",
                "Available data range is from ",floor_date_sensor," to ",ceiling_date_sensor,". Please adjust your date range accordingly."))
  }
  
  # Check if SEXTANT data are being requested and correct time step to daily if required
  if(dl_product == "SEXTANT" & dl_time_step != "day"){
    message("SEXTANT data product only contains daily data. 'dl_time_step' adjusted accordingly.")
    dl_time_step <- "day"
  }
  
  # TODO: 
  # Create a more intelligent logic gate that chooses the best product, sensor, and correction based on the variable and download dates chosen
  # It must create messages for the user letting them know what has been chosen and why
  # This would address the TODO item of choosing the 'best' chl-a etc. data given the bbox etc.
  
  # Iterate over each date in the range
  current_date <- start_date
  while(current_date <= end_date){
    
    # Prep date strings
    # NB: Will intentionally be overwritten below if requesting 8-day data
    dl_date <- current_date
    dl_date_flat <- gsub("-", "", dl_date)
    dl_year <- substr(dl_date, start = 1, stop = 4)
    dl_month <- substr(dl_date, start = 6, stop = 7)
    dl_day <- substr(dl_date, start = 9, stop = 10)
    url_year_doy <- paste0(substr(dl_date, start = 1, stop = 4),"/",strftime(as.Date(dl_date), format = "%j"))
    
    # Get URL specifics for SEXTANT
    if(toupper(dl_product) == "SEXTANT"){
      
      # Base URL
      # TODO: See if the data are stored anywhere with OpeNDAP access
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
      if(is.null(dl_bbox)){
        url_base <- "https://tds%40odatis-ocean.fr:odatis@tds-odatis.aviso.altimetry.fr/thredds/fileServer/dataset-l3-ocean-color-odatis-mr-v_1_0.xml/FRANCE"
      } else {
        url_base <- "https://tds%40odatis-ocean.fr:odatis@tds-odatis.aviso.altimetry.fr/thredds/dodsC/dataset-l3-ocean-color-odatis-mr-v_1_0.xml/FRANCE"
      }
      
      # Prep sensor strings
      dl_sensor_flat <- tolower(gsub("-", "", dl_sensor))
      if(grepl("OLCI", dl_sensor)){
        dl_sensor_chunk <- gsub("CI-", "", dl_sensor)
      } else {
        dl_sensor_chunk <- substr(dl_sensor, start = 1, stop = 3)
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
      
      # Prep the variable string
      if(dl_var %in% c("CDOM")){
        dl_var_chunk <- dl_var
      } else if(dl_var %in% c("CHL", "CHLA")){
        dl_var_chunk <- "CHL-OC5"
      } else if(dl_var %in% c("RRS", "NRRS")){
        if(dl_sensor == "MODIS"){
          dl_var_chunk <- "NRRS555"
        } else {
          dl_var_chunk <- "NRRS560"
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
      
      # Prep time-step string
      if(dl_time_step == "day"){
        dl_time_step_chunk <- "DAY"
      } else if(dl_time_step == "8-day"){
        dl_time_step_chunk <- "8D"
        
        # Get directory content
        url_8day_dir <- paste("https://tds-odatis.aviso.altimetry.fr/thredds/catalog/dataset-l3-ocean-color-odatis-mr-v_1_0.xml/FRANCE", 
                              dl_correction_flat, dl_sensor_flat, dl_time_step, dl_year, dl_month, "catalog.html", sep = "/")
        url_8day_catalog <- grep("\\d{2}/catalog\\.html", readLines(url_8day_dir), value = TRUE)
        url_8day_days <- unique(unlist(regmatches(url_8day_catalog, gregexpr("\\d+", url_8day_catalog))))
        url_8day_dates <- as.Date(paste(dl_year, dl_month, url_8day_days, sep = "/"))
        if(length(url_8day_dates) < 3 | length(url_8day_dates) > 5) stop("Something has gone wrong with the 8-day HTML scraping.")
        
        # Select closest date
        date_diffs <- abs(as.numeric(url_8day_dates - current_date))
        dl_date <- url_8day_dates[which.min(date_diffs)]
        dl_day <- substr(dl_date, start = 9, stop = 10)
        dl_date_flat <- paste0(format(dl_date, "%Y%m%d"), "-", format(dl_date+7, "%Y%m%d"))
        
        # Overwrite current date to match the 8-day value
        current_date <- dl_date
        
        # Double check that the product ceiling date has not been surpassed and correct accordingly
        # Note that the last file for every year in the ODATIS-MR product catalog is not 8 days of data
        # It is whatever the start date is, up to the last date of the year, or the last date of the product itself
        if(year(dl_date+7) > year(dl_date)){
          dl_date_flat <- paste0(format(dl_date, "%Y%m%d"), "-", format(as.Date(paste(dl_year, 12, 31, sep = "/")), "%Y%m%d"))
        }
        if(dl_date+7 > ceiling_date_sensor){
          dl_date_flat <- paste0(format(dl_date, "%Y%m%d"), "-", format(ceiling_date_sensor, "%Y%m%d"))
        }
        
      } else if(dl_time_step == "month"){
        dl_time_step_chunk <- "MO"
        dl_date_flat <- paste0(format(floor_date(ymd(dl_date), "month"), "%Y%m%d"), "-", format(ceiling_date(ymd(dl_date), "month") - days(1), "%Y%m%d"))
      } else {
        stop("Please check the value given for 'dl_time_step'")
      }
      
      # Product URL
      url_product <- paste(dl_correction_flat, dl_sensor_flat, dl_time_step, dl_year, dl_month, dl_day, sep = "/")
      file_name <- paste0("L3m_",dl_date_flat,"__FRANCE_03_",dl_sensor_chunk,"_",dl_var_chunk,"-",dl_correction_chunk,"_",dl_time_step_chunk,"_00.nc")
      # nc_file <- file.path(output_dir, file_name)
      
    } else {
      stop("Please check the value used for 'dl_product'")
    }

    # Assemble final URL
    url_final <- paste(url_base, url_product, file_name, sep = "/")
    file_name_full <- file.path(output_dir, file_name)
    
    # Fetch file
    if(file.exists(file_name_full) & !overwrite){
      
      message(paste0(file_name_full," already exists. Set 'overwrite = TRUE' to force the download."))
      
    } else {
      
      # Ensure output directory exists
      if(!dir.exists(output_dir)){
        dir.create(output_dir, recursive = TRUE)
      }
      
      # Check if file already exists and delete it if necessary
      if(file.exists(file_name_full)){
        tryCatch({ file.remove(file_name_full) 
          message(paste("File removed :", file_name_full))
        }, error = function(e) {
            message(paste("Impossible to remove file :", e$message)) 
        })
      }
      
      # Subset data by bounding box if possible/requested
      if(!is.null(dl_bbox) & dl_product == "ODATIS-MR"){
        
        # Open the connection
        nc_data <- nc_open(url_final)
        # names(nc_data)
        
        # Extract the longitude and latitude variables
        lon <- ncvar_get(nc_data, "lon")
        lat <- ncvar_get(nc_data, "lat")
        
        # Find the indices for the bounding box
        lon_indices <- which(lon >= dl_bbox[1] & lon <= dl_bbox[2])
        lat_indices <- which(lat >= dl_bbox[3] & lat <= dl_bbox[4])
        
        # Subset lon/lat
        lon_subset <- lon[lon_indices] 
        lat_subset <- lat[lat_indices]
        
        # Determine the start and count for reading the subset
        lonlat_start <- c(min(lon_indices), min(lat_indices))  # Assuming time is the third dimension
        lonlat_count <- c(length(lon_indices), length(lat_indices))  # -1 means read all along this dimension
        
        # Define the lon/lat dimensions
        lon_dim <- ncdim_def(name = "lon", units = "degrees_east", vals = lon_subset)
        lat_dim <- ncdim_def(name = "lat", units = "degrees_north", vals = lat_subset)
        
        # Get global attributes from the original file
        global_attrs <- ncatt_get(nc_data, varid = 0)
        
        # Initialize empty lists
        var_data_list <- list()
        var_def_list <- list()
        
        # Loop to assign values and definitions
        for (i in 1:length(names(nc_data$var))) {
          var_info_i <- nc_data$var[i][[1]]
          var_data_list[[paste0("var_data_", i)]] <- ncvar_get(nc = nc_data, varid = var_info_i$name, 
                                                               start = lonlat_start, count = lonlat_count)
          var_def_list[[paste0("var_def_", i)]] <- ncvar_def(name = var_info_i$name, units = var_info_i$units,
                                                             dim = list(lon = lon_dim, lat = lat_dim),
                                                             prec = var_info_i$prec, missval = var_info_i$missval, 
                                                             longname = var_info_i$longname)
        }; rm(i, var_info_i)

        # Close the NetCDF file
        nc_close(nc_data)
        
        # Create a new NetCDF file for the subset
        output_nc_data <- nc_create(file_name_full, vars = var_def_list, force_v4 = TRUE)
        
        # Add the data
        for (i in 1:length(var_def_list)) {
          ncvar_put(nc = output_nc_data, varid = var_def_list[[i]], vals = var_data_list[[i]])
        }; rm(i)
        
        # Additional attributes
        # ncatt_put(output_nc_data, "lon", "axis", "X")
        # ncatt_put(output_nc_data, "lat", "axis", "Y")
        
        # Define dimensions in the new NetCDF file
        # ncvar_def(output_nc_data, "lon", "double", lon_indices, missval = NA, longname = "Longitude", units = "degrees_east")
        # ncvar_def(output_nc_data, "lat", "double", lat_indices, missval = NA, longname = "Latitude", units = "degrees_north")
        # ncvar_def(output_nc_data, "time", "double", 0, missval = NA, longname = "Time", units = "days since 1900-01-01")
        
        # Define the variable in the new NetCDF file
        # ncvar_def(output_nc_data, var_name, "double", list(lon = lon_indices, lat = lat_indices), missval = NA)
        
        # Copy global attributes to the new file
        for (attr in names(global_attrs)) {
          ncatt_put(output_nc_data, 0, attr, global_attrs[[attr]])
        }
        
        # Update lon/lat global attributes
        # test1 <- ncdf4::ncatt_get(output_nc_data, varid = 0)
        ncatt_put(output_nc_data, varid = 0, attname = "northernmost_latitude", attval = max(lat_subset), definemode = TRUE)
        ncatt_put(output_nc_data, varid = 0, attname = "southernmost_latitude", attval = min(lat_subset), definemode = TRUE)
        ncatt_put(output_nc_data, varid = 0, attname = "easternmost_longitude", attval = max(lon_subset), definemode = TRUE)
        ncatt_put(output_nc_data, varid = 0, attname = "westernmost_longitude", attval = min(lon_subset), definemode = TRUE)
        ncatt_put(output_nc_data, varid = 0, attname = "spatial_subset", definemode = TRUE,
                  attval = paste("Spatially subset on", Sys.Date(), "by", 
                                 "https://github.com/RiOMar-projet/sat_access/blob/main/sat_access_script.R"))
        
        # Write the subset data to the new NetCDF file
        # ncvar_put(output_nc_data, "longitude", lon[lon_indices])
        # ncvar_put(output_nc_data, "latitude", lat[lat_indices])
        # ncvar_put(output_nc_data, var_name, var_subset)
        
        # Close the new NetCDF file
        nc_close(output_nc_data)
        
        # Print a message indicating the process is complete
        message("Subsetting complete. The results are saved as : ", file_name_full)
        
      } else {
        # Download
        tryCatch({
          curl::curl_download(url_final, destfile = file_name_full)
          message(paste0("File downloaded at: ",file_name_full))
        }, error = function(e) {
          if(grepl("Server denied you to change to the given directory", e$message[1])){
            e$message <- " not found on server."
          }
          message(paste0(file_name," : ",e$message[1]))
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
    }
    
    # Move to next time step depending on day, 8-day, or month
    if(dl_time_step == "day"){
      current_date <- current_date + 1
    } else if(dl_time_step == "8-day"){
      current_date <- current_date + 8
    } else if(dl_time_step == "month"){
      current_date <- current_date %m+% months(1)
    }
  }
}


# The plotting function ---------------------------------------------------

plot_nc <- function(nc_file, bbox = NULL, 
                    plot_width = NULL, plot_height = NULL, 
                    output_dir) {
  
  # Check for file
  if(!file.exists(nc_file)){
    return(message(paste0(nc_file," cannot be found.")))
  }
  
  # Peek inside the file
  # ncdump::NetCDF(nc_file)
  # nc_var_info <- ncdump::NetCDF(nc_file)$variable
  # nc_var_atts <- ncdump::NetCDF(nc_file)$attribute$global
  
  # Get the relevant snippet from the file
  file_snippets <- unlist(strsplit(basename(nc_file), "_"))
  # SEXTANT files
  if(length(file_snippets) == 1){
    if(grepl("SPIM-ATL", nc_file)){
      nc_var_name <- "analysed_spim"
    } else if(grepl("CHL-ATL", nc_file)){
      nc_var_name <- "analysed_chl_a"
    } else {
      stop("File structure cannot be inferred from file name. Please ensure the file was downloaded via the function found in this script.")
    }
  # ODATIS-MR files
  } else if(length(file_snippets) == 9){
    nc_var_name <- paste0(file_snippets[7],"_mean")
  } else {
    stop("File structure cannot be inferred from file name. Please ensure the file was downloaded via the function found in this script.")
  }
  
  # Get plotting label based on variable name
  if(grepl("spim|SPIM|SPM", nc_var_name)){
    var_label <- "SPM [g m-3]"
  } else if(grepl("chl|CHL", nc_var_name)){
    var_label <- "Chl a [mg m-3]"
  } else if(grepl("CDOM", nc_var_name)){
    var_label <- "CDOM [m-1]"
  } else if(grepl("SST", nc_var_name)){
    var_label <- "SST [°C]"
  } else if(grepl("T-FNU", nc_var_name)){
    var_label <- "Turbidity [FNU]"
  } else if(grepl("NRRS560", nc_var_name)){
    var_label <- "Rrs at 560 nm [sr-1]"
  } else if(grepl("NRRS555", nc_var_name)){
    var_label <- "Rrs at 555 nm [sr-1]"
  } else {
    stop("Variable label cannot be inferred from variable name.")
  }
  
  # Open the NetCDF file
  nc_data <- nc_open(nc_file)
  
  # Extract longitude, latitude, and the specified variable
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat")
  var <- ncvar_get(nc_data, nc_var_name)
  
  # Get time attribute based on file structure
  if(grepl("SPIM-ATL|CHL-ATL", nc_file)){
    time <- ncvar_get(nc_data, "time")
    time_step <- "daily"
    plot_date <- as.Date(as.POSIXct(time, origin = "1998-01-01"))
  } else if(grepl("L3m_", nc_file)){
    start_date <- as.Date(ncatt_get(nc_data, varid = 0)[["period_start_day"]], format = "%Y%m%d")
    end_date <- as.Date(ncatt_get(nc_data, varid = 0)[["period_end_day"]], format = "%Y%m%d")
    time_step <- ncatt_get(nc_data, varid = 0)[["product_type"]]
    plot_date <- as.Date(ncatt_get(nc_data, varid = 0)[["period_start_day"]], format = "%Y%m%d")
  } else {
    stop("Date value cannot be inferred from file structure.")
  }
  
  # Create plot title based on time step
  if(time_step == "daily"){
    plot_title <- paste("Map of", nc_var_name, "on", plot_date)
    title_size <- 10
  } else {
    plot_title <- paste("Map of", nc_var_name, "from", start_date, "to", end_date)
    title_size <- 8
  }
  
  # Close the NetCDF file
  nc_close(nc_data)
  
  # Reshape data for ggplot
  var_df <- melt(var)
  names(var_df) <- c("lon_idx", "lat_idx", "value")
  var_df$lon <- lon[var_df$lon_idx]
  var_df$lat <- lat[var_df$lat_idx]
  
  # Set bounding box if none exists
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
  var_df_sub <- var_df_sub[!is.na(var_df_sub$value),]
  
  # Set plot name
  plot_name <- file.path(output_dir, paste0(nc_var_name,"_",plot_date,".png"))
  
  # Plot using ggplot2
  p <- ggplot(var_df_sub, aes(x = lon, y = lat, fill = value)) +
    annotation_borders(regions = c("France", "Spain", "Germany", "UK", "Italy", "Andorra", 
                                   "Portugal", "Belgium", "Netherlands", "Switzerland", 
                                   "Austria", "Ireland", "Greece", "Luxembourg", "Liechtenstein"), 
                       fill = "grey80") +
    geom_tile() +
    scale_fill_viridis_c(var_label) +
    labs(title = plot_title,
         subtitle = paste0("Source : ",nc_file),
         x = "Longitude [°E]", y = "Latitude [°N]") +
    coord_quickmap(xlim = bbox[1:2], ylim = bbox[3:4]) +
    theme_minimal() +
    theme(legend.position = "bottom", 
          plot.title = element_text(size = title_size),
          plot.subtitle = element_text(size = 6),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6))
  
  # Determine plot dimensions, save, and exit
  if(is.null(plot_width)) plot_width <- 6
  if(is.null(plot_height)) plot_height <- 5
  ggsave(filename = plot_name, plot = p, width = plot_width, height = plot_height)
  message(paste0("Image saved at: ",plot_name))
}


# Examples ----------------------------------------------------------------

## Downloading ------------------------------------------------------------

# download_nc() will attempt to pick the correct sensors etc based on the requested variable and dates

# NB: It is only necessary to give the function a value for 'dl_var', it will figure out the rest
# But if you want a specific product etc., these are the options
## Possible values for 'dl_var' : "SPM", "CHL", "CDOM", "RRS", "T", "SST" 
## Possible values for 'dl_product' : "SEXTANT", "ODATIS-MR"
## Possible values for 'dl_sensor' : "MODIS", "MERIS", "OLCI-A", "OLCI-B"
## Possible values for 'dl_correction' : "polymer", "nirswir" 
## Possible values for 'dl_time_step' : "day", "8-day", "month"

# Download a few days of SPM data
# NB: If no product is chosen, CHL and SPM default to SEXTANT
download_nc(
  dl_var = "SPM",
  dl_dates = c("2025-09-01", "2025-09-05"),
  output_dir = "~/Downloads/SEXTANT", # Change as desired/required
  overwrite = FALSE # Change to TRUE to force downloads
)

# Download a few days of CHL data
# NB: If no product is chosen, CHL and SPM default to SEXTANT
download_nc(
  dl_var = "CHL",
  dl_dates = c("2025-09-01", "2025-09-05"),
  output_dir = "~/Downloads", # Change as desired/required
  overwrite = FALSE # Change to TRUE to force downloads
)

# Download one day of SPM data from MODIS for a given bounding box
download_nc(
  dl_var = "SPM",
  dl_dates = c("2008-12-25"),
  dl_product = "ODATIS-MR",
  dl_sensor = "MODIS",
  dl_bbox = c(3, 4, 42.5, 44),
  output_dir = "~/Downloads/MODIS", # Change as desired/required
  overwrite = TRUE # Change to TRUE to force downloads
)

# Download one day of Chl a data from MERIS
download_nc(
  dl_var = "CHL",
  dl_dates = c("2008-12-25"),
  dl_product = "ODATIS-MR",
  dl_sensor = "MERIS",
  output_dir = "~/Downloads", # Change as desired/required
  overwrite = FALSE # Change to TRUE to force downloads
)

# Download a few days of CDOM data
download_nc(
  dl_var = "CDOM",
  dl_dates = c("2022-09-01", "2022-09-05"),
  output_dir = "~/Downloads",
  overwrite = FALSE
)

# Download one day of SST data
download_nc(
  dl_var = "SST",
  dl_dates = "2014-01-01",
  output_dir = "~/Downloads",
  overwrite = FALSE
)

# Download one day of turbidity data from OLCI-B
download_nc(
  dl_var = "T",
  dl_dates = "2019-01-01",
  dl_product = "ODATIS-MR",
  dl_sensor = "OLCI-B",
  output_dir = "~/Downloads",
  overwrite = FALSE
)

# Download RRS data from MODIS (W_nm 555)
download_nc(
  dl_var = "RRS",
  dl_dates = "2019-01-01",
  dl_product = "ODATIS-MR",
  dl_sensor = "MODIS",
  output_dir = "~/Downloads",
  overwrite = FALSE
)

# Download RRS data from OLCI-A (W_nm 560)
download_nc(
  dl_var = "RRS",
  dl_dates = "2019-01-01",
  dl_product = "ODATIS-MR",
  dl_sensor = "OLCI-A",
  output_dir = "~/Downloads",
  overwrite = FALSE
)

# Download a few 8-day averages of Chl a data
download_nc(
  dl_var = "CHL",
  dl_dates = c("2018-12-01", "2019-01-25"),
  dl_product = "ODATIS-MR",
  dl_sensor = "OLCI-A",
  dl_time_step = "8-day",
  output_dir = "~/Downloads",
  overwrite = FALSE
)

# Download a few months of Turbidity data from OLCI-A
download_nc(
  dl_var = "T",
  dl_dates = c("2019-01-01", "2019-03-01"),
  dl_product = "ODATIS-MR",
  dl_sensor = "OLCI-A",
  dl_time_step = "month",
  output_dir = "~/Downloads",
  overwrite = FALSE
)


## Plotting ---------------------------------------------------------------

# NB: plot_nc() will automagically detect the necessary .nc structure and variable from the file name
# Note that it is only designed to work with files downloaded via download_nc()

# Plot an SPM NetCDF file
plot_nc(
  nc_file = "~/Downloads/MODIS/L3m_20081225__FRANCE_03_MOD_SPM-G-NS_DAY_00.nc", 
  bbox = c(3, 4, 42.5, 44),
  plot_width = 5, 
  plot_height = 5, 
  output_dir = "~/Downloads"
)

# Plot a Chl a NetCDF file
plot_nc(
  nc_file = "~/Downloads/20251225-EUR-L4-CHL-ATL-v01-fv01-OI.nc", 
  bbox = c(-3.5, -0.5, 44, 48), 
  plot_width = 5, 
  plot_height = 9, 
  output_dir = "~/Downloads"
)

# Plot a CDOM NetCDF file
# NB: Providing no 'bbox' argument will plot the full extent of the data in the file
plot_nc(
  nc_file = "~/Downloads/L3m_20220901__FRANCE_03_OLA_CDOM-PO_DAY_00.nc", 
  plot_width = 7, 
  plot_height = 6, 
  output_dir = "~/Downloads"
)

# Plot an SST NetCDF file
plot_nc(
  nc_file = "~/Downloads/L3m_20140101__FRANCE_03_MOD_SST-NIGHT-NS_DAY_00.nc",
  bbox = c(-5, 10, 40, 55), 
  plot_width = 8, 
  plot_height = 6, 
  output_dir = "~/Downloads"
)

# Plot a turbidity NetCDF file
plot_nc(
  nc_file = "~/Downloads/L3m_20190101__FRANCE_03_OLB_T-FNU-PO_DAY_00.nc",
  bbox = c(2, 6, 41, 44),
  plot_width = 6,
  plot_height = 5,
  output_dir = "~/Downloads"
)

# Plot an RRS NetCDF file from MODIS
plot_nc(
  nc_file = "~/Downloads/L3m_20190101__FRANCE_03_MOD_NRRS555-NS_DAY_00.nc",
  bbox = c(2, 5, 41, 45),
  plot_width = 7,
  plot_height = 6,
  output_dir = "~/Downloads"
)

# Plot an RRS NetCDF file from OLCI-A
plot_nc(
  nc_file = "~/Downloads/L3m_20190101__FRANCE_03_OLA_NRRS560-PO_DAY_00.nc",
  plot_width = 7,
  plot_height = 6,
  output_dir = "~/Downloads"
)

# Plot an 8-day average of Chl a data
plot_nc(
  nc_file = "~/Downloads/L3m_20181219-20181226__FRANCE_03_OLA_CHL-OC5-PO_8D_00.nc",
  # bbox = c(2, 6, 41, 44),
  plot_width = 5,
  plot_height = 5,
  output_dir = "~/Downloads"
)

# Plot one month of turbidity data
plot_nc(
  nc_file = "~/Downloads/L3m_20190201-20190228__FRANCE_03_OLA_T-FNU-PO_MO_00.nc",
  # bbox = c(2, 6, 41, 44),
  plot_width = 5,
  plot_height = 5,
  output_dir = "~/Downloads"
)

