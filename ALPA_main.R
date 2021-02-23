#
# ###################################################################### #
# ###################################################################### #
#                                                                        #
#                                ALPA                                    #
#                                                                        #
#      An R-Script for for automatic analysis of landslide profile       #
#                                                                        #
#                             Version 2.0                                #
#                                                                        #
#                          February 3rd, 2021                            #
#                                                                        #
#                       Langping LI, Hengxing LAN                        #
#                                                                        #
#                           LREIS, IGSNRR, CAS                           #
#                                                                        #
#               Email: lilp@lreis.ac.cn, lanhx@lreis.ac.cn               #
#                                                                        #
# ###################################################################### #
# ###################################################################### #
#
#
#
# ########## commonly it is NOT SUGGESTED to change this script ##########
#
#
#
ALPA <- function(FileDEM, FileLandslides, MinGrpPntCnt = 3, MinGrpAcrDst = 0, MinEndAptRto = 0, MaxEndDvtAgl = 180, FilePrefix = "", OutputTemp = F) {
  #
  # check input
  if (dendextend::is.natural.number(MinGrpPntCnt) == F) {
    #
    # warning
    warning("Error: MinGrpPntCnt is not a natural number.")
    #
    # return
    return(NA)
  }
  #
  # check input
  if (MinGrpAcrDst < 0) {
    #
    # warning
    warning("Error: MinGrpAcrDst is negative.")
    #
    # return
    return(NA)
  }
  #
  # import vector layer
  shp_lasld <- rgdal::readOGR(FileLandslides)
  #
  # get polygons
  list_polygons <- shp_lasld@polygons
  #
  # get count
  count_polygons <- length(list_polygons)
  #
  # initial list_spldf_profile
  list_spldf_profile <- list()
  #
  # get CRS
  raster_dem <- raster::raster(FileDEM)
  raster_CRS <- sp::proj4string(raster_dem)
  #
  # get file
  file_out_profile <- paste(FilePrefix, "_lasld_profile2d", sep = "", collapse = NULL)
  file_out_xlsx <- paste(FilePrefix, "_lasld_profile2d.xlsx", sep = "", collapse = NULL)
  #
  # split
  for (i in 1:count_polygons) {
    #
    # get time
    time_lasld_start <- Sys.time()
    # #
    # show progress
    # cat( paste("Process the ", as.character(i), "/", as.character(count_polygons), " landslide start.", "\n", sep = "") )
    #
    # get Polygons
    pPolygons <- list_polygons[[i]]
    #
    # split using on Polygons
    spldf_profile <- f_lasld_split(raster_dem, pPolygons, MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, FilePrefix, i, OutputTemp)
    #
    # update list
    if (is.null(spldf_profile) == F) {
      #
      # append list
      list_spldf_profile[[length(list_spldf_profile)+1]] <- spldf_profile
      #
      # list2shp
      data_df <- f_save_list2spldf(list_spldf_profile, raster_CRS, file_out_profile)
      #
      # save xlsx
      f_pnts_save_xls(data_df, file_out_xlsx)
      #
      # show progress
      cat( paste("Process the ", as.character(i), "/", as.character(count_polygons), " landslide done.", "\n", sep = "") )
    }
    else {
      #
      # show progress
      cat( paste("Process the ", as.character(i), "/", as.character(count_polygons), " landslide failed.", "\n", sep = "") )
    }
    #
    # get time
    time_lasld_end <- Sys.time()
    #
    # show time
    cat(paste(as.character(as.numeric(time_lasld_end - time_lasld_start,  units = "secs")), " seconds used.", "\n", sep = "", collapse = NULL))
  }
}
#
#
#
f_lasld_split <- function(pRasterDEM, pPolygons, MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, FilePrefix, FileID, OutputTemp) {
  # #
  # # show step
  # time_step_start <- f_show_step(paste("Read pnts start.", "\n", sep = ""))
  #
  # initial
  spldf_profile <- NULL
  #
  # get CRS
  raster_CRS <- sp::proj4string(pRasterDEM)
  #
  # get cell size
  raster_cellsize <- (raster::res(pRasterDEM))[1]
  #
  # get SpatialPolygons
  pSpatialPolygons <- sp::SpatialPolygons(list(pPolygons), proj4string = sp::CRS(raster_CRS))
  #
  # rasterize polygon
  list_polygon_rasterized <- f_polygon_rasterize(pRasterDEM, pSpatialPolygons)
  #
  # get
  shp_lasld <- list_polygon_rasterized[[1]]
  # plot(shp_lasld, main = "polygon", axes = TRUE, border = "blue")
  #
  # set possible file_out (for the line of rasterized landslide plolygon)
  file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_lasld_pln2d", sep = "", collapse = NULL)
  #
  # get
  pRasterDEM_masked <- list_polygon_rasterized[[2]]
  #
  # get
  spxdf_lasld <- as(pRasterDEM_masked, "SpatialPixelsDataFrame")
  #
  # get pnts
  pnts <- f_pnts_read(pRasterDEM_masked, shp_lasld, file_out)
  #
  # if read pnts fail
  if (max(pnts[, "IDB"]) == 0) {
    #
    # error
    str_failed <- paste("Error: Read pnts failed for the ", as.character(FileID), "th landslide.", sep = "", collapse = NULL)
    warning(str_failed)
    #
    # save pnts
    file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts", sep = "", collapse = NULL)
    f_pnts_save_points(pnts, raster_CRS, file_out)
    #
    # return
    return(NULL)
  }
  # #
  # # show step time
  # f_show_time(paste("Read pnts done.", "\n", sep = ""), time_step_start)
  #
  # save the original pnts
  # if (OutputTemp) {
  if (T) {
    #
    # save pnts
    file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts", sep = "", collapse = NULL)
    f_pnts_save_points(pnts, raster_CRS, file_out)
  }
  #
  # initial pnts_split
  pnts_split <- pnts
  #
  # if normal is vertical (horizontal plane)
  if (pnts_split[1, "Dx"] == 0 && pnts_split[1, "Dy"] == 0 && pnts_split[1, "Dz"] == 0) { 
    #
    # error, if when reading pnts
    warning("Error: The first fitted plane for the landslide is horizontal.")
    #
    # set SMS, to be "LSH" (horizontal landslide)
    pnts_split[, "SMS"] <- "LSH"
    #
    # save pnts
    file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts_split, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp", sep = "", collapse = NULL)
    f_pnts_save_points(pnts_split, raster_CRS, file_out)
    #
    # return
    return(NULL)
  }
  #
  # initial list_pnts
  list_pnts <- list(pnts_split)
  #
  # initial list_pntu, pntu is a pnt marking the position of the upper stream, for each pnts
  # list_pntu <- list(pnts_split[which(pnts_split[, "Cz"] == max(pnts_split[, "Cz"])), ])
  list_pntu <- list(pnts_split[which(pnts_split[, "zip"] == max(pnts_split[, "zip"])), ])
  #
  # initial list_grpm, grpm is a marker indicating subgrouping state, for each pnts (0 means need group)
  list_grpm <- list(c(0))
  #
  # get the maximum ID of the boundary points of the landslide
  IDB_max <- max(pnts_split[, "IDB"])
  #
  # subgrouping
  while (TRUE) {
    #
    # save for every split
    if (OutputTemp) {
      #
      # save split
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), ".xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split, raster_CRS, file_out)
    }
    #
    # get centers
    pnts_split_anchors <- f_pnts_anchors(list_pnts, IDB_max)
    #
    # save for every split
    if (OutputTemp) {
      #
      # save anchors
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_anchors.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split_anchors, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_anchors", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split_anchors, raster_CRS, file_out)
      #
      # save shp pln
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_anchors_pln2d", sep = "", collapse = NULL)
      spldf_profile <- f_pnts_save_proflle(pnts_split_anchors, spxdf_lasld, file_out, FileID, c(MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl))
      #
      # save paras
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_anchors_pln2d.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(spldf_profile@data, file_out)
    }
    #
    # break, if no group need split
    if (min(data.frame(list_grpm)) == 1) {
      # 
      # set SMS, to be "SMG", stop because all subgroup stop split
      pnts_split[, "SMS"] <- "SMG"
      #
      # save pnts
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split, raster_CRS, file_out)
      #
      # save anchors
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp_anchors.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split_anchors, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp_anchors", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split_anchors, raster_CRS, file_out)
      #
      # save shp pln
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp_anchors_pln2d", sep = "", collapse = NULL)
      spldf_profile <- f_pnts_save_proflle(pnts_split_anchors, spxdf_lasld, file_out, FileID, c(MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl))
      #
      # save paras
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grp_anchors_pln2d.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(spldf_profile@data, file_out)
      #
      # break
      break
    }
    #
    # initial list_pnts_new
    list_pnts_new <- list()
    #
    # initial list_pntu_new
    list_pntu_new <- list()
    #
    # initial list_grpm_new
    list_grpm_new <- list()
    #
    # get the maximum grp of pnts_split
    grp_max <- max(pnts_split[, "grp"])
    #
    # initial grp (group ID) of the upper stream subgroup (sub- pnts)
    grpu <- 1
    #
    # for initial pnts
    pnts_initial <- list_pnts[[1]]
    pnts_initial_bnd <- pnts_initial[which(pnts_initial[, "IDB"] != 0), ]
    #
    # for distal pnts
    pnts_distal <- list_pnts[[length(list_pnts)]]
    pnts_distal_bnd <- pnts_distal[which(pnts_distal[, "IDB"] != 0), ]
    #
    # subgroup all pnts
    for (i in 1:length(list_pnts)) {
      #
      # get pnts_i
      pnts_i <- list_pnts[[i]]
      #
      # get upper stream pnt
      pntu_i <- list_pntu[[i]]
      #
      # get group marker
      grpm_i <- list_grpm[[i]]
      #
      # if still need subgrouping
      if (grpm_i == 0) {
        #
        # get grpm
        grpm <- grpm_i
        # #
        # # show step
        # time_step_start <- f_show_step(paste("Split pnts start.", "\n", sep = ""))
        #
        # get intercept of the split plane
        pnts_i_intercept <- f_pnts_intercept(pnts_i)
        #
        # get list_pnts_i_split
        list_pnts_i_split <- f_pnts_split(pnts_i, pnts_i_intercept)
        #
        # get pnts
        pnts1 <- list_pnts_i_split[[1]]
        pnts2 <- list_pnts_i_split[[2]]
        #
        # get pnts
        pnts_i_upper <- pnts1
        pnts_i_lower <- pnts2
        #
        # upper or lower
        if ( pracma::isempty( which(pnts1[, "ID"] == pntu_i[1, "ID"]) ) ) {
          #
          # switch upper and lower
          pnts_i_upper <- pnts2
          pnts_i_lower <- pnts1
        }
        # #
        # # show step time
        # f_show_time(paste("Split pnts done.", "\n", sep = ""), time_step_start)
        # #
        # # show step
        # time_step_start <- f_show_step(paste("Check spatial contiguity start.", "\n", sep = ""))
        #
        # get clusters, check spatial contiguities
        list_clusters_upper <- f_pnts_clustering(pnts_i_upper, raster_cellsize, raster_CRS)
        list_clusters_lower <- f_pnts_clustering(pnts_i_lower, raster_cellsize, raster_CRS)
        #
        # save original
        pnts_i_upper_original <- pnts_i_upper
        pnts_i_lower_original <- pnts_i_lower
        #
        # allow for both upper and lower subset
        # possibly both list_clusters_upper and list_clusters_lower >= 2
        # can process reassign together
        if (length(list_clusters_upper) >= 2 || length(list_clusters_lower) >= 2) {
          # #
          # # debugging
          # f_pnts_save_points(pnts_i_upper, raster_CRS, "pnts_i_upper")
          # f_pnts_save_points(pnts_i_lower, raster_CRS, "pnts_i_lower")
          #
          # initial
          list_clusters_upper_for_lower <- list()
          #
          # update the upper subset
          for (k in 1:length(list_clusters_upper)) {
            #
            # get pnts
            pnts_cluster <- list_clusters_upper[[k]]
            #
            # assign
            if ( pracma::isempty( which(pnts_cluster[, "ID"] == pntu_i[1, "ID"]) ) ) {
              #
              # to lower
              # pnts_i_lower <- rbind(pnts_i_lower, pnts_cluster)
              list_clusters_upper_for_lower[[length(list_clusters_upper_for_lower) + 1]] <- pnts_cluster
            }
            else {
              #
              # to upper
              pnts_i_upper <- pnts_cluster
            }
          }
          #
          # initial
          # find the cluster in the lower subset with the largest count of pnts
          cluster_largest_ID <- 1
          cluster_largest_count <- nrow(list_clusters_lower[[cluster_largest_ID]])
          #
          # get the largest cluster, assign to the new lower subset
          for (k in 1:length(list_clusters_lower)) {
            #
            # get count
            cluster_count <- nrow(list_clusters_lower[[k]])
            #
            # update
            if (cluster_count > cluster_largest_count) {
              #
              # update
              cluster_largest_ID <- k
              cluster_largest_count <- cluster_count
            }
          }
          #
          # initial
          list_clusters_lower_for_upper <- list()
          #
          # update the lower subset
          for (k in 1:length(list_clusters_lower)) {
            #
            # get pnts
            pnts_cluster <- list_clusters_lower[[k]]
            #
            # assign
            if (k == cluster_largest_ID) {
              #
              # to lower
              pnts_i_lower <- pnts_cluster
            }
            else {
              #
              # to upper
              # pnts_i_upper <- rbind(pnts_i_upper, pnts_cluster)
              list_clusters_lower_for_upper[[length(list_clusters_lower_for_upper) + 1]] <- pnts_cluster
            }
          }
          #
          # update the upper subset again
          if (length(list_clusters_lower_for_upper) >= 1) {
            #
            # update
            for (k in 1:length(list_clusters_lower_for_upper)) {
              #
              # to upper
              pnts_i_upper <- rbind(pnts_i_upper, list_clusters_lower_for_upper[[k]])
            }
          }
          #
          # update the lower subset again
          if (length(list_clusters_upper_for_lower) >= 1) {
            #
            # update
            for (k in 1:length(list_clusters_upper_for_lower)) {
              #
              # to lower
              pnts_i_lower <- rbind(pnts_i_lower, list_clusters_upper_for_lower[[k]])
            }
          }
        }
        #
        # debugging
        if ( (nrow(pnts_i_upper) + nrow(pnts_i_lower) ) != nrow(pnts_i)) {
          #
          # warning
          warning("Error: count of pnts in clustering is not consistent.")
          #
          # return
          return(NULL)
        }
        #
        # get clusters, check spatial contiguities, again
        list_clusters_upper <- f_pnts_clustering(pnts_i_upper, raster_cellsize, raster_CRS)
        list_clusters_lower <- f_pnts_clustering(pnts_i_lower, raster_cellsize, raster_CRS)
        #
        # no spatial contiguities, even after combination
        # possibly both list_clusters_upper and list_clusters_lower >= 2, ???
        if (length(list_clusters_upper) >= 2 || length(list_clusters_lower) >= 2) {
          #
          # set grpm
          grpm <- 1
          #
          # set pnts_i
          pnts_i[, "SMG"] <- "SGT"
          pnts_i[, "grp"] <- grpu
          #
          # update list
          list_pnts_new[[grpu]] <- pnts_i
          list_pntu_new[[grpu]] <- pntu_i
          list_grpm_new[[grpu]] <- grpm
          #
          # update the upper group
          grpu <- grpu + 1
          #
          # continue
          next
          # #
          # # debugging
          # f_pnts_save_points(pnts_i_upper_original, raster_CRS, "pnts_i_upper_original")
          # f_pnts_save_points(pnts_i_lower_original, raster_CRS, "pnts_i_lower_original")
          # #
          # # debugging
          # f_pnts_save_points(pnts_i_upper, raster_CRS, "pnts_i_upper")
          # f_pnts_save_points(pnts_i_lower, raster_CRS, "pnts_i_lower")
          # #
          # # warning
          # warning("Error: subsets can be clustered after combination.")
          # #
          # # return
          # return(NULL)
        }
        # #
        # # show step time
        # f_show_time(paste("Check spatial contiguity done.", "\n", sep = ""), time_step_start)
        #
        # count of pnt in sub- pnts is smaller than a threshold
        if (nrow(pnts_i_upper) < MinGrpPntCnt || nrow(pnts_i_lower) < MinGrpPntCnt) {
          #
          # set grpm
          grpm <- 1
          #
          # set pnts_i
          pnts_i[, "SMG"] <- "SGP"
          pnts_i[, "grp"] <- grpu
          #
          # update list
          list_pnts_new[[grpu]] <- pnts_i
          list_pntu_new[[grpu]] <- pntu_i
          list_grpm_new[[grpu]] <- grpm
          #
          # update the upper group
          grpu <- grpu + 1
          #
          # continue
          next
        }
        #
        # get plane
        pnts_i_upper <- f_pnts_plane(pnts_i_upper)
        pnts_i_lower <- f_pnts_plane(pnts_i_lower)
        #
        # if can not get plane (collinear pnts)
        if (is.na(pnts_i_upper) || is.na(pnts_i_lower)) {
          #
          # set grpm
          grpm <- 1
          #
          # set pnts_i
          pnts_i[, "SMG"] <- "SGC"
          pnts_i[, "grp"] <- grpu
          #
          # update list
          list_pnts_new[[grpu]] <- pnts_i
          list_pntu_new[[grpu]] <- pntu_i
          list_grpm_new[[grpu]] <- grpm
          #
          # update the upper group
          grpu <- grpu + 1
          #
          # continue
          next
        }
        #
        # vertical normal (horizontal plane)
        if ((pnts_i_upper[1, "Dx"] == 0 && pnts_i_upper[1, "Dy"] == 0 && pnts_i_upper[1, "Dz"] == 0) ||
            (pnts_i_lower[1, "Dx"] == 0 && pnts_i_lower[1, "Dy"] == 0 && pnts_i_lower[1, "Dz"] == 0)) {
          #
          # set grpm
          grpm <- 1
          #
          # set pnts_i
          pnts_i[, "SMG"] <- "SGH"
          pnts_i[, "grp"] <- grpu
          #
          # update list
          list_pnts_new[[grpu]] <- pnts_i
          list_pntu_new[[grpu]] <- pntu_i
          list_grpm_new[[grpu]] <- grpm
          #
          # update the upper group
          grpu <- grpu + 1
          #
          # continue
          next
        }
        #
        # get boolEND
        boolEND_upper <- FALSE
        if (pnts_i_upper[1, "grp"] == 1) { boolEND_upper <- TRUE }
        boolEND_lower <- FALSE
        if (pnts_i_lower[1, "grp"] == grp_max) { boolEND_lower <- TRUE }
        # #
        # # show step
        # time_step_start <- f_show_step(paste("Check spatial impartiality start.", "\n", sep = ""))
        #
        # check if sub- pnts have boundary points on two sides
        grpm_i_upper <- f_pnts_bndcorners_check(pnts_i_upper, IDB_max, pnts_initial_bnd, pnts_distal_bnd, boolEND_upper)
        grpm_i_lower <- f_pnts_bndcorners_check(pnts_i_lower, IDB_max, pnts_initial_bnd, pnts_distal_bnd, boolEND_lower)
        #
        # if sub- pnts do not have boundary points on two sides
        if (grpm_i_upper == 1 || grpm_i_lower == 1) {
          # #
          # # debugging
          # f_pnts_save_points(pnts_i_upper, raster_CRS, "pnts_i_upper")
          # f_pnts_save_points(pnts_i_lower, raster_CRS, "pnts_i_lower")
          #
          # set grpm
          grpm <- 1
          #
          # set pnts_i
          pnts_i[, "SMG"] <- "SGB"
          pnts_i[, "grp"] <- grpu
          #
          # update list
          list_pnts_new[[grpu]] <- pnts_i
          list_pntu_new[[grpu]] <- pntu_i
          list_grpm_new[[grpu]] <- grpm
          #
          # update the upper group
          grpu <- grpu + 1
          #
          # continue
          next
        }
        #
        # initial or distal sub- pnts is not correct
        if (grpm_i_upper == 3 || grpm_i_lower == 3) {
          #
          # set grpm
          grpm <- 1
          #
          # set pnts_i
          pnts_i[, "SMG"] <- "SGE"
          pnts_i[, "grp"] <- grpu
          #
          # update list
          list_pnts_new[[grpu]] <- pnts_i
          list_pntu_new[[grpu]] <- pntu_i
          list_grpm_new[[grpu]] <- grpm
          #
          # update the upper group
          grpu <- grpu + 1
          #
          # continue
          next
        }
        #
        # if sub- pnts can be clustered, have more than two sub- culsters
        # theoretically, this situation will not appear???
        if (grpm_i_upper == 2 || grpm_i_lower == 2) {
          # #
          # # set grpm
          # grpm <- 1
          # #
          # # set pnts_i
          # pnts_i[, "SMG"] <- "SGT"
          # pnts_i[, "grp"] <- grpu
          # #
          # # update list
          # list_pnts_new[[grpu]] <- pnts_i
          # list_pntu_new[[grpu]] <- pntu_i
          # list_grpm_new[[grpu]] <- grpm
          # #
          # # update the upper group
          # grpu <- grpu + 1
          # #
          # # continue
          # next
          #
          # debugging
          f_pnts_save_points(pnts_i_upper_original, raster_CRS, "pnts_i_upper_original")
          f_pnts_save_points(pnts_i_lower_original, raster_CRS, "pnts_i_lower_original")
          #
          # debugging
          f_pnts_save_points(pnts_i_upper, raster_CRS, "pnts_i_upper")
          f_pnts_save_points(pnts_i_lower, raster_CRS, "pnts_i_lower")
          #
          # warning
          warning("Error: spatial in-contiguities are not expected here.")
          #
          # return
          return(NULL)
        }
        #
        # if check fail
        if (is.na(grpm_i_upper) || is.na(grpm_i_lower)) {
          #
          # set grpm
          grpm <- 1
          #
          # set pnts_i
          pnts_i[, "SMG"] <- "SGF"
          pnts_i[, "grp"] <- grpu
          #
          # update list
          list_pnts_new[[grpu]] <- pnts_i
          list_pntu_new[[grpu]] <- pntu_i
          list_grpm_new[[grpu]] <- grpm
          #
          # update the upper group
          grpu <- grpu + 1
          #
          # continue
          next
        }
        # #
        # # show step time
        # f_show_time(paste("Check spatial impartiality done.", "\n", sep = ""), time_step_start)
        #
        # make a copy
        list_pnts_new_anchor <- list_pnts_new
        #
        # update list_pnts_new
        list_pnts_new_anchor[[grpu]] <- pnts_i_upper
        list_pnts_new_anchor[[grpu + 1]] <- pnts_i_lower
        #
        # if still have pnts, append
        if ((i+1) <= length(list_pnts)) {
          #
          # append
          for (k in (i+1):length(list_pnts)) {
            #
            # get pnts_i
            pnts_k <- list_pnts[[k]]
            list_pnts_new_anchor[[length(list_pnts_new_anchor) + 1]] <- pnts_k
          }
        }
        #
        # get distances (before split)
        pnts_split_distances_before <- f_pnts_distances(pnts_split_anchors[, c("Cx", "Cy", "Cz")])
        #
        # get centers
        pnts_split_anchors <- f_pnts_anchors(list_pnts_new_anchor, IDB_max)
        #
        # get distances
        pnts_split_distances <- f_pnts_distances(pnts_split_anchors[, c("Cx", "Cy", "Cz")])
        #
        # initial
        pnts_split_distances_min <- Inf
        #
        # get distance minimum
        for (k in 1:nrow(pnts_split_distances)) {
          #
          # sometimes, centerP and centerL will be the same
          if (pnts_split_distances[k] != 0 && pnts_split_distances[k] < pnts_split_distances_min) {
            #
            # update
            pnts_split_distances_min <- pnts_split_distances[k]
          }
        }
        #
        # if minimum measurement scale
        if (pnts_split_distances_min <= MinGrpAcrDst) {
          #
          # set grpm
          grpm <- 1
          #
          # set pnts_i
          pnts_i[, "SMG"] <- "SGS"
          pnts_i[, "grp"] <- grpu
          #
          # update list
          list_pnts_new[[grpu]] <- pnts_i
          list_pntu_new[[grpu]] <- pntu_i
          list_grpm_new[[grpu]] <- grpm
          #
          # update the upper group
          grpu <- grpu + 1
          #
          # continue
          next
        }
        #
        # handle end initial
        if (boolEND_upper == TRUE) {
          #
          # get aspect ratio for end group before split
          EndAptRto_initial_before <- f_pnts_end_aspect_ratio(nrow(pnts_i) * raster_cellsize^2, pnts_split_distances_before[1] + pnts_split_distances_before[2])
          #
          # if adequate small
          if (EndAptRto_initial_before < MinEndAptRto) {
            #
            # get deviation angle for end group after split
            EndDvtAgl_initial <- f_pnts_anchors_deviation_angle(pnts_split_anchors, 1)
            #
            # if adequate large
            if (EndDvtAgl_initial > MaxEndDvtAgl) {
              #
              # set grpm
              grpm <- 1
              #
              # set pnts_i
              pnts_i[, "SMG"] <- "SGD"
              pnts_i[, "grp"] <- grpu
              #
              # update list
              list_pnts_new[[grpu]] <- pnts_i
              list_pntu_new[[grpu]] <- pntu_i
              list_grpm_new[[grpu]] <- grpm
              #
              # update the upper group
              grpu <- grpu + 1
              #
              # continue
              next
            }
          }
        }
        #
        # handle end distal
        if (boolEND_lower == TRUE) {
          #
          # get count
          cnt <- nrow(pnts_split_distances_before)
          #
          # get aspect ratio for end group before split
          EndAptRto_distal_before <- f_pnts_end_aspect_ratio(nrow(pnts_i) * raster_cellsize^2, pnts_split_distances_before[cnt] + pnts_split_distances_before[cnt-1])
          #
          # if adequate small
          if (EndAptRto_distal_before < MinEndAptRto) {
            #
            # get deviation angle for end group after split
            EndDvtAgl_distal <- f_pnts_anchors_deviation_angle(pnts_split_anchors, 2)
            #
            # if adequate large
            if (EndDvtAgl_distal > MaxEndDvtAgl) {
              #
              # set grpm
              grpm <- 1
              #
              # set pnts_i
              pnts_i[, "SMG"] <- "SGD"
              pnts_i[, "grp"] <- grpu
              #
              # update list
              list_pnts_new[[grpu]] <- pnts_i
              list_pntu_new[[grpu]] <- pntu_i
              list_grpm_new[[grpu]] <- grpm
              #
              # update the upper group
              grpu <- grpu + 1
              #
              # continue
              next
            }
          }
        }
        #
        # update grp
        pnts_i_upper[, "grp"] <- grpu
        pnts_i_lower[, "grp"] <- grpu + 1
        #
        # update list_pnts_new
        list_pnts_new[[grpu]] <- pnts_i_upper
        list_pnts_new[[grpu + 1]] <- pnts_i_lower
        #
        # get pntu_i_lower
        pnts_i_dip <- c(pnts_i[1, "Dx"], pnts_i[1, "Dy"], pnts_i[1, "Dz"])
        pntu_i_lower <- f_pnts_MinDist_intercept(pnts_i_lower, pnts_i_dip, pnts_i_intercept)
        #
        # update list_pntu_new
        list_pntu_new[[grpu]] <- pntu_i
        list_pntu_new[[grpu + 1]] <- pntu_i_lower
        #
        # update list_grpm_new
        list_grpm_new[[grpu]] <- grpm
        list_grpm_new[[grpu + 1]] <- grpm
        #
        # update the upper group
        grpu <- grpu + 2
      }
      else {
        #
        # set pnts_i
        pnts_i[, "grp"] <- grpu
        #
        # update list
        list_pnts_new[[grpu]] <- pnts_i
        list_pntu_new[[grpu]] <- pntu_i
        list_grpm_new[[grpu]] <- grpm_i
        #
        # update the upper group
        grpu <- grpu + 1
      }
    }
    #
    # update list_pnts
    list_pnts <- list_pnts_new
    #
    # update list_pnts
    list_pntu <- list_pntu_new
    #
    # update list_pnts
    list_grpm <- list_grpm_new
    #
    # list2pnts
    pnts_split <- f_list2pnts(list_pnts)
  }
  #
  # return
  return(spldf_profile)
}
#
#
#
f_polygon_rasterize <- function(raster_module, shp_pgn) {
  #
  # # plot raster_module
  # plot(raster_module, col = rev(terrain.colors(50)))
  #
  # # plot shapefile
  # # notice that you use add = T to add a layer on top of an existing plot in R.
  # plot(shp_pgn, main = "polygon", axes = TRUE, border = "blue")
  #
  # must first crop
  # OR, for large raster file
  # very large temporary raster files will be saved in C disk
  # C disk will become full quickly
  raster_module_crop <- raster::crop(raster_module, shp_pgn)
  #
  # # solve "Error in plot.new() : figure margins too large"
  # par(mar = c(1, 1, 1, 1))
  # #
  # # plot
  # plot(raster_module_crop, main = "cropped raster")
  #
  # mask the raster using the polygon
  raster_module_crop_mask <- raster::mask(raster_module_crop, shp_pgn)
  #
  # # plot
  # plot(raster_module_crop_mask, main = "masked raster")
  # plot(shp_pgn, add = TRUE)
  #
  # set one
  raster_module_crop_mask_set1 <- raster::calc(raster_module_crop_mask, function(x) x / x)
  #
  # # plot
  # plot(raster_module_crop_mask_set1, main = "masked raster set1")
  # plot(shp_pgn, add = TRUE)
  #
  # raster to polygon
  raster_module_crop_mask_set1_pgn <- raster::rasterToPolygons(raster_module_crop_mask_set1, fun = NULL, na.rm = TRUE, dissolve = TRUE)
  #
  # # plot shapefile
  # # notice that you use add = T to add a layer on top of an existing plot in R.
  # plot(raster_module_crop_mask_set1_pgn, main = "polygon set1", axes = TRUE, border = "blue")
  # #
  # # write shapefile
  # rgdal::writeOGR(raster_module_crop_mask_set1_pgn, ".", "raster_module_crop_mask_set1_pgn", driver = "ESRI Shapefile")
  #
  # return
  # return(raster_module_crop_mask_set1_pgn)
  return(list(raster_module_crop_mask_set1_pgn, raster_module_crop_mask))
}
#
#
#
f_pnts_read <- function(raster_dem, shp_lasld, file_out = "lasld_pln") {
  #
  # mask the raster using the vector extent
  raster_dem_mask <- raster::mask(raster_dem, shp_lasld)
  #
  # get pnts from raster
  pnts_from_raster <- raster::rasterToPoints(raster_dem_mask)
  #
  # get pnts
  pnts <- pnts_from_raster
  #
  # colnames
  colnames(pnts) <- c("Cx", "Cy", "Cz")
  #
  # get count
  count_pnt <- nrow(pnts)
  #
  # get line from polygon
  shp_lasld_pln <- as(shp_lasld, 'SpatialLinesDataFrame')
  # rgdal::writeOGR(shp_lasld_pln, ".", "shp_lasld_pln", driver = "ESRI Shapefile")
  #
  # initial
  column_IDB <- c()
  #
  # maybe more than 1 (polygon) lines, difficult to handle
  if (length(shp_lasld_pln@lines) > 1) {
    #
    # warning
    warning("Error: shp_lasld_pln has more than 1 lines, read pnts failed.")
    #
    # save
    rgdal::writeOGR(shp_lasld_pln, ".", file_out, overwrite_layer = T, driver = "ESRI Shapefile")
    #
    # get column_IDB
    column_IDB <- matrix(0, count_pnt, 1)
  }
  else if (length(shp_lasld_pln@lines[[1]]@Lines) > 1) {
    #
    # warning
    warning("Error: lines in shp_lasld_pln has more than 1 lines, read pnts failed.")
    #
    # save
    rgdal::writeOGR(shp_lasld_pln, ".", file_out, overwrite_layer = T, driver = "ESRI Shapefile")
    #
    # get column_IDB
    column_IDB <- matrix(0, count_pnt, 1)
  }
  else {
    #
    # get coords of line from polygon
    shp_lasld_pln_coords <- shp_lasld_pln@lines[[1]]@Lines[[1]]@coords
    #
    # get count
    count_coords <- nrow(shp_lasld_pln_coords)
    #
    # get cell size
    raster_cellsize <- (raster::res(raster_dem))[1]
    #
    # initial column_IDB, IDB (ID of boundary pnts)
    column_IDB <- matrix(0, count_pnt, 1)
    #
    # initial IDB
    IDB <- 0
    #
    # get column_IDB
    for (i in 1:(count_coords-1)) {
      #
      # get pnt, the first and the last are the same
      pnt1 <- shp_lasld_pln_coords[i, ]
      pnt2 <- shp_lasld_pln_coords[i+1, ]
      #
      # get pnt_IDB
      pnt_IDB <- f_pnt_find_IDB(pnt1, pnt2, pnts, raster_cellsize)
      #
      # debugging
      if (is.na(pnt_IDB)) {
        #
        # warning
        warning("Error: pnt_IDB is NA, read pnts failed.")
        #
        # save
        rgdal::writeOGR(shp_lasld_pln, ".", file_out, overwrite_layer = T, driver = "ESRI Shapefile")
        #
        # get column_IDB
        column_IDB <- matrix(0, count_pnt, 1)
        #
        # break
        break
      }
      #
      # assign IDB
      if (column_IDB[pnt_IDB] == 0) { 
        #
        # assign IDB
        IDB <- IDB + 1
        column_IDB[pnt_IDB] = IDB
      }
    }
  }
  #
  # add column
  column_ID <- c(1:count_pnt)
  pnts <- cbind(column_ID, column_IDB, pnts)
  #
  # add column
  column_grp <- matrix(1, count_pnt, 1)
  pnts <- cbind(pnts, column_grp)
  #
  # add columns of zeros
  pnts <- cbind(pnts, matrix(0, nrow(pnts), 7))
  #
  # colnames
  colnames(pnts) <- c("ID", "IDB", "Cx", "Cy", "Cz", "grp", "Nx", "Ny", "Nz", "zip", "Dx", "Dy", "Dz")
  #
  # set data frame
  pnts <- data.frame(pnts)
  #
  # add columns of SMG and SMS (stop marker)
  pnts_SM <- cbind(matrix("SMS", nrow(pnts), 1), matrix("NA", nrow(pnts), 1))
  colnames(pnts_SM) <- c("SMG", "SMS")
  pnts_SM <- data.frame(pnts_SM)
  pnts_SM[, "SMG"] <- as.character(pnts_SM[, "SMG"])
  pnts_SM[, "SMS"] <- as.character(pnts_SM[, "SMS"])
  #
  # get pnts
  pnts <- cbind(pnts, pnts_SM)
  #
  # get plane
  pnts <- f_pnts_plane(pnts)
  #
  # return
  return(pnts)
}
#
#
#
f_pnt_find_IDB <- function(pnt1, pnt2, pnts, raster_cellsize) {
  #
  # get coords
  x1 <- pnt1[1]
  y1 <- pnt1[2]
  x2 <- pnt2[1]
  y2 <- pnt2[2]
  #
  # if horizontal
  if (y1 == y2) {
    #
    # get bnd of x
    xmin <- (x1 + x2) / 2 - raster_cellsize / 4
    xmax <- (x1 + x2) / 2 + raster_cellsize / 4
    #
    # get bnd of y
    ymin <- y1 - raster_cellsize / 4 * 3
    ymax <- y1 + raster_cellsize / 4 * 3
    #
    # get IDB
    pnt_IDB <- f_pnt_find_ID(c(xmin, xmax), c(ymin, ymax), pnts)
  }
  #
  # if vertical
  else {
    #
    # get bnd of x
    xmin <- x1 - raster_cellsize / 4 * 3
    xmax <- x1 + raster_cellsize / 4 * 3
    #
    # get bnd of y
    ymin <- (y1 + y2) / 2 - raster_cellsize / 4
    ymax <- (y1 + y2) / 2 + raster_cellsize / 4
    #
    # get IDB
    pnt_IDB <- f_pnt_find_ID(c(xmin, xmax), c(ymin, ymax), pnts)
  }
  #
  # return
  return(pnt_IDB)
}
#
#
#
f_pnt_find_ID <- function(bndx, bndy, pnts) {
  #
  # get count
  count_pnt <- nrow(pnts)
  #
  # get bnds
  xmin <- bndx[1]
  xmax <- bndx[2]
  ymin <- bndy[1]
  ymax <- bndy[2]
  #
  # get 
  for (i in 1:count_pnt) {
    #
    # get pnt
    pntx <- pnts[i, "Cx"]
    pnty <- pnts[i, "Cy"]
    #
    # if pnt in bnds
    if (pntx > xmin && pntx < xmax && pnty > ymin && pnty < ymax) {
      #
      # return
      return(i)
    }
  }
  #
  # return
  return(NA)
}
#
#
#
f_pnts_plane <- function(pnts) {
  #
  # get lm
  pnts_lm <- lm(Cz ~ Cx + Cy, data = data.frame(pnts))
  pnts_lm_coeffs <- pnts_lm[[1]]
  #
  # get normal
  # pnts_normal <- c(pnts_lm_coeffs[2], pnts_lm_coeffs[3], -1)
  pnts_normal <- c(-pnts_lm_coeffs[2], -pnts_lm_coeffs[3], 1)
  pnts_normal <- f_vector_unit(pnts_normal)
  #
  # debugging
  if (is.na(pnts_normal[1]) || is.na(pnts_normal[2]) || is.na(pnts_normal[3])) {
    #
    # return
    return(NA)
  }
  #
  # get normal
  pnts[, "Nx"] <- data.frame(matrix(pnts_normal[1], nrow(pnts), 1))
  pnts[, "Ny"] <- data.frame(matrix(pnts_normal[2], nrow(pnts), 1))
  pnts[, "Nz"] <- data.frame(matrix(pnts_normal[3], nrow(pnts), 1))
  #
  # get zip (z in plane)
  pnts[, "zip"] <-  pnts[, "Cx"] * pnts_lm_coeffs[2] + pnts[, "Cy"] * pnts_lm_coeffs[3] + pnts_lm_coeffs[1]
  #
  # get dip
  pnts_dip <- f_pnts_dip(pnts)
  pnts[, "Dx"] <- data.frame(matrix(pnts_dip[1], nrow(pnts), 1))
  pnts[, "Dy"] <- data.frame(matrix(pnts_dip[2], nrow(pnts), 1))
  pnts[, "Dz"] <- data.frame(matrix(pnts_dip[3], nrow(pnts), 1))
  #
  # return
  return(pnts)
}
#
#
#
f_pnts_dip <- function(pnts) {
  #
  # get normal
  pnts_normal <- c(pnts[1, "Nx"], pnts[1, "Ny"], pnts[1, "Nz"])
  #
  # debugging
  if (is.na(pnts_normal[1]) || is.na(pnts_normal[2]) || is.na(pnts_normal[3])) {
    #
    # return
    return(NA)
  }
  #
  # get plumb
  vector_plumb <- c(0, 0, -1)
  #
  # get strike
  pnts_strike <- pracma::cross(pnts_normal, vector_plumb)
  pnts_strike <- f_vector_unit(pnts_strike)
  #
  # vertical normal (horizontal plane)
  if (pnts_strike[1] ==0 && pnts_strike[2] == 0 && pnts_strike[3] ==0) { 
    #
    # waring
    warning("Warning: A fitted plane is horizontal.")
    #
    # use the dip of mother plane
    return(c(pnts[1, "Dx"], pnts[1, "Dy"], pnts[1, "Dz"])) 
  }
  #
  # get dip
  pnts_dip <- pracma::cross(pnts_strike, pnts_normal)
  pnts_dip <- f_vector_unit(pnts_dip)
  # 
  # return
  return(pnts_dip)
}
#
#
#
f_vector_unit <- function(vector) {
  #
  # get magnitude
  Magnitue_vector <- sqrt(vector[1] * vector[1] + vector[2] * vector[2] + vector[3] * vector[3])
  #
  # get normalized
  vector <- vector / Magnitue_vector
  # 
  # return
  return(vector)
}
#
#
#
f_pnts_intercept <- function(pnts) {
  #
  # get count
  count_pnt <- nrow(pnts)
  #
  # get normal of cut plane (pnts_dip)
  pnts_dip <- c(pnts[1, "Dx"], pnts[1, "Dy"], pnts[1, "Dz"])
  #
  # initial column
  column_intercept <- c()
  #
  # get intercepts
  for (i in 1:count_pnt) {
    #
    # get pnt
    pnt <- c(pnts[i, "Cx"], pnts[i, "Cy"], pnts[i, "Cz"])
    #
    # get intercepts
    column_intercept <- rbind(column_intercept, f_pnt_intercept(pnt, pnts_dip))
  }
  #
  # return
  return((min(column_intercept) + max(column_intercept)) / 2)
}
#
#
#
f_pnt_intercept <- function(pnt, normal) {
  #
  # get intercept
  intercept <- -(normal[1] * pnt[1] + normal[2] * pnt[2] + normal[3] * pnt[3])
  #
  # return
  return(intercept)
}
#
#
#
f_pnts_split <- function(pnts, pnts_intercept) {
  #
  # get dip
  pnts_dip <- c(pnts[1, "Dx"], pnts[1, "Dy"], pnts[1, "Dz"])
  #
  # new
  pnts1 <- c()
  pnts2 <- c()
  #
  # get count
  count_pnt <- nrow(pnts)
  #
  # get split
  for (i in 1:count_pnt) {
    #
    # get direction
    direction = pnts_dip[1] * pnts[i, "Cx"] + pnts_dip[2] * pnts[i, "Cy"] + pnts_dip[3] * pnts[i, "Cz"] + pnts_intercept
    #
    # split
    if (direction <= 0) { pnts1 <- rbind(pnts1, pnts[i, ]) }
    else { pnts2 <- rbind(pnts2, pnts[i, ]) }
  }
  #
  # return
  return(list(pnts1, pnts2))
}
#
#
#
f_pnts_clustering <- function(pnts, raster_cellsize, raster_CRS) {
  #
  # initial list
  list_clusters <- list()
  #
  # get polygons
  pnts_polygon <- f_pnts2shp_polygon_rasterized(pnts, raster_cellsize, raster_CRS)
  #
  # convert
  pnts_polygon_sf <- as(pnts_polygon, "sf")
  #
  # multiparts to singleparts
  pnts_polygon_sf_sparts <- sf::st_cast(pnts_polygon_sf, "POLYGON")
  #
  # convert
  pnts_polygon_sf_sparts_sp <- as(pnts_polygon_sf_sparts, "Spatial")
  #
  # ID of polygon
  pnts_polygon_sf_sparts_sp_ID <- 0
  #
  # debugging
  if (F) {
    #
    # pause
    if (length(pnts_polygon_sf_sparts_sp@polygons) > 1) {
      #
      # warning
      warning("pause.")
    }
  }
  #
  # list polygons
  for (i in 1:length(pnts_polygon_sf_sparts_sp@polygons)) {
    #
    # list Polygons
    for (j in 1:length(pnts_polygon_sf_sparts_sp@polygons[[i]]@Polygons)) {
      #
      # get ID
      pnts_polygon_sf_sparts_sp_ID <- pnts_polygon_sf_sparts_sp_ID + 1
      #
      # get Polygon
      pnts_polygon_sf_sparts_sp_pgn <- pnts_polygon_sf_sparts_sp[pnts_polygon_sf_sparts_sp_ID, ]
      #
      # debugging
      if (F) {
        #
        # plot polygon
        plot(pnts_polygon_sf_sparts_sp_pgn, main = "polygon", axes = TRUE, border = "red", add = TRUE)
      }
      #
      # get pnts shp
      pnts_shp <- f_pnts2shp_points(pnts, raster_CRS)
      #
      # get over
      pnts_shp_over <- sp::over(pnts_shp, pnts_polygon_sf_sparts_sp_pgn)
      #
      # update
      list_clusters[[length(list_clusters) + 1]] <- pnts[!is.na(pnts_shp_over), ]
    }
  }
  #
  # return
  return(list_clusters)
}
#
#
#
f_pnts_bndcorners_check <- function(pnts, IDB_max, pnts_initial_bnd, pnts_distal_bnd, boolEND) {
  #
  # get pnts_bnd
  pnts_bnd <- pnts[which(pnts[, "IDB"] != 0), ]
  #
  # get bndcorners
  bc_pnts_bnd <- f_pnts_bndcorners(pnts_bnd, IDB_max)
  #
  # check spatial impartiality
  if (is.null(bc_pnts_bnd)) {
    #
    # spatial impartiality
    return(1)
  }
  #
  # can not be 1???
  if (nrow(bc_pnts_bnd) == 1) {
    #
    # warning
    warning("Error: bndcorners count of initial or distal pnts is 1.")
    #
    # return
    return(NA)
  }
  #
  # if initial or distal pnts, check
  if (boolEND == TRUE) {
    #
    # ok or not
    if (nrow(bc_pnts_bnd) > 2) {
      #
      # not ok for end pnts
      return(3)
    }
    else if (nrow(bc_pnts_bnd) == 2) {
      #
      # ok
      return(0)
    }
    else {
      #
      # warning
      warning("Error: bndcorners count of initial or distal pnts is 1.")
      return(NA)
    }
  }
  #
  # for debugging
  pnts_initial_bnd_input <- pnts_initial_bnd
  pnts_distal_bnd_input <- pnts_distal_bnd
  #
  # get count
  count_pnt <- nrow(pnts_bnd)
  #
  # remove duplicates
  # end pnts excluding the pnts being checked
  # it is possible that the pnts being checked is splitted from end pnts
  for (i in 1:count_pnt) {
    #
    # remove duplicates
    ID_start <- which(pnts_initial_bnd[, "ID"] == pnts_bnd[i, "ID"])
    if (length(ID_start) != 0) {
      #
      pnts_initial_bnd <- pnts_initial_bnd[-ID_start, ]
    }
    #
    # remove duplicates
    ID_end <- which(pnts_distal_bnd[, "ID"] == pnts_bnd[i, "ID"])
    if (length(ID_end) != 0) {
      #
      pnts_distal_bnd <- pnts_distal_bnd[-ID_end, ]
    }
  }
  #
  # get bndcorners of initial and distal pnts
  bc_pnts_initial_bnd <- f_pnts_bndcorners(pnts_initial_bnd, IDB_max)
  bc_pnts_distal_bnd <- f_pnts_bndcorners(pnts_distal_bnd, IDB_max)
  #
  # check end pnts after splitting
  if (nrow(bc_pnts_initial_bnd) > 2 || nrow(bc_pnts_distal_bnd) > 2) {
    #
    # not ok for end pnts
    return(3)
  }
  else if (nrow(bc_pnts_initial_bnd) < 2 || nrow(bc_pnts_distal_bnd) < 2) {
    #
    # warning
    warning("Error: bndcorners count of initial or distal pnts is 1.")
    return(NA)
  }
  #
  # get IDBs of initial and distal pnts
  IDB_initial1 <- bc_pnts_initial_bnd[1, "IDB"]
  IDB_initial2 <- bc_pnts_initial_bnd[2, "IDB"]
  IDB_distal1 <- bc_pnts_distal_bnd[1, "IDB"]
  IDB_distal2 <- bc_pnts_distal_bnd[2, "IDB"]
  #
  # circle from initial1
  IDB_initial1_initial2 <- IDB_initial2 - IDB_initial1
  if (IDB_initial1_initial2 < 0) { IDB_initial1_initial2 <- IDB_initial1_initial2 + IDB_max }
  IDB_initial1_distal1 <- IDB_distal1 - IDB_initial1
  if (IDB_initial1_distal1 < 0) { IDB_initial1_distal1 <- IDB_initial1_distal1 + IDB_max }
  #
  # if get initial2 first
  if (IDB_initial1_initial2 < IDB_initial1_distal1) {
    #
    # replace
    IDB_tmp <- IDB_initial1
    IDB_initial1 <- IDB_initial2
    IDB_initial2 <- IDB_tmp
  }
  #
  # circle from initial1
  IDB_initial1_distal1 <- IDB_distal1 - IDB_initial1
  if (IDB_initial1_distal1 < 0) { IDB_initial1_distal1 <- IDB_initial1_distal1 + IDB_max }
  IDB_initial1_distal2 <- IDB_distal2 - IDB_initial1
  if (IDB_initial1_distal2 < 0) { IDB_initial1_distal2 <- IDB_initial1_distal2 + IDB_max }
  #
  # if get distal2 first
  if (IDB_initial1_distal2 < IDB_initial1_distal1) {
    #
    # replace
    IDB_tmp <- IDB_distal1
    IDB_distal1 <- IDB_distal2
    IDB_distal2 <- IDB_tmp
  }
  #
  # circle from initial1
  IDB_initial1_initial2 <- IDB_initial2 - IDB_initial1
  if (IDB_initial1_initial2 < 0) { IDB_initial1_initial2 <- IDB_initial1_initial2 + IDB_max }
  IDB_initial1_distal1 <- IDB_distal1 - IDB_initial1
  if (IDB_initial1_distal1 < 0) { IDB_initial1_distal1 <- IDB_initial1_distal1 + IDB_max }
  IDB_initial1_distal2 <- IDB_distal2 - IDB_initial1
  if (IDB_initial1_distal2 < 0) { IDB_initial1_distal2 <- IDB_initial1_distal2 + IDB_max }
  #
  # initial
  count_bc_pnts_bnd_0 <- 0
  count_bc_pnts_bnd_1 <- 0
  count_bc_pnts_bnd_2 <- 0
  #
  # get count
  count_bc_pnts_bnd <- nrow(bc_pnts_bnd)
  #
  # get sides
  for (i in 1:count_bc_pnts_bnd) {
    #
    # get IDB
    IDB_bc <- bc_pnts_bnd[i, "IDB"]
    #
    # circle from initial1
    IDB_initial1_bc <- IDB_bc - IDB_initial1
    if (IDB_initial1_bc < 0) { IDB_initial1_bc <- IDB_initial1_bc + IDB_max }
    #
    # check sides
    if (IDB_initial1_bc < IDB_initial1_distal1) {
      #
      # 1st side
      count_bc_pnts_bnd_1 <- count_bc_pnts_bnd_1 + 1
    }
    else if (IDB_initial1_bc > IDB_initial1_distal2 && IDB_initial1_bc < IDB_initial1_initial2) {
      #
      # 2nd side
      count_bc_pnts_bnd_2 <- count_bc_pnts_bnd_2 + 1
    }
    else {
      #
      # neither side
      count_bc_pnts_bnd_0 <- count_bc_pnts_bnd_0 + 1
    }
  }
  #
  # check sides
  if (count_bc_pnts_bnd_0 != 0) {
    #
    # warning
    warning("Error: bndcorners do not belong to either side.")
    return(NA)
  }
  #
  # check sides
  if (count_bc_pnts_bnd_1 == 0 || count_bc_pnts_bnd_2 == 0) {
    #
    # spatial impartiality
    return(1)
  }
  #
  # check sides
  if (count_bc_pnts_bnd_1 > 2 || count_bc_pnts_bnd_2 > 2) {
    #
    # spatial discontiguity
    return(2)
  }
  #
  # check sides
  if (count_bc_pnts_bnd_1 == 2 && count_bc_pnts_bnd_2 == 2) {
    #
    # ok
    return(0)
  }
  #
  # warning
  warning("Error: bndcorners count is not right.")
  return(NA)
}
#
#
#
f_pnts_bndcorners <- function(pnts_bnd, IDB_max) {
  #
  # the initial and distal pnts only have two bnd corners?
  #
  # get count
  count_pnt <- nrow(pnts_bnd)
  #
  # on bnd
  if (count_pnt == 0) { return(NULL) }
  #
  # only one bnd
  if (count_pnt == 1) {
    #
    # return
    return(rbind(pnts_bnd, pnts_bnd))
  }
  #
  # get sorted
  pnts_bnd_ascending <- pnts_bnd[order(pnts_bnd[, "IDB"]), ]
  #
  # initial
  bndcorners <- c()
  #
  # get corners
  for (i in 1:(count_pnt-1)) {
    #
    # get IDB step
    IDB_step <- pnts_bnd_ascending[i+1, "IDB"] - pnts_bnd_ascending[i, "IDB"]
    #
    # debugging
    if (is.na(IDB_step)) {
      #
      # warning
      warning("IDB_step is NA")
    }
    #
    # adjust IDB step
    if (IDB_step < 0) { IDB_step <- IDB_step + IDB_max }
    #
    # bndcorners
    if (IDB_step != 1) {
      #
      # append
      bndcorners <- rbind(bndcorners, pnts_bnd_ascending[i, ], pnts_bnd_ascending[i+1, ])
    }
  }
  #
  # get IDB step
  IDB_step <- pnts_bnd_ascending[1, "IDB"] - pnts_bnd_ascending[count_pnt, "IDB"]
  #
  # debugging
  if (is.na(IDB_step)) {
    #
    # warning
    warning("IDB_step is NA")
  }
  #
  # adjust IDB step
  if (IDB_step < 0) { IDB_step <- IDB_step + IDB_max }
  #
  # bndcorners
  if (IDB_step != 1) {
    #
    # append
    bndcorners <- rbind(bndcorners, pnts_bnd_ascending[1, ], pnts_bnd_ascending[count_pnt, ])
  }
  #
  # debugging
  if (is.null(bndcorners) || is.na(bndcorners)) {
    #
    # warning
    warning("bndcorners is Null or NA.")
  }
  #
  # return
  return(bndcorners)
}
#
#
#
f_pnts_MinDist_intercept <- function(pnts, pnts_dip, pnts_intercept) {
  #
  # get count
  count_pnt <- nrow(pnts)
  # #
  # # get normal of cut plane (pnts_dip)
  # pnts_dip <- c(pnts[1, "Dx"], pnts[1, "Dy"], pnts[1, "Dz"])
  #
  # initial
  ID_MinDiffIntercept <- NA
  MinDiffIntercept <- Inf
  #
  # get intercept
  for (i in 1:count_pnt) {
    #
    # get pnt
    pnt <- c(pnts[i, "Cx"], pnts[i, "Cy"], pnts[i, "Cz"])
    #
    # get intercept
    pnt_intercept <- f_pnt_intercept(pnt, pnts_dip)
    #
    # get DiffIntercept
    DiffIntercept <- abs(pnt_intercept - pnts_intercept)
    #
    # find the nearest
    if (DiffIntercept < MinDiffIntercept) {
      #
      # update i
      ID_MinDiffIntercept <- i
      MinDiffIntercept <- DiffIntercept
    }
  }
  #
  # return
  return(pnts[ID_MinDiffIntercept, ])
}
#
#
#
f_pnts_anchors <- function(list_pnts, IDB_max) {
  #
  # get count
  count_pnts <- length(list_pnts)
  #
  # if only 1 pnts, i.e. the initial elevation point group
  if (count_pnts == 1) {
    #
    # get pnts
    pnts <- list_pnts[[1]]
    #
    # get anchor initial, the point with maximum Cz
    anchor_initial <- pnts[which.max(pnts[, "Cz"]), c("Cx", "Cy", "Cz")]
    #
    # get anchors
    list_anchors <- f_pnts_anchor_end(pnts, anchor_initial[, c("Cx", "Cy")])
    #
    anchor_distal <- list_anchors[[1]]
    anchor_mid <- list_anchors[[2]]
    #
    # get type
    column_type <- data.frame("Initial group anchor", "Group center", "Distal group anchor")
    column_type <- t( as.data.frame(column_type, col.names = F) )
    #
    # get grp
    column_grp <- data.frame(0, 1, 0)
    column_grp <- t( as.data.frame(column_grp, col.names = F) )
    #
    # get column_anchors
    column_anchors <- rbind(anchor_initial, anchor_mid, anchor_distal)
    column_anchors <- data.frame(column_anchors)
    #
    # cbind
    column_anchors <- cbind(column_anchors, column_type, column_grp)
    column_anchors <- data.frame(column_anchors)
    #
    # colnames
    colnames(column_anchors) <- c("Cx", "Cy", "Cz", "type", "grp")
    #
    # return
    return(column_anchors)
  }
  #
  # initial
  column_centerL <- c()
  #
  # for initial pnts
  pnts_initial <- list_pnts[[1]]
  #
  # get center of split "line"
  centerL_initial <- f_pnts_centerL_end(pnts_initial, IDB_max)
  column_centerL <- rbind(centerL_initial, column_centerL)
  #
  # only if more than 2
  if (count_pnts >= 3) {
    #
    # get centers
    for (i in 2:(count_pnts-1)) {
      #
      # get pnts_i
      pnts_i <- list_pnts[[i]]
      #
      # get pnts_i_bnd
      pnts_i_bnd <- pnts_i[which(pnts_i[, "IDB"] != 0), ]
      #
      # get bnd corners()
      pnts_i_bndcorners <- f_pnts_bndcorners(pnts_i_bnd, IDB_max)
      #
      # get center
      if (nrow(pnts_i_bndcorners) == 4) {
        #
        # get center xy
        centerP_x1 <- (pnts_i_bndcorners[2, "Cx"] + pnts_i_bndcorners[3, "Cx"]) / 2
        centerP_y1 <- (pnts_i_bndcorners[2, "Cy"] + pnts_i_bndcorners[3, "Cy"]) / 2
        centerP_x2 <- (pnts_i_bndcorners[4, "Cx"] + pnts_i_bndcorners[1, "Cx"]) / 2
        centerP_y2 <- (pnts_i_bndcorners[4, "Cy"] + pnts_i_bndcorners[1, "Cy"]) / 2
        #
        # get center xyz
        centerL_x1 <- (pnts_i_bndcorners[1, "Cx"] + pnts_i_bndcorners[2, "Cx"]) / 2
        centerL_y1 <- (pnts_i_bndcorners[1, "Cy"] + pnts_i_bndcorners[2, "Cy"]) / 2
        centerL_z1 <-  f_pnts_zip(pnts_i, centerL_x1, centerL_y1)
        centerL1 <- c(centerL_x1, centerL_y1, centerL_z1)
        #
        # get center xyz
        centerL_x2 <- (pnts_i_bndcorners[3, "Cx"] + pnts_i_bndcorners[4, "Cx"]) / 2
        centerL_y2 <- (pnts_i_bndcorners[3, "Cy"] + pnts_i_bndcorners[4, "Cy"]) / 2
        centerL_z2 <-  f_pnts_zip(pnts_i, centerL_x2, centerL_y2)
        centerL2 <- c(centerL_x2, centerL_y2, centerL_z2)
        #
        # the last centerL in the stack
        centerL_last_x <- column_centerL[nrow(column_centerL), 1]
        centerL_last_y <- column_centerL[nrow(column_centerL), 2]
        #
        # get distances
        dist1 <- sqrt((centerL_x1 - centerL_last_x)^2 + (centerL_y1 - centerL_last_y)^2)
        dist2 <- sqrt((centerL_x2 - centerL_last_x)^2 + (centerL_y2 - centerL_last_y)^2)
        #
        # get center of line
        if (dist1 > dist2) { column_centerL <- rbind(column_centerL, centerL2, centerL1) }
        else { column_centerL <- rbind(column_centerL, centerL1, centerL2) }
      }
      else {
        #
        # error
        warning("Error: pnts do not have 4 bnd corners")
        #
        # return
        return(NULL)
      }
    }
  }
  #
  # for distal pnts
  pnts_distal <- list_pnts[[count_pnts]]
  #
  # get center of split "line"
  centerL_distal <- f_pnts_centerL_end(pnts_distal, IDB_max)
  column_centerL <- rbind(column_centerL, centerL_distal)
  # #
  # # output column_centerL for show
  # file_out <- paste("column_centerL_grp", as.character(pnts_distal[nrow(pnts_distal), "grp"]), ".xlsx", sep = "", collapse = NULL)
  # f_pnts_save_xls(column_centerL, file_out)
  #
  # initial
  column_centerG <- c()
  #
  # initial type
  column_type <- c()
  #
  # initial grp
  column_grp <- c()
  #
  # get count
  count_centerL <- nrow(column_centerL)
  #
  # update centerG
  for (i in 1:(count_centerL-1)) {
    #
    # get xyz
    centerL1x <- column_centerL[i, 1]
    centerL1y <- column_centerL[i, 2]
    centerL1z <- column_centerL[i, 3]
    # get xyz
    centerL2x <- column_centerL[i+1, 1]
    centerL2y <- column_centerL[i+1, 2]
    centerL2z <- column_centerL[i+1, 3]
    # get xyz
    centerLx <- (centerL1x + centerL2x) / 2
    centerLy <- (centerL1y + centerL2y) / 2
    centerLz <- (centerL1z + centerL2z) / 2
    #
    # replace
    column_centerG <- rbind(column_centerG, c(centerLx, centerLy, centerLz))
    #
    # type and grp
    if (i%%2 == 1) {
      #
      # odd
      column_type <- rbind(column_type, "Inter-group center")
      column_grp <- rbind(column_grp, 0)
    }
    else {
      #
      # even
      column_type <- rbind(column_type, "Group center")
      column_grp <- rbind(column_grp, i/2 + 1)
    }
  }
  #
  # must do this, or numeric will become character
  column_centerG <- data.frame(column_centerG)
  column_type <- data.frame(column_type)
  column_grp <- data.frame(column_grp)
  #
  # initial
  column_anchors <- cbind(column_centerG, column_type, column_grp)
  column_anchors <- data.frame(column_anchors)
  #
  # colnames
  # colnames(column_anchors) <- c("Cx", "Cy", "zip")
  # colnames(column_anchors) <- c("Cx", "Cy", "Cz")
  colnames(column_anchors) <- c("Cx", "Cy", "Cz", "type", "grp")
  #
  # get anchors, for initial
  list_anchors <- f_pnts_anchor_end(pnts_initial, c(column_anchors[1, 1], column_anchors[1, 2]))
  anchor_end_initial_edge <- list_anchors[[1]]
  #
  # is possible that f_pnts_anchor_end return Null?
  if (is.null(anchor_end_initial_edge)) {
    #
    # return
    return(NULL)
  }
  #
  # get anchor_end_intial
  anchor_end_initial_edge <- cbind(anchor_end_initial_edge, "Initial group anchor", 0)
  colnames(anchor_end_initial_edge) <- c("Cx", "Cy", "Cz", "type", "grp")
  anchor_end_initial_center <- cbind(list_anchors[[2]], "Initial group center", 1)
  colnames(anchor_end_initial_center) <- c("Cx", "Cy", "Cz", "type", "grp")
  #
  # update anchors
  column_anchors <- rbind(anchor_end_initial_edge, anchor_end_initial_center, column_anchors)
  #
  # get anchors, for distal
  list_anchors <- f_pnts_anchor_end(pnts_distal, c(column_anchors[nrow(column_anchors), 1], column_anchors[nrow(column_anchors), 2]))
  anchor_end_distal_edge <- list_anchors[[1]]
  #
  # is possible that f_pnts_anchor_end return NULL?
  if (is.null(anchor_end_distal_edge)) {
    #
    # return
    return(NULL)
  }
  #
  # get anchor_end_distal
  anchor_end_distal_edge <- cbind(anchor_end_distal_edge, "Distal group anchor", 0)
  colnames(anchor_end_distal_edge) <- c("Cx", "Cy", "Cz", "type", "grp")
  anchor_end_distal_center <- cbind(list_anchors[[2]], "Distal group center", (list_pnts[[count_pnts]])[1, "grp"])
  colnames(anchor_end_distal_center) <- c("Cx", "Cy", "Cz", "type", "grp")
  #
  # update anchors
  column_anchors <- rbind(column_anchors, anchor_end_distal_center, anchor_end_distal_edge)
  # #
  # # for debugging
  # file_out <- paste("column_anchors.xlsx", sep = "", collapse = NULL)
  # f_pnts_save_xls(column_anchors, file_out)
  # 
  # return
  return(data.frame(column_anchors))
}
#
#
#
f_pnts_centerL_end <- function(pnts, IDB_max) {
  #
  # get pnts_bnd
  pnts_bnd <- pnts[which(pnts[, "IDB"] != 0), ]
  #
  # get bndcorners
  pnts_bndcorners <- f_pnts_bndcorners(pnts_bnd, IDB_max)
  #
  # get center
  if (nrow(pnts_bndcorners) == 2) {
    #
    # get center xy
    center_x <- (pnts_bndcorners[1, "Cx"] + pnts_bndcorners[2, "Cx"]) / 2
    center_y <- (pnts_bndcorners[1, "Cy"] + pnts_bndcorners[2, "Cy"]) / 2
    #
    # get center z (z in plane)
    center_z <-  f_pnts_zip(pnts, center_x, center_y)
    #
    # get center
    center <- c(center_x, center_y, center_z)
    #
    # return
    return(center)
  }
  else {
    #
    # error
    warning("Error: (end) pnts do not have 2 bnd corners")
    #
    # return
    return(NA)
  }
}
#
#
#
f_pnts_anchor_end <- function(pnts, lpnt0) {
  #
  # for the two end pnts (inital and distal)
  # and also
  # for the initial landslide elevation points
  #
  # get anchor end
  # the elevation of anchor end is the elevation of the nearest bnd pnt
  # not calculated using the fitted plane
  anchor_end <- f_pnts_split_equal_anchor(pnts, lpnt0)
  #
  # get x and y, anchor middle
  anchor_mid_x <- (lpnt0[1] + anchor_end[1]) / 2
  anchor_mid_y <- (lpnt0[2] + anchor_end[2]) / 2
  #
  # get z
  anchor_mid_z <- f_pnts_zip(pnts, as.numeric(anchor_mid_x), as.numeric(anchor_mid_y))
  #
  # get anchor middle
  anchor_mid <- data.frame(cbind(anchor_mid_x, anchor_mid_y, anchor_mid_z))
  colnames(anchor_mid) <- c("Cx", "Cy", "Cz")
  #
  # return
  return(list(anchor_end, anchor_mid))
}
#
#
#
f_pnts_split_equal_anchor <- function(pnts, lpnt0) {
  #
  # initial
  split_equal_angle <- f_pnts_split_equal(pnts, lpnt0)
  #
  # get line
  lX <- 1 * (pi/2 - abs(split_equal_angle)) / abs(pi/2 - abs(split_equal_angle))
  lY <- lX * tan(split_equal_angle)
  if (pi/2 - abs(split_equal_angle) == 0) { lX <- 0 }
  #
  # get limits
  Xmin <- pnts[which.min(pnts[, "Cx"]), "Cx"]
  Xmax <- pnts[which.max(pnts[, "Cx"]), "Cx"]
  Ymin <- pnts[which.min(pnts[, "Cy"]), "Cy"]
  Ymax <- pnts[which.max(pnts[, "Cy"]), "Cy"]
  #
  # get coords
  Xmin_Y <- lpnt0[2] + (Xmin - lpnt0[1]) * lY / lX
  Xmax_Y <- lpnt0[2] + (Xmax - lpnt0[1]) * lY / lX
  #
  # get pnts_curve
  pnts_curve_segment <- rbind(as.numeric(c(Xmin, Xmin_Y)), as.numeric(c(Xmax, Xmax_Y)))
  #
  # if vertical
  if (lX == 0) { pnts_curve_segment <- rbind(as.numeric(c(lpnt0[1], Ymin)), as.numeric(c(lpnt0[1], Ymax))) }
  #
  # get
  colnames(pnts_curve_segment) <- c("Cx", "Cy")
  pnts_curve_segment <- data.frame(pnts_curve_segment)
  #
  # get pnts_bnd
  pnts_bnd <- pnts[which(pnts[, "IDB"] != 0), ]
  count_pnts_bnd <- nrow(pnts_bnd)
  #
  # get sorted
  pnts_bnd_ascending <- pnts_bnd[order(pnts_bnd[, "IDB"]), ]
  #
  # initial
  ID_break <- 1
  #
  # get
  for (i in 2:count_pnts_bnd) {
    #
    # get IDB_step
    IDB_step <- pnts_bnd_ascending[i, "IDB"] - pnts_bnd_ascending[i-1, "IDB"]
    #
    # get ID
    if (IDB_step > 1) {
      #
      # assign
      ID_break <- i
      #
      # break
      break
    }
  }
  #
  # check
  if (ID_break != 1) {
    #
    # sort
    pnts_bnd_ascending <- rbind(pnts_bnd_ascending[ID_break:count_pnts_bnd, ], pnts_bnd_ascending[1:ID_break-1, ])
  }
  #
  # get pnts_curve
  pnts_curve_bnd <- pnts_bnd_ascending[, c("Cx", "Cy")]
  #
  # for initial landslide pnts, you must create a closed circle
  pnts_curve_bnd <- rbind(pnts_curve_bnd, pnts_curve_bnd[1, ])
  #
  # get intersects
  ints <- f_intersection_curves(pnts_curve_bnd, pnts_curve_segment)
  #
  # debugging
  if (is.null(ints)) {
    # #
    # # debugging
    # f_pnts_save_points(pnts, NA, "pnts")
    # #
    # # warning
    # warning("Warning: no intscts gotten, in defining anchors in an end group.")
    #
    # return
    return(NULL)
  }
  #
  # initial
  ints_dists <- c()
  #
  for (i in 1:nrow(ints)) {
    #
    # get dist
    int_dist <- 0
    int_dist <- (ints[i, 1] - lpnt0[1])^2 + int_dist
    int_dist <- (ints[i, 2] - lpnt0[2])^2 + int_dist
    int_dist <- sqrt(int_dist)
    #
    # update
    ints_dists <- c(ints_dists, int_dist)
  }
  #
  # get anchor
  anchorX <- ints[which.max(ints_dists), 1]
  anchorY <- ints[which.max(ints_dists), 2]
  #
  # initial
  distance_min <- Inf
  distance_min_ID <- c()
  #
  # get nearest bnd pnt
  for (i in 1:count_pnts_bnd) {
    #
    # get distance
    distance_i <- sqrt((anchorX - pnts_bnd[i, "Cx"])^2 + (anchorY - pnts_bnd[i, "Cy"])^2)
    #
    # check
    if (distance_i < distance_min) {
      #
      # update
      distance_min <- distance_i
      distance_min_ID <- i
    }
  }
  #
  # get anchor z, just use the elevation of bnd pnt
  anchorZ <- pnts_bnd[distance_min_ID, "Cz"]
  #
  # get anchor
  anchor <- data.frame(cbind(anchorX, anchorY, anchorZ))
  colnames(anchor) <- c("Cx", "Cy", "Cz")
  #
  # return
  return(anchor)
}
#
#
#
f_pnts_split_equal <- function(pnts, lpnt0) {
  #
  # remove lpnt0 frome pnts, might be redundant
  pnts_plnt0 <- abs(pnts[, "Cx"] - as.numeric(lpnt0[1])) + abs(pnts[, "Cy"] - as.numeric(lpnt0[2]))
  ID_remove <- which(data.frame(pnts_plnt0) == 0)
  if (length(ID_remove) >= 1) { pnts <- pnts[-ID_remove, ] }
  #
  # initial
  pnts_angles <- c()
  #
  # get angle
  for (i in 1:nrow(pnts)) {
    #
    # get xy
    lx <- as.numeric(pnts[i, "Cx"] - lpnt0[1])
    ly <- as.numeric(pnts[i, "Cy"] - lpnt0[2])
    #
    # if the same
    if (lx != 0 || ly != 0) {
      #
      # get angle
      pnts_angles <- rbind(pnts_angles, atan2(ly, lx))
    }
  }
  #
  # get circular
  pnts_angles_circular <- circular::circular(pnts_angles)
  #
  # get median
  pnts_angles_circular_median <- circular::median.circular(pnts_angles_circular)
  #
  # get attribute, not used
  pnts_angles_circular_median_attr <- attr(pnts_angles_circular_median, "medians")
  #
  # initial
  pnts_angles_sign <- pnts_angles - pnts_angles_circular_median
  #
  # update
  pnts_angles_sign[pnts_angles_sign > pi] <- pnts_angles_sign[pnts_angles_sign > pi] - pi*2
  pnts_angles_sign[pnts_angles_sign < -pi] <- pnts_angles_sign[pnts_angles_sign < -pi] + pi*2
  #
  # get median again, using ordinary median, not circular median
  pnts_angles_sign_median <- stats::median(pnts_angles_sign)
  #
  # do not need to check if the median value is right
  # if the median value appears many time
  # the number of left and right values will be not equal
  # for example, if number of left and right values differ 1
  # the number of median value could be as large as possible
  # like, c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2)
  # the median value is just chose from the 15 zeros
  #
  # return
  return(pnts_angles_sign_median + pnts_angles_circular_median)
}
#
#
#
f_intersection_curves <- function(pnts1, pnts2) {
  #
  # get x and y
  x1 <- pnts1[, 1]
  y1 <- pnts1[, 2]
  # #
  # # do not add the first point to the end to form a circle
  # x1 <- c(x1, pnts1[1, 1])
  # y1 <- c(y1, pnts1[1, 2])
  #
  # get x and y
  x2 <- pnts2[, 1]
  y2 <- pnts2[, 2]
  #
  # convert to a sp object (spatial lines)
  l1 <-  sp::Line(matrix(c(x1, y1), nc = 2, byrow = F))
  l2 <-  sp::Line(matrix(c(x2, y2), nc = 2, byrow = F))
  ll1 <- sp::Lines(list(l1), ID = "1")
  ll2 <- sp::Lines(list(l2), ID = "1")
  sl1 <- sp::SpatialLines(list(ll1))
  sl2 <- sp::SpatialLines(list(ll2))
  #
  # debugging
  if (F) {
    #
    # save
    f_pnts_save_points(pnts1, NA, "pnts1")
    f_pnts_save_points(pnts2, NA, "pnts2")
  }
  #
  # Calculate locations where spatial lines intersect
  int.pts <- rgeos::gIntersection(sl1, sl2, byid = TRUE)
  #
  # debugging
  # if (is.na(int.pts) || is.nan(int.pts) || is.null(int.pts) || pracma::isempty(int.pts)) {
  if (is.null(int.pts) || pracma::isempty(int.pts)) {
    # #
    # # warning
    # warning("Warning: invalid intersection.")
    #
    # return
    return(NULL)
  }
  #
  # initial
  int.coords <- 0
  #
  # get class
  int.pts.class <- class(int.pts)
  int.pts.class <- int.pts.class[1]
  #
  # if point
  if (int.pts.class == "SpatialPoints") {
    #
    # get coords
    int.coords <- int.pts@coords
  }
  #
  # if line
  else if (int.pts.class == "SpatialLines") {
    #
    # get coords of intersection line
    int.coords.line <- int.pts@lines[[1]]@Lines[[1]]@coords
    #
    # get average
    int.coords.line.x <- (int.coords.line[1, 1] + int.coords.line[2, 1])/2
    int.coords.line.y <- (int.coords.line[1, 2] + int.coords.line[2, 2])/2
    #
    # get coords
    int.coords <- data.frame(t(c(int.coords.line.x, int.coords.line.y)))
  }
  #
  # if collections
  else if (int.pts.class == "SpatialCollections") {
    #
    # get coords of intersection point
    int.coords <- int.pts@pointobj@coords
    #
    # get coords of intersection line
    int.coords.line <- int.pts@lineobj@lines[[1]]@Lines[[1]]@coords
    #
    # get average
    int.coords.line.x <- (int.coords.line[1, 1] + int.coords.line[2, 1])/2
    int.coords.line.y <- (int.coords.line[1, 2] + int.coords.line[2, 2])/2
    #
    # get coords
    int.coords <- rbind(int.coords, data.frame(t(c(int.coords.line.x, int.coords.line.y))))
  }
  #
  # other type?
  else {
    #
    # return
    return(NA)
  }
  #
  # debugging
  if (F) {
    #
    # Plot line data and points of intersection
    plot(x1, y1, type = "l")
    lines(x2, y2, type = "l", col = "red")
    points(int.coords[, 1], int.coords[, 2], pch = 20, col = "blue")
  }
  #
  # return
  return(int.coords)
}
#
#
#
f_dist_pnt2line <- function(pnt, paras) {
  #
  # get xy
  x <- pnt[1, "Cx"]
  y <- pnt[1, "Cy"]
  #
  # get paras
  a <- paras[1]
  b <- paras[2]
  c <- paras[3]
  #
  # get dist
  dist <- abs(a * x + b * y + c)
  dist <- dist / sqrt(a^2 + b^2)
  #
  # return
  return(dist)
}
#
#
#
f_pnts_zip <- function(pnts, CoordX, CoordY) {
  #
  # get lm
  pnts_lm <- lm(Cz ~ Cx + Cy, data = data.frame(pnts))
  pnts_lm_coeffs <- pnts_lm[[1]]
  #
  # get center z (z in plane)
  CoordZ <-  CoordX * pnts_lm_coeffs[2] + CoordY * pnts_lm_coeffs[3] + pnts_lm_coeffs[1]
  #
  # debugging
  if (is.na(CoordZ) || is.nan(CoordZ) || is.null(CoordZ) || pracma::isempty(CoordZ)) {
    #
    # return
    retur(NA)
  }
  #
  # return
  return(CoordZ)
}
#
#
#
f_pnts_distances <- function(column_centers) {
  #
  # get count
  count_center <- nrow(column_centers)
  #
  # if only one anchor, return Inf
  if (count_center == 1) {
    #
    # return
    return(Inf)
  }
  #
  # initial
  column_distances <- c() 
  #
  # get distances
  for (i in 1:(count_center-1)) {
    #
    # get distance
    distance <- 0
    distance <- distance + (column_centers[i+1, 1] - column_centers[i, 1])^2
    distance <- distance + (column_centers[i+1, 2] - column_centers[i, 2])^2
    distance <- distance + (column_centers[i+1, 3] - column_centers[i, 3])^2
    distance <- sqrt(distance)
    #
    # update distances
    column_distances <- rbind(column_distances, distance)
  }
  # 
  # return
  return(column_distances)
}
#
#
#
f_pnts_end_aspect_ratio <- function(end_area, end_length) {
  #
  # get aspect ratio
  end_aspect_ratio <- end_length / (end_area / end_length)
  #
  # return
  return (end_aspect_ratio)
}
#
#
#
f_pnts_anchors_deviation_angle <- function(column_centers, ends) {
  #
  # initial vectors
  vector1 <- c()
  vector2 <- c()
  # 
  # initial end
  if (ends == 1) {
    #
    # get vectors
    vector1 <- c(column_centers[3, 1] - column_centers[2, 1], column_centers[3, 2] - column_centers[2, 2])
    vector2 <- c(column_centers[4, 1] - column_centers[3, 1], column_centers[4, 2] - column_centers[3, 2])
  }
  #
  # distal end
  if (ends == 2) {
    #
    # get count
    column_centers_nrow <- nrow(column_centers)
    #
    # get vectors
    vector1 <- c(column_centers[column_centers_nrow-2, 1] - column_centers[column_centers_nrow-3, 1], column_centers[column_centers_nrow-2, 2] - column_centers[column_centers_nrow-3, 2])
    vector2 <- c(column_centers[column_centers_nrow-1, 1] - column_centers[column_centers_nrow-2, 1], column_centers[column_centers_nrow-1, 2] - column_centers[column_centers_nrow-2, 2])
  }
  #
  # get deviation angle
  vectors_dvt_agl <- f_get_vectors2d_angle(vector1, vector2)
  #
  # return
  return(vectors_dvt_agl)
}
#
#
#
f_get_vectors2d_angle <- function(vector1, vector2) {
  #
  # get deviation angle cosine
  vectors_dot <- crossprod(vector1, vector2)
  vectors_cos <- vectors_dot / sqrt(vector1[1]^2 + vector1[2]^2) / sqrt(vector2[1]^2 + vector2[2]^2)
  #
  # numerical error
  if (abs(vectors_cos) > 1) {
    #
    # adjust
    vectors_cos <- vectors_cos / abs(vectors_cos)
  }
  #
  # get deviation angle
  vectors_agl <- acos(vectors_cos) * 180 / pi
  #
  # for debugging
  if (is.nan(vectors_agl)) {
    #
    # warning
    warning("Warning: vectors angle is NaN.")
  }
  #
  # return
  return (vectors_agl)
}
#
#
#
f_list2pnts <- function(list_pnts) {
  #
  # initial
  pnts <- c()
  #
  # convert list to 
  for (i in 1:length(list_pnts)) {
    #
    # get pnts_i
    pnts_i <- list_pnts[[i]]
    #
    # rbind
    pnts <- rbind(pnts, pnts_i)
  }
  # 
  # return
  return(pnts)
}
#
#
#
f_pnts_save_xls <- function(pnts, file_out) {
  #
  # write
  xlsx::write.xlsx(pnts, file = file_out, sheetName = "pnts", row.names = FALSE, col.names = TRUE)
}
#
#
#
f_pnts_save_points <- function(pnts, CRS_out, file_out) {
  #
  # get SpatialPointsDataFrame
  PntsCoords_sptdf <- f_pnts2shp_points(pnts, CRS_out)
  #
  # write
  rgdal::writeOGR(PntsCoords_sptdf, dsn = ".", layer = file_out, overwrite = T, driver="ESRI Shapefile")
  # #
  # # write 3D, will produce extra Cx, Cy and Cz columns in shapefile???
  # # even, "Warning message: writePointsShape is deprecated; use rgdal::writeOGR or sf::st_write"
  # rgdal::writePointsShape(PntsCoords_sptdf, file_out, factor2char = TRUE, max_nchar = 254)
  # #
  # # indeed, writePointsShape() does not write the projection file
  # # but, using function showWKT from rgdal, you can also create one like that
  # cat(showWKT(sp::proj4string(PntsCoords_sptdf)), file = paste(file_out, ".prj", sep = "", collapse = NULL))
  #
  # return
  return(PntsCoords_sptdf)
}
#
#
#
f_pnts_save_proflle <- function(pnts_prfl, spxdf_lasld, file_out, LineID = 1, paras_input = c(3, 0, 0, 180)) {
  #
  # # solve "Error in plot.new() : figure margins too large"
  # par(mar = c(1, 1, 1, 1))
  # #
  # # plot
  # plot(spxdf_lasld, main = "lasld DEM")
  # # 
  # # write
  # rgdal::writeGDAL(spxdf_lasld, "spxdf_lasld.tiff")
  #
  # get areas
  area3d <- sp::surfaceArea(spxdf_lasld)
  area2d <- length(spxdf_lasld@grid.index) * spxdf_lasld@grid@cellsize[1] * spxdf_lasld@grid@cellsize[2]
  #
  # get parameters for landslide profile
  paras <- f_pnts2paras_profile(pnts_prfl, c(area3d, area2d))
  #
  # get data.frame
  paras_input <- data.frame(iMGPC = paras_input[1], iMGAD = paras_input[2], iMEAR = paras_input[3], iMEDA = paras_input[4])
  #
  # get data_df, inputs leed paras
  data_df <- cbind(paras_input, data.frame(ID = LineID), paras)
  #
  # get spldf
  spldf <- f_pnts2shp_polyline(pnts_prfl, spxdf_lasld@proj4string, data_df)
  #
  # write
  rgdal::writeOGR(spldf, dsn = ".", layer = file_out, overwrite = T, driver="ESRI Shapefile")
  #
  # return
  return(spldf)
}
#
#
#
f_pnts2shp_points <- function(pnts, CRS_out) {
  #
  # get coords
  PntsCoords <- pnts
  # sp::coordinates(PntsCoords) = ~ Cx + Cy + Cz
  sp::coordinates(PntsCoords) = ~ Cx + Cy
  #
  # get SpatialPointsDataFrame
  PntsCoords_sptdf <- sp::SpatialPointsDataFrame(PntsCoords, pnts)
  #
  # get prj
  sp::proj4string(PntsCoords_sptdf) = CRS_out
  # #
  # # write, in another function
  # rgdal::writeOGR(PntsCoords_sptdf, dsn = ".", layer = "PntsCoords_sptdf", overwrite = T, driver="ESRI Shapefile")
  #
  # return
  return(PntsCoords_sptdf)
}
#
#
#
f_pnts2shp_polyline <- function(pnts, CRS_out, data_df = NULL) {
  #
  # get data_spl
  # no R library supports 3d line now?
  # data_spl <- cbind(pnts[, "Cx"], pnts[, "Cy"], pnts[, "Cz"])
  data_spl <- cbind(pnts[, "Cx"], pnts[, "Cy"])
  # rownames(data_spl) <- c(1:nrow(data_spl))
  #
  # check
  if (is.null(data_df)) {
    #
    # get data_df
    data_df <- data.frame(c(1))
    #
    # rownames
    rownames(data_df) <- c("a")
  }
  #
  # get spldf
  spldf <- f_data2spldf(data_spl, data_df)
  #
  # get prj
  sp::proj4string(spldf) = CRS_out
  # #
  # # write
  # rgdal::writeOGR(spldf, dsn = ".", layer = "spldf", overwrite = T, driver="ESRI Shapefile")
  #
  # return
  return(spldf)
}
#
#
#
f_pnts2shp_polygon <- function(pnts, CRS_out, data_df = NULL) {
  #
  # create a Polygon
  p <- sp::Polygon(pnts)
  #
  # wrap that into a Polygons object
  ps <- sp::Polygons(list(p), 1)
  #
  # wrap that into a SpatialPolygons object
  sps <- sp::SpatialPolygons(list(ps))
  #
  # check
  if (is.null(data_df)) {
    #
    # get data_df
    data_df <- data.frame(c(1))
    #
    # rownames
    rownames(data_df) <- c("1")
  }
  #
  # get spgdf
  spgdf <- sp::SpatialPolygonsDataFrame(sps, data_df)
  #
  # get prj
  sp::proj4string(spgdf) = CRS_out
  # #
  # # write
  # rgdal::writeOGR(spgdf, dsn = ".", layer = "spgdf", overwrite = T, driver="ESRI Shapefile")
  #
  # return
  return(spgdf)
}
#
#
#
f_pnts2shp_polygon_rasterized <- function(pnts, raster_cellsize, raster_CRS) {
  #
  # initialize
  pnts_raster0 <- raster::raster()
  #
  # get extent
  x_min <- min(pnts[, "Cx"]) - raster_cellsize/2
  x_max <- max(pnts[, "Cx"]) + raster_cellsize/2
  #
  y_min <- min(pnts[, "Cy"]) - raster_cellsize/2
  y_max <- max(pnts[, "Cy"]) + raster_cellsize/2
  #
  raster_extent <- raster::extent(x_min, x_max, y_min, y_max)
  #
  # get extent
  raster::extent(pnts_raster0) <- raster::extent(raster_extent)
  #
  # get cnts
  cnt_col <- (x_max - x_min) / raster_cellsize
  cnt_row <- (y_max - y_min) / raster_cellsize
  #
  # get cnts
  raster::ncol(pnts_raster0) <- cnt_col
  raster::nrow(pnts_raster0) <- cnt_row
  #
  # sum of the values associated with the points
  pnts_raster <- raster::rasterize(pnts[, c("Cx", "Cy")], pnts_raster0, pnts[, c("Cz")], fun = sum)
  #
  # assign
  raster::crs(pnts_raster) <- raster_CRS
  #
  # set one
  pnts_raster_set1 <- raster::calc(pnts_raster, function(x) x / x)
  #
  # raster to polygon
  pnts_raster_set1_pgn <- raster::rasterToPolygons(pnts_raster_set1, fun = NULL, na.rm = TRUE, dissolve = TRUE)
  #
  # debugging
  if (F) {
    #
    # get
    pnts_raster_spxdf <- as(pnts_raster, "SpatialPixelsDataFrame")
    #
    # write raster
    rgdal::writeGDAL(pnts_raster_spxdf, "pnts_raster.tiff")
    #
    # plot raster
    plot(pnts_raster_set1)
    #
    # plot shapefile
    # notice that you use add = T to add a layer on top of an existing plot in R.
    plot(pnts_raster_set1_pgn, main = "polygon set1", axes = TRUE, border = "blue", add = TRUE)
    #
    # write shapefile
    rgdal::writeOGR(pnts_raster_set1_pgn, ".", "pnts_raster_set1_pgn", overwrite = T, driver = "ESRI Shapefile")
  }
  #
  # return
  return(pnts_raster_set1_pgn)
}
#
#
#
f_pnts2paras_profile <- function(pnts, areas) {
  #
  # get count
  count_pnts = nrow(pnts)
  #
  # initial, displacement for the whole profile
  l_profile = 0
  h_profile = 0
  s_profile = 0
  #
  # initial, points in the profile
  l_series = 0
  h_series = 0
  #
  # get profile
  for (i in 2:count_pnts) {
    #
    x1 = pnts[i-1, 1]
    y1 = pnts[i-1, 2]
    z1 = pnts[i-1, 3]
    #
    x2 = pnts[i, 1]
    y2 = pnts[i, 2]
    z2 = pnts[i, 3]
    #
    l = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    h = z2 - z1
    s = sqrt(l^2 + h^2)
    #
    l_profile = l_profile + l
    h_profile = h_profile + h
    s_profile = s_profile + s
    #
    l_series = rbind(l_series, l_series[i-1] + l)
    h_series = rbind(h_series, h_series[i-1] + h)
  }
  #
  # get apparent miu
  miuapp_profile = -h_profile / l_profile
  #
  # initial
  g_series <- 0
  #
  # get residual height
  for (i in 2:count_pnts) {
    #
    # get residual height
    g_series <- rbind(g_series, (-h_series[i]) - miuapp_profile * l_series[i])
  }
  #
  # initial
  H_integral <- 0
  G_integral <- 0
  #
  # get integral
  for (i in 2:count_pnts) {
    #
    # get l
    l <- l_series[i] - l_series[i-1]
    #
    # get area
    Area_H <- l * (-h_series[i] + -h_series[i-1]) / 2
    Area_G <- l * (g_series[i] + g_series[i-1]) / 2
    #
    # get int
    H_integral <- H_integral + Area_H
    G_integral <- G_integral + Area_G
  }
  #
  # get tortuosity overall
  d_overall <- sqrt((pnts[count_pnts, 1] - pnts[1, 1])^2 + 
                    (pnts[count_pnts, 2] - pnts[1, 2])^2 +
                    (pnts[count_pnts, 3] - pnts[1, 3])^2)
  T_overall <- s_profile / d_overall
  #
  # get tortuosity horizontal
  d_horizontal <- sqrt((pnts[count_pnts, 1] - pnts[1, 1])^2 + 
                       (pnts[count_pnts, 2] - pnts[1, 2])^2)
  T_horizontal <- l_profile / d_horizontal
  #
  # get tortuosity longitudinal
  d_longitudinal <- sqrt((l_profile)^2 +
                         (h_profile)^2)
  T_longitudinal <- s_profile / d_longitudinal
  #
  # get paras
  paras <- data.frame(Dall = d_overall, Lall = s_profile,
                      Dhrz = d_horizontal, Lhrz = l_profile,
                      Dlng = d_longitudinal,
                      Hfnl = -h_profile,
                      Hmax = -min(h_series), Hint = H_integral,
                      Gmax = max(g_series), Gint = G_integral,
                      Tall = T_overall, Thrz = T_horizontal, Tlng = T_longitudinal,
                      miuapp = miuapp_profile, 
                      Aall = areas[1], Wall = areas[1] / s_profile, epsall = s_profile^2 / areas[1],
                      Ahrz = areas[2], Whrz = areas[2] / l_profile, epshrz = l_profile^2 / areas[2])
  #
  # return
  return(paras)
}
#
#
#
f_data2spldf <- function(data_spl, data_df) {
  #
  # get spl
  pLine <- sp::Line(data_spl)
  pLines <- sp::Lines(list(pLine), ID = "a")
  pSpatialLines <- sp::SpatialLines(list(pLines))
  #
  # get df
  rownames(data_df) <- c("a")
  #
  # get
  spldf <- sp::SpatialLinesDataFrame(pSpatialLines, data = data_df)
  # #
  # plot
  # plot(spldf, col = c("red"))
  #
  # return
  return(spldf)
}
#
#
#
f_save_list2spldf <- function(list_spldf, CRS_out, file_out) {
  #
  # get count
  count_spldf <- length(list_spldf)
  #
  # initial
  list_SpatialLines <- list()
  #
  # initial
  data_df <- data.frame()
  #
  # get stack
  for (i in 1:count_spldf) {
    #
    # get spldf
    spldf_i <- list_spldf[[i]]
    #
    # get Lines
    pLine <- sp::Line(spldf_i@lines[[1]]@Lines[[1]]@coords)
    pLines <- sp::Lines(list(pLine), ID = i)
    #
    # update list
    list_SpatialLines[[i]] <- pLines
    #
    # get data
    data_df <- rbind(data_df, spldf_i@data)
  }
  #
  # get SpatialLines
  pSpatialLines = sp::SpatialLines(list_SpatialLines)
  #
  # set names
  rownames(data_df) <- c(1:nrow(data_df))
  #
  # get spldf
  spldf <-  sp::SpatialLinesDataFrame(pSpatialLines, data = data_df)
  # plot(spldf, col = c("red"))
  #
  # get prj
  sp::proj4string(spldf) = CRS_out
  #
  # write
  rgdal::writeOGR(spldf, dsn = ".", layer = file_out, overwrite = T, driver="ESRI Shapefile")
  #
  # return
  return(data_df)
}
#
#
#
f_show_step <- function(logging) {
  #
  # show progress
  cat(logging)
  #
  # get time start
  time_step_start <- Sys.time()
  #
  # return
  return(time_step_start)
}
#
#
#
f_show_time <- function(logging, time_step_start) {
  #
  # show progress
  cat(logging)
  #
  # get time end
  time_step_end <- Sys.time()
  #
  # show time
  cat(paste(as.character(as.numeric(time_step_end - time_step_start,  units = "secs")), " seconds used.", "\n", sep = "", collapse = NULL))
}