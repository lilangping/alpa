#
# ###################################################################### #
# ###################################################################### #
#                                                                        #
#                                ALPA                                    #
#                                                                        #
#      An R-Script for for automatic analysis of landslide profile       #
#                                                                        #
#                             Version 2.1                                #
#                                                                        #
#                            May 5th, 2021                               #
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
ALPA <- function(FileDEM, FileLandslides, MinGrpPntCnt = 3, MinGrpAcrDst = 0, MinEndAptRto = 0, MaxEndDvtAgl = 180, MinStpHrzLen = 30, FilePrefix = "", OutputTemp = F) {
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
    # get polygons
    pPolygons <- list_polygons[[i]]
    #
    # split using one polygons
    spldf_profile <- f_lasld_split(raster_dem, pPolygons, MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, MinStpHrzLen, FilePrefix, i, OutputTemp)
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
f_lasld_split <- function(pRasterDEM, pPolygons, MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, MinStpHrzLen, FilePrefix, FileID, OutputTemp) {
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
    # get path (centers and bnd sides)
    list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons)
    # get centers
    pnts_split_anchors <- list_pnts_path[[1]]
    #
    # if save
    if (OutputTemp || (min(data.frame(list_grpm)) == 1)) {
      #
      # get strips
      list_strips <- f_strips(pRasterDEM, pPolygons, list_pnts_path, MinStpHrzLen)
      # get 
      strip_sp <- list_strips[[1]]
      strip_df <- list_strips[[2]]
      # get 
      strip_nodes <- list_strips[[3]]
      #
      # get paras
      paras_strip <- c(sum(strip_df[, "Lhrz"]), sum(strip_df[, "Ahrz"]), sum(strip_df[, "Lall"]), sum(strip_df[, "Aall"]))
    }
    #
    # save for every split
    if (OutputTemp) {
      #
      # save pnts
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), ".xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split, raster_CRS, file_out)
      #
      # save anchors
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_anchors.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split_anchors, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_anchors", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split_anchors, raster_CRS, file_out)
      #
      # save strips
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_df, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips", sep = "", collapse = NULL)
      f_save_list2sp(strip_sp, strip_df, raster_CRS, file_out)
      #
      # save nodes
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_nodes.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_nodes, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_nodes", sep = "", collapse = NULL)
      f_pnts_save_points(strip_nodes, raster_CRS, file_out)
      #
      # save shp pln
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_nodes_pln2d", sep = "", collapse = NULL)
      spldf_profile <- f_pnts_save_proflle(strip_nodes, paras_strip, c(MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, MinStpHrzLen), raster_CRS, file_out, FileID)
      #
      # save paras
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_nodes_pln2d.xlsx", sep = "", collapse = NULL)
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
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split, raster_CRS, file_out)
      #
      # save anchors
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps_anchors.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split_anchors, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps_anchors", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split_anchors, raster_CRS, file_out)
      #
      # save strips
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps_strips.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_df, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps_strips", sep = "", collapse = NULL)
      f_save_list2sp(strip_sp, strip_df, raster_CRS, file_out)
      #
      # save nodes
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps_strips_nodes.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_nodes, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps_strips_nodes", sep = "", collapse = NULL)
      f_pnts_save_points(strip_nodes, raster_CRS, file_out)
      #
      # save shp pln
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps_strips_nodes_pln2d", sep = "", collapse = NULL)
      spldf_profile <- f_pnts_save_proflle(strip_nodes, paras_strip, c(MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, MinStpHrzLen), raster_CRS, file_out, FileID)
      #
      # save paras
      file_out <- paste(FilePrefix, "_l", sprintf("%05d", FileID), "_pnts_grps_strips_nodes_pln2d.xlsx", sep = "", collapse = NULL)
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
        #
        # set bool end
        bool_initial <- pnts_i[1, "grp"] == 1
        bool_distal <- pnts_i[1, "grp"] == grp_max
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
        # #
        # # show step
        # time_step_start <- f_show_step(paste("Check spatial impartiality start.", "\n", sep = ""))
        #
        # get sides
        pnts_i_upper_bnd <- pnts_i_upper[which(pnts_i_upper[, "IDB"] != 0), ]
        list_bnd_sides_upper <- f_pnts_bnd_sides(pnts_i_upper_bnd, IDB_max)
        #
        # if sub- pnts do not have boundary points on two sides
        if ((bool_initial && length(list_bnd_sides_upper) != 1) ||
            (!bool_initial && length(list_bnd_sides_upper) != 2)) {
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
        # get sides
        pnts_i_lower_bnd <- pnts_i_lower[which(pnts_i_lower[, "IDB"] != 0), ]
        list_bnd_sides_lower <- f_pnts_bnd_sides(pnts_i_lower_bnd, IDB_max)
        #
        # if sub- pnts do not have boundary points on two sides
        if ((bool_distal && length(list_bnd_sides_lower) != 1) ||
            (!bool_distal && length(list_bnd_sides_lower) != 2)) {
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
        # #
        # # show step time
        # f_show_time(paste("Check spatial impartiality done.", "\n", sep = ""), time_step_start)
        #
        # check ability to prolong for end group, initial
        if (bool_initial) {
          #
          # get launch anchor for end groups
          pnts_end_LaunchAnchor <- f_pnts_end_LaunchAnchor(pnts_i_upper)
          #
          # get bnd sides
          list_pnts_i_upper_bnd_sides <- f_pnts_bnd_sides(pnts_i_upper, IDB_max, pnts_end_LaunchAnchor)
          #
          # check
          if (is.na(list_pnts_i_upper_bnd_sides)) {
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
        }
        #
        # check ability to prolong for end group, distal
        if (bool_distal) {
          #
          # get launch anchor for end groups
          pnts_end_LaunchAnchor <- f_pnts_end_LaunchAnchor(pnts_i_lower)
          #
          # get bnd sides
          list_pnts_i_lower_bnd_sides <- f_pnts_bnd_sides(pnts_i_lower, IDB_max, pnts_end_LaunchAnchor)
          #
          # check
          if (is.na(list_pnts_i_lower_bnd_sides)) {
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
        }
        #
        # get anchors and distances (before split)
        pnts_split_anchors_before_split <- pnts_split_anchors
        # pnts_split_distances_before_split <- f_pnts_distances(pnts_split_anchors[, c("Cx", "Cy", "Cz")])
        #
        # make a copy
        list_pnts_new_for_distances <- list_pnts_new
        #
        # update list_pnts_new
        list_pnts_new_for_distances[[grpu]] <- pnts_i_upper
        list_pnts_new_for_distances[[grpu + 1]] <- pnts_i_lower
        #
        # if still have pnts, append
        if ((i+1) <= length(list_pnts)) {
          #
          # append
          for (k in (i+1):length(list_pnts)) {
            #
            # get pnts_i
            pnts_k <- list_pnts[[k]]
            list_pnts_new_for_distances[[length(list_pnts_new_for_distances) + 1]] <- pnts_k
          }
        }
        #
        # get path (centers and bnd sides)
        list_pnts_path <- f_pnts_path(list_pnts_new_for_distances, IDB_max, pPolygons)
        # get centers
        pnts_split_anchors <- list_pnts_path[[1]]
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
          # sometimes, group center anchor and inter-group center anchor will be the same?
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
        # check deviation angle for end groups
        if (length(list_pnts) > 1) {
          #
          # do not check deviation angle for the very initial landslide pnts
          # the very initial landslide pnts is both initial end and distal end
          # sometimes, length of the profile of the very initial landslide pnts
          # is very short because path crossing bnd is not allowed
          # this will cause very small aspect ratio
          #
          # and, do not use the generated length of end groups
          # use the possible most prolong length along the moving direction
          # the length of generated segment could be short when path is narrow
          #
          # handle end initial
          if (bool_initial) {
            #
            # get area and length of the split pnts
            end_area <- nrow(pnts_i) * raster_cellsize^2
            end_length <- f_pnts_end_prolong(pnts_i, pnts_split_anchors_before_split[c(3, 1), ])
            #
            # get aspect ratio for end group after split
            EndAptRto_initial <- end_length^2 / end_area
            #
            # if adequate small
            if (EndAptRto_initial < MinEndAptRto) {
              #
              # get deviation angle for end group after split
              EndDvtAgl_initial <- f_get_vectors2d_deviation_angle(pnts_split_anchors[c(4,3,2), ])
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
          if (bool_distal) {
            #
            # get area and length of the split pnts
            end_area <- nrow(pnts_i) * raster_cellsize^2
            end_length <- f_pnts_end_prolong(pnts_i, pnts_split_anchors_before_split[c(nrow(pnts_split_anchors_before_split)-2, nrow(pnts_split_anchors_before_split)), ])
            #
            # get aspect ratio for end group after split
            EndAptRto_distal <- end_length^2 / end_area
            #
            # if adequate small
            if (EndAptRto_distal < MinEndAptRto) {
              #
              # get deviation angle for end group after split
              EndDvtAgl_distal <- f_get_vectors2d_deviation_angle(pnts_split_anchors[(nrow(pnts_split_anchors)-3):(nrow(pnts_split_anchors)-1), ])
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
  # shp_lasld is organized clockwise
  # change from clockwise to counterclockwise
  IDsB <- which(pnts[, "IDB"] != 0)
  pnts[IDsB, "IDB"] <- max(pnts[IDsB, "IDB"]) + 1 - pnts[IDsB, "IDB"]
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
f_pnts_end_prolong <- function(pnts, vectors) {
  #
  # get vector1
  vector1_x <- vectors[1, "Cx"]
  vector1_y <- vectors[1, "Cy"]
  #
  # get vector, normalized
  vector_x <- vectors[2, "Cx"] - vector1_x
  vector_y <- vectors[2, "Cy"] - vector1_y
  vector_length <- sqrt(vector_x^2 + vector_y^2)
  vector_x <- vector_x / vector_length
  vector_y <- vector_y / vector_length
  #
  # get count
  count_pnt <- nrow(pnts)
  #
  # initial
  column_length <- c()
  #
  # get prolong length
  for (i in 1:count_pnt) {
    #
    # get vector
    vector_i_x <- pnts[i, "Cx"] - vector1_x
    vector_i_y <- pnts[i, "Cy"] - vector1_y
    #
    # append
    column_length <- c(column_length, vector_i_x * vector_x + vector_i_y * vector_y)
  }
  #
  # return
  return(max(column_length))
}
#
#
#
f_get_vectors2d_deviation_angle <- function(pnts) {
  #
  # get vectors
  vector1 <- c(pnts[2, "Cx"] - pnts[1, "Cx"], pnts[2, "Cy"] - pnts[1, "Cy"])
  vector2 <- c(pnts[3, "Cx"] - pnts[2, "Cx"], pnts[3, "Cy"] - pnts[2, "Cy"])
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
f_pnts_end_LaunchAnchor <- function(pnts) {
  #
  # get pnts_bnd
  pnts_bnd <- pnts[which(pnts[, "IDB"] != 0), ]
  #
  # get count
  count_pnt_bnd <- nrow(pnts_bnd)
  #
  # no bnd
  if (count_pnt_bnd == 0) { return(NULL) }
  #
  # only one bnd
  if (count_pnt_bnd == 1) {
    #
    # return
    return(pnts_bnd[1, c("Cx", "Cy")])
  }
  #
  # get sorted
  pnts_bnd_ascending <- f_pnts_bnd_ascending(pnts_bnd)
  #
  # get
  LaunchAnchorX <- (pnts_bnd_ascending[1, "Cx"] + pnts_bnd_ascending[count_pnt_bnd, "Cx"]) / 2
  LaunchAnchorY <- (pnts_bnd_ascending[1, "Cy"] + pnts_bnd_ascending[count_pnt_bnd, "Cy"]) / 2
  #
  # return
  return(c(LaunchAnchorX, LaunchAnchorY))
}
#
#
#
f_pnts_end_LaunchAnchor_initial <- function(list_bnd_sides, IDB_max) {
  #
  # get bnd sides
  bnd_sides_right <- list_bnd_sides[[1]]
  bnd_sides_left <- list_bnd_sides[[2]]
  #
  # get IDB start
  IDB_start_right <- bnd_sides_right[1, "IDB"]
  IDB_start_left <- bnd_sides_left[1, "IDB"]
  #
  # get IDB end
  IDB_end_right <- bnd_sides_right[nrow(bnd_sides_right), "IDB"]
  IDB_end_left <- bnd_sides_left[nrow(bnd_sides_left), "IDB"]
  #
  # get IDB difference
  IDB_start_right_end_left <- IDB_start_right - IDB_end_left
  IDB_start_left_end_right <- IDB_start_left - IDB_end_right
  #
  # check
  if (IDB_start_right_end_left < 0) { IDB_start_right_end_left <- IDB_start_right_end_left + IDB_max }
  if (IDB_start_left_end_right < 0) { IDB_start_left_end_right <- IDB_start_left_end_right + IDB_max }
  #
  # return
  if (IDB_start_right_end_left <= 1) { return(bnd_sides_right[1, ]) }
  if (IDB_start_left_end_right <= 1) { return(bnd_sides_left[1, ]) }
}
#
#
#
f_pnts_bnd_ascending <- function(pnts_bnd) {
  #
  # get sorted
  pnts_bnd_ascending <- pnts_bnd[order(pnts_bnd[, "IDB"]), ]
  #
  # initial
  ID_start <- 1
  #
  # get
  for (i in 2:nrow(pnts_bnd_ascending)) {
    #
    # get IDB_step
    IDB_step <- pnts_bnd_ascending[i, "IDB"] - pnts_bnd_ascending[i-1, "IDB"]
    #
    # get ID
    if (IDB_step > 1) {
      #
      # assign
      ID_start <- i
      #
      # break
      break
    }
  }
  #
  # check
  if (ID_start != 1) {
    #
    # sort
    pnts_bnd_ascending <- rbind(pnts_bnd_ascending[ID_start:nrow(pnts_bnd_ascending), ], pnts_bnd_ascending[1:ID_start-1, ])
  }
  #
  # return
  return(pnts_bnd_ascending)
}
#
#
#
f_pnts_bnd_sides <- function(pnts, IDB_max, LaunchAnchor = c()) {
  #
  # provide launch anchor for end group or not?
  if (nrow(data.frame(LaunchAnchor)) != 0) {
    #
    # for the very initial pnts
    # i.e., the overall landslide polygon is here
    # the launch pnt (anchor), i.e. the highest pnt must be removed
    # to form a gap in the IDB pnts
    #
    # get bnd sides
    list_bnd_sides <- f_pnts_split_equal_fromApnt(pnts, LaunchAnchor)
    #
    # return
    return(list_bnd_sides)
  }
  #
  # get pnts_bnd
  pnts_bnd <- pnts[which(pnts[, "IDB"] != 0), ]
  #
  # get count
  count_pnt <- nrow(pnts_bnd)
  #
  # no bnd
  if (count_pnt == 0) { return(NULL) }
  #
  # only one bnd
  if (count_pnt == 1) {
    #
    # return
    return(list(pnts_bnd))
  }
  #
  # get sorted
  pnts_bnd_ascending <- pnts_bnd[order(pnts_bnd[, "IDB"]), ]
  #
  # initial
  list_bnd_sides <- list()
  count_sides <- 0
  #
  # initial
  bnd_sides <- pnts_bnd_ascending[1, ]
  #
  # get sides
  for (i in 2:count_pnt) {
    #
    # get IDB step
    IDB_step <- pnts_bnd_ascending[i, "IDB"] - pnts_bnd_ascending[i-1, "IDB"]
    #
    # debugging
    if (is.na(IDB_step)) {
      #
      # warning
      warning("IDB_step is NA")
    }
    #
    # consecutive or not
    if (IDB_step == 1) {
      #
      # append
      bnd_sides <- rbind(bnd_sides, pnts_bnd_ascending[i, ])
    }
    else {
      #
      # update
      count_sides <- count_sides + 1
      list_bnd_sides[[count_sides]] <- bnd_sides
      #
      # update
      bnd_sides <- pnts_bnd_ascending[i, ]
    }
  }
  #
  # if consecutive
  if (length(list_bnd_sides) == 0) {
    #
    # update
    count_sides <- count_sides + 1
    list_bnd_sides[[count_sides]] <- bnd_sides
  }
  #
  # if cross IDB_max
  else if (pnts_bnd_ascending[1, "IDB"] == 1 && pnts_bnd_ascending[count_pnt, "IDB"] == IDB_max) {
    #
    # update
    bnd_sides <- rbind(bnd_sides, list_bnd_sides[[1]])
    list_bnd_sides[[1]] <- bnd_sides
  }
  #
  # if not
  else {
    #
    # update
    count_sides <- count_sides + 1
    list_bnd_sides[[count_sides]] <- bnd_sides
  }
  #
  # return
  return(list_bnd_sides)
}
#
#
#
f_pnts_bnd_sides_switch2 <- function(list_bnd_sides, LaunchAnchor, IDB_max) {
  #
  # if null
  if (is.null(list_bnd_sides)) { return(NULL) }
  #
  # get bnd sides
  bnd_sides_right <- list_bnd_sides[[1]]
  bnd_sides_left <- list_bnd_sides[[2]]
  #
  # get IDB start
  IDB_start_right <- bnd_sides_right[1, "IDB"]
  IDB_start_left <- bnd_sides_left[1, "IDB"]
  #
  # get IDB launch anchor
  IDB_launch_anchor <- LaunchAnchor[1, "IDB"]
  #
  # get IDB diff
  IDB_diff_right <- IDB_start_right - IDB_launch_anchor
  IDB_diff_left <- IDB_start_left - IDB_launch_anchor
  #
  # check
  if (IDB_diff_right < 0) { IDB_diff_right <- IDB_diff_right + IDB_max }
  if (IDB_diff_left < 0) { IDB_diff_left <- IDB_diff_left + IDB_max }
  #
  # return
  if (IDB_diff_right < IDB_diff_left) { return(list_bnd_sides) }
  if (IDB_diff_right > IDB_diff_left) { return(list(list_bnd_sides[[2]], list_bnd_sides[[1]])) }
  #
  # return, possibly all are zeros
  return(list_bnd_sides)
}
#
#
#
f_pnts_bnd_sides_corners <- function(list_bnd_sides) {
  #
  # if null
  if (is.null(list_bnd_sides)) { return(NULL) }
  #
  # initial
  bnd_sides_corners <- c()
  #
  # get count
  count_bnd_sides <- length(list_bnd_sides)
  #
  # get corners
  for (i in 1:count_bnd_sides) {
    #
    # get bnd_sides
    bnd_sides <- list_bnd_sides[[i]]
    #
    # append
    bnd_sides_corners <- rbind(bnd_sides_corners, bnd_sides[1, ], bnd_sides[nrow(bnd_sides), ])
  }
  #
  # return
  return(bnd_sides_corners)
}
#
#
#
f_pnts_bnd_sides_centerLs <- function(list_bnd_sides) {
  #
  # if null
  if (is.null(list_bnd_sides)) { return(NULL) }
  #
  # get bnd sides corners()
  bnd_sides_corners <- f_pnts_bnd_sides_corners(list_bnd_sides)
  #
  # get center xy
  centerL_x1 <- (bnd_sides_corners[4, "Cx"] + bnd_sides_corners[1, "Cx"]) / 2
  centerL_y1 <- (bnd_sides_corners[4, "Cy"] + bnd_sides_corners[1, "Cy"]) / 2
  centerL1 <- c(centerL_x1, centerL_y1)
  #
  # get center xy
  centerL_x2 <- (bnd_sides_corners[3, "Cx"] + bnd_sides_corners[2, "Cx"]) / 2
  centerL_y2 <- (bnd_sides_corners[3, "Cy"] + bnd_sides_corners[2, "Cy"]) / 2
  centerL2 <- c(centerL_x2, centerL_y2)
  #
  # get center middle
  centerL_x_mid <- (centerL_x1 + centerL_x2) / 2
  centerL_y_mid <- (centerL_y1 + centerL_y2) / 2
  centerL_mid <- c(centerL_x_mid, centerL_y_mid)
  #
  # return
  return(rbind(centerL1, centerL2, centerL_mid))
}
#
#
#
f_pnts_split_equal_fromApnt <- function(pnts, lpnt0) {
  #
  # do not use circular statistics
  # just loop all the bnd pnts, to see which bnd pnt or bnd pair splits most equally
  # also, check if the most equally bnd pnt is blocked by other bnd pnt
  # straight flow path that crosses the bnd is not accepted
  # it might be slightly slower than using circular statistics
  # but much more robust and concise
  #
  # situation is very complex
  # so, check unblocked for all bnd pnts
  # although might cost more time, but more robust and concise
  #
  # initial, remove lpnt0 from pnts, might be redundant
  pnts_lpnt0_rm <- c()
  #
  # initial
  pnts_angles <- c()
  #
  # get angles
  for (i in 1:nrow(pnts)) {
    #
    # get xy
    lx <- as.numeric(pnts[i, "Cx"] - lpnt0[1])
    ly <- as.numeric(pnts[i, "Cy"] - lpnt0[2])
    #
    # if not the same
    if (lx != 0 || ly != 0) {
      #
      # append
      pnts_lpnt0_rm <- rbind(pnts_lpnt0_rm, pnts[i, ])
      #
      # get angle
      pnts_angles <- c(pnts_angles, atan2(ly, lx))
    }
  }
  #
  # get bnd pnts
  pnts_lpnt0_rm_bnd <- pnts_lpnt0_rm[which(pnts_lpnt0_rm[, "IDB"] != 0), ]
  #
  # get sorted
  pnts_lpnt0_rm_bnd_ascending <- f_pnts_bnd_ascending(pnts_lpnt0_rm_bnd)
  #
  # get pnts_curve, for bnd pnts
  pnts_curve_bnd <- pnts_lpnt0_rm_bnd_ascending[, c("Cx", "Cy")]
  #
  # get count
  count_pnt_bnd <- nrow(pnts_curve_bnd)
  #
  # initial
  pnts_angles_bnd <- c()
  #
  # get angles, for bnd pnts
  for (i in 1:count_pnt_bnd) {
    #
    # get xy
    lx <- as.numeric(pnts_curve_bnd[i, "Cx"] - lpnt0[1])
    ly <- as.numeric(pnts_curve_bnd[i, "Cy"] - lpnt0[2])
    #
    # if not the same, theoretically will not happen here
    if (lx != 0 || ly != 0) {
      #
      # get angle
      pnts_angles_bnd <- c(pnts_angles_bnd, atan2(ly, lx))
    }
  }
  #
  # initial
  column_count_diff <- c()
  #
  # get difference
  for (i in 1:count_pnt_bnd) {
    #
    # get angle for bnd pnt
    pnts_angles_bnd_i <- pnts_angles_bnd[i]
    #
    # initial
    pnts_angles_sign <- pnts_angles - pnts_angles_bnd_i
    #
    # update
    pnts_angles_sign[pnts_angles_sign > pi] <- pnts_angles_sign[pnts_angles_sign > pi] - pi*2
    pnts_angles_sign[pnts_angles_sign < -pi] <- pnts_angles_sign[pnts_angles_sign < -pi] + pi*2
    #
    # update
    pnts_angles_sign[pnts_angles_sign == pi] <- 0
    pnts_angles_sign[pnts_angles_sign == -pi] <- 0
    #
    # get count
    count_right <- sum(pnts_angles_sign < 0)
    count_left <- sum(pnts_angles_sign > 0)
    #
    # get count diff
    column_count_diff <- c(column_count_diff, count_right - count_left)
  }
  #
  # initial
  column_unblocked <- c()
  #
  # check blocked or not
  for (i in 1:count_pnt_bnd) {
    #
    # get bnd pnt
    pnt_bnd <- pnts_curve_bnd[i, ]
    #
    # get intersects
    ints <- f_pnts_split_equal_fromApnt_ints(pnts_curve_bnd, lpnt0, pnt_bnd)
    #
    # if, or not, blocked by another bnd pnt
    if (nrow(ints) == 1) { column_unblocked <- c(column_unblocked, 1) }
    else { column_unblocked <- c(column_unblocked, 0) }
  }
  #
  # initial ID split
  ID_split_right <- NA
  ID_split_left <- NA
  #
  # get minimum absolute difference
  diff_abs_min <- min(abs(column_count_diff[column_unblocked == 1]))
  diff_abs_min_ID <- which((abs(column_count_diff) == diff_abs_min) & (column_unblocked == 1))
  #
  # check
  if (length(diff_abs_min_ID) == 0) {
    #
    # possibly
    # mostly because end group is too small
    # for example, two bnd pnts, one internal pnts
    # mark stop marker for group as "SGE"
    #
    # return
    return(NA)
  }
  #
  # check
  if (column_count_diff[diff_abs_min_ID[1]] == 0) {
    #
    # check
    if (length(diff_abs_min_ID) > 1) {
      #
      # possibly
      # mostly because end group is too small
      # for example, lpnt0 is circled by pnts
      # then, any segment can split equally
      # mark stop marker for group as "SGE"
      #
      # return
      return(NA)
    }
    #
    # get
    ID_split_right <- diff_abs_min_ID
    ID_split_left <- diff_abs_min_ID
  }
  else {
    #
    # initial
    column_count_diff_around_zero <- c()
    #
    # get
    for (i in 1:(count_pnt_bnd-1)) { 
      #
      # append
      if ((column_unblocked[i] == 1) && (column_unblocked[i+1] == 1) && (column_count_diff[i]*column_count_diff[i+1] < 0)) { 
        #
        # append
        column_count_diff_around_zero <- c(column_count_diff_around_zero, 1) 
        #
        # get
        ID_split_right <- i
        ID_split_left <- i+1
      }
      else {
        #
        # append
        column_count_diff_around_zero <- c(column_count_diff_around_zero, 0)
      }
    }
    #
    # debugging
    if (sum(column_count_diff_around_zero == 1) > 1) {
      #
      # get pnts_curve
      pnts_curve_segment <- rbind(as.numeric(lpnt0))
      #
      # get
      colnames(pnts_curve_segment) <- c("Cx", "Cy")
      pnts_curve_segment <- data.frame(pnts_curve_segment)
      #
      # debugging
      f_pnts_save_points(pnts_curve_bnd, NA, "pnts_curve_bnd")
      f_pnts_save_points(pnts_curve_segment, NA, "pnts_curve_segment")
      #
      # warning
      warning("Warning: more than one pnt pair with difference around zero.")
      #
      # return
      return(NULL)
    }
    #
    # do not have one pair around zero, use min
    if (sum(column_count_diff_around_zero == 1) < 1) {
      #
      # get
      ID_split_right <- diff_abs_min_ID
      ID_split_left <- diff_abs_min_ID
    }
  }
  #
  # get pnts bnd sides
  # both sides keep split pnt?
  # not necessarily real left and right, just for a mark here
  pnts_bnd_sides_right <- pnts_lpnt0_rm_bnd_ascending[1:ID_split_right, ]
  pnts_bnd_sides_left <- pnts_lpnt0_rm_bnd_ascending[ID_split_left:nrow(pnts_lpnt0_rm_bnd_ascending), ]
  #
  # return
  return(list(pnts_bnd_sides_right, pnts_bnd_sides_left))
}
#
#
#
f_pnts_split_equal_fromApnt_ints <- function(pnts_curve_bnd, lpnt0, pnt_bnd) {
  #
  # get pnts_curve
  pnts_curve_segment <- rbind(as.numeric(lpnt0), as.numeric(pnt_bnd[1, c("Cx", "Cy")]))
  # #
  # # stop at the bnd pnt
  # # only check blocked, do not check overreached?
  # #
  # # get limits
  # Xmin <- min(pnts_curve_bnd[, "Cx"])
  # Xmax <- max(pnts_curve_bnd[, "Cx"])
  # Ymin <- min(pnts_curve_bnd[, "Cy"])
  # Ymax <- max(pnts_curve_bnd[, "Cy"])
  # #
  # # initial, landing pnt of segment
  # landing_x <- lpnt0[1]
  # landing_y <- lpnt0[2]
  # #
  # # get bnd pnt, xy
  # pnt_bnd_x <- pnt_bnd[1, c("Cx")]
  # pnt_bnd_y <- pnt_bnd[1, c("Cy")]
  # #
  # # get, landing pnt
  # if (pnt_bnd_x == lpnt0[1]) {
  #   #
  #   # update
  #   landing_x <- lpnt0[1]
  #   #
  #   # update
  #   if (pnt_bnd_y < lpnt0[2]) { landing_y <- Ymin }
  #   else { landing_y <- Ymax }
  #   #
  # }
  # else {
  #   #
  #   # update
  #   if (pnt_bnd_x < lpnt0[1]) { landing_x <- Xmin }
  #   else { landing_x <- Xmax }
  #   #
  #   # update
  #   landing_y <- (pnt_bnd_y - lpnt0[2]) / (pnt_bnd_x - lpnt0[1]) * (landing_x - lpnt0[1]) + lpnt0[2]
  # }
  # #
  # # get pnts_curve
  # pnts_curve_segment <- rbind(as.numeric(lpnt0), as.numeric(c(landing_x, landing_y)))
  #
  # get
  colnames(pnts_curve_segment) <- c("Cx", "Cy")
  pnts_curve_segment <- data.frame(pnts_curve_segment)
  #
  # get intersects
  ints <- f_intersection_curves(pnts_curve_bnd, pnts_curve_segment)
  #
  # debugging
  if (is.null(ints)) {
    #
    # debugging
    f_pnts_save_points(pnts_curve_bnd, NA, "pnts_curve_bnd")
    f_pnts_save_points(pnts_curve_segment, NA, "pnts_curve_segment")
    #
    # warning
    warning("Warning: no intscts gotten, in equally splitting pnts.")
    #
    # return
    return(NULL)
  }
  #
  # return
  return(ints)
}
#
#
#
f_pnts_path <- function(list_pnts, IDB_max, pPolygons) {
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
    # get pnts_bnd
    pnts_bnd <- pnts[which(pnts[, "IDB"] != 0), ]
    #
    # get launch anchor
    LaunchAnchor <- pnts_bnd[which.max(pnts_bnd[, "Cz"]), ]
    #
    # remove launch anchor
    pnts_rm_LaunchAnchor <- pnts[-which(pnts[, "ID"] == LaunchAnchor[, "ID"]), ]
    #
    # get bnd sides
    list_bnd_sides <- f_pnts_bnd_sides(pnts_rm_LaunchAnchor, IDB_max, LaunchAnchor[, c("Cx", "Cy")])
    #
    # switch, do not need here?
    list_bnd_sides <- f_pnts_bnd_sides_switch2(list_bnd_sides, LaunchAnchor, IDB_max)
    #
    # get centerLs
    centerLs <- f_pnts_bnd_sides_centerLs(list_bnd_sides)
    #
    # get anchor initial
    anchor_initial <- c(centerLs[1, ], f_pnts_zip(pnts, centerLs[1, 1], centerLs[1, 2]))
    #
    # get anchor distal
    anchor_distal <- c(centerLs[2, ], f_pnts_zip(pnts, centerLs[2, 1], centerLs[2, 2]))
    #
    # get anchor mid
    anchor_mid <- c(centerLs[3, ], f_pnts_zip(pnts, centerLs[3, 1], centerLs[3, 2]))
    #
    # get type
    column_type <- data.frame("Anchor: Initial group anchor", "Anchor: Group center", "Anchor: Distal group anchor")
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
    # get bnd sides
    bnd_sides_right <- list_bnd_sides[[1]]
    bnd_sides_left <- list_bnd_sides[[2]]
    #
    # update initial and distal anchors
    # let them connect to polygon boundary
    #
    # set a switch
    if (T) {
      #
      # check
      if (nrow(pnts_bnd) > 3) {
        #
        # update
        anchor_xy_initial <- f_pnts_path_extending(column_anchors[c(2, 1), c("Cx", "Cy")], pPolygons)
        column_anchors[1, c("Cx", "Cy", "Cz")] <- c(anchor_xy_initial, f_pnts_zip(pnts, anchor_xy_initial[1], anchor_xy_initial[2]))
        #
        # update
        anchor_xy_distal <- f_pnts_path_extending(column_anchors[c(2, 3), c("Cx", "Cy")], pPolygons)
        column_anchors[3, c("Cx", "Cy", "Cz")] <- c(anchor_xy_distal, f_pnts_zip(pnts, anchor_xy_distal[1], anchor_xy_distal[2]))
      }
      else {
        #
        # update
        anchor_xy_initial <- f_pnts_path_extending(rbind(column_anchors[2, c("Cx", "Cy")], LaunchAnchor[1, c("Cx", "Cy")]), pPolygons)
        column_anchors[1, c("Cx", "Cy", "Cz")] <- c(anchor_xy_initial, f_pnts_zip(pnts, anchor_xy_initial[1], anchor_xy_initial[2]))
        #
        # update
        anchor_xy_distal <- f_pnts_path_extending(rbind(LaunchAnchor[1, c("Cx", "Cy")], column_anchors[2, c("Cx", "Cy")]), pPolygons)
        column_anchors[3, c("Cx", "Cy", "Cz")] <- c(anchor_xy_distal, f_pnts_zip(pnts, anchor_xy_distal[1], anchor_xy_distal[2]))
        #
        # update
        anchor_xy_mid <- c((anchor_xy_initial[1]+anchor_xy_distal[1])/2, (anchor_xy_initial[2]+anchor_xy_distal[2])/2)
        column_anchors[2, c("Cx", "Cy", "Cz")] <- c(anchor_xy_mid, f_pnts_zip(pnts, anchor_xy_mid[1], anchor_xy_mid[2]))
        #
        # update
        # do not update bnd_sides_right and bnd_sides_left?
      }
    }
    #
    # return
    return(list(column_anchors, list(bnd_sides_right), list(bnd_sides_left)))
  }
  #
  # get pnts
  pnts <- f_list2pnts(list_pnts)
  #
  # initial, for distinguishing right and left sides
  LaunchAnchor <- c()
  #
  # initial
  list_bnd_sides_right <- list()
  list_bnd_sides_left <- list()
  #
  # initial
  column_centerL <- c()
  #
  # get centers
  for (i in 1:count_pnts) {
    #
    # get pnts_i
    pnts_i <- list_pnts[[i]]
    #
    # get launch anchor for end groups
    pnts_end_LaunchAnchor <- c()
    if (i == 1 || i == count_pnts) { pnts_end_LaunchAnchor <- f_pnts_end_LaunchAnchor(pnts_i) }
    #
    # get bnd sides
    list_pnts_i_bnd_sides <- f_pnts_bnd_sides(pnts_i, IDB_max, pnts_end_LaunchAnchor)
    #
    # get launch anchor
    if (i == 1) {
      #
      # update launch anchor
      LaunchAnchor <- f_pnts_end_LaunchAnchor_initial(list_pnts_i_bnd_sides, IDB_max)
    }
    #
    # switch
    list_pnts_i_bnd_sides <- f_pnts_bnd_sides_switch2(list_pnts_i_bnd_sides, LaunchAnchor, IDB_max)
    # #
    # # debugging
    # f_pnts_save_points(LaunchAnchor, NA, "LaunchAnchor")
    # f_pnts_save_points(list_pnts_i_bnd_sides[[1]], NA, "pnts_i_bnd_sides_right")
    # f_pnts_save_points(list_pnts_i_bnd_sides[[2]], NA, "pnts_i_bnd_sides_left")
    #
    # get center
    if (length(list_pnts_i_bnd_sides) == 2) {
      #
      # get centerLs
      centerLs <- f_pnts_bnd_sides_centerLs(list_pnts_i_bnd_sides)
      centerL1 <- c(centerLs[1, ], f_pnts_zip(pnts_i, centerLs[1, 1], centerLs[1, 2]))
      centerL2 <- c(centerLs[2, ], f_pnts_zip(pnts_i, centerLs[2, 1], centerLs[2, 2]))
      #
      # get center of line
      column_centerL <- rbind(column_centerL, centerL1, centerL2)
    }
    else {
      #
      # error
      warning("Error: pnts do not have 2 bnd sides.")
      #
      # return
      return(NULL)
    }
    #
    # append
    list_bnd_sides_right[[i]] <- list_pnts_i_bnd_sides[[1]]
    list_bnd_sides_left[[i]] <- list_pnts_i_bnd_sides[[2]]
  }
  #
  # initial
  column_anchors <- c()
  #
  # initial type
  column_type <- c()
  #
  # initial grp
  column_grp <- c()
  #
  # update centerG
  for (i in 1:(nrow(column_centerL)-1)) {
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
    column_anchors <- rbind(column_anchors, c(centerLx, centerLy, centerLz))
    #
    # type and grp
    if (i%%2 == 1) {
      #
      # odd
      column_type <- rbind(column_type, "Anchor: Group center")
      column_grp <- rbind(column_grp, (i-1)/2 + 1)
    }
    else {
      #
      # even
      column_type <- rbind(column_type, "Anchor: Inter-group center")
      column_grp <- rbind(column_grp, 0)
    }
  }
  #
  # must do this, or numeric will become character
  column_anchors <- data.frame(column_anchors)
  column_type <- data.frame(column_type)
  column_grp <- data.frame(column_grp)
  #
  # initial
  column_anchors <- cbind(column_anchors, column_type, column_grp)
  column_anchors <- data.frame(column_anchors)
  #
  # colnames
  colnames(column_anchors) <- c("Cx", "Cy", "Cz", "type", "grp")
  #
  # get anchor_end_initial
  anchor_end_initial <- cbind(data.frame(column_centerL)[1, ], "Anchor: Initial group anchor", 0)
  colnames(anchor_end_initial) <- c("Cx", "Cy", "Cz", "type", "grp")
  #
  # update anchors
  column_anchors <- rbind(anchor_end_initial, column_anchors)
  #
  # get anchor_end_distal
  anchor_end_distal <- cbind(data.frame(column_centerL)[nrow(column_centerL), ], "Anchor: Distal group anchor", 0)
  colnames(anchor_end_distal) <- c("Cx", "Cy", "Cz", "type", "grp")
  #
  # update anchors
  column_anchors <- rbind(column_anchors, anchor_end_distal)
  #
  # update initial and distal anchors
  # let them connect to polygon boundary
  #
  # set a switch
  if (T) {
    #
    # update
    anchor_xy_initial <- f_pnts_path_extending(column_anchors[c(2, 1), c("Cx", "Cy")], pPolygons)
    column_anchors[1, c("Cx", "Cy", "Cz")] <- c(anchor_xy_initial, f_pnts_zip(list_pnts[[1]], anchor_xy_initial[1], anchor_xy_initial[2]))
    #
    # get count
    count_anchor <- nrow(column_anchors)
    #
    # update
    anchor_xy_distal <- f_pnts_path_extending(column_anchors[c(count_anchor-1, count_anchor), c("Cx", "Cy")], pPolygons)
    column_anchors[count_anchor, c("Cx", "Cy", "Cz")] <- c(anchor_xy_distal, f_pnts_zip(list_pnts[[count_pnts]], anchor_xy_distal[1], anchor_xy_distal[2]))
  }
  # #
  # # for debugging
  # file_out <- paste("column_anchors.xlsx", sep = "", collapse = NULL)
  # f_pnts_save_xls(column_anchors, file_out)
  # 
  # return
  return(list(data.frame(column_anchors), list_bnd_sides_right, list_bnd_sides_left))
}
#
#
#
f_pnts_path_extending <- function(pnts, pPolygons) {
  #
  # get coords of polygons
  pPolygons_coords <- pPolygons@Polygons[[1]]@coords
  #
  # get limits of polygons
  XMin <- min(pPolygons_coords[, 1])
  XMax <- max(pPolygons_coords[, 1])
  YMin <- min(pPolygons_coords[, 2])
  YMax <- max(pPolygons_coords[, 2])
  #
  # get
  pnt1_x <- pnts[1, 1]
  pnt1_y <- pnts[1, 2]
  pnt2_x <- pnts[2, 1]
  pnt2_y <- pnts[2, 2]
  #
  # initial
  pnt0_x <- c()
  pnt0_y <- c()
  #
  # check
  if (pnt1_x == pnt2_x) {
    #
    # get x
    pnt0_x <- pnt1_x
    #
    # get y
    if (pnt2_y > pnt1_y) { pnt0_y <- YMax}
    else { pnt0_y <- YMin }
  }
  else {
    #
    # get x
    if (pnt2_x > pnt1_x) { pnt0_x <- XMax}
    else { pnt0_x <- XMin }
    #
    # get y
    pnt0_y <- pnt1_y + (pnt0_x - pnt1_x) * (pnt2_y - pnt1_y) / (pnt2_x - pnt1_x)
  }
  #
  # get pnts, line
  pnts_line <- matrix(rbind(c(pnt1_x, pnt1_y), c(pnt0_x, pnt0_y)), nc = 2, byrow = F)
  #
  # get ints
  ints <- f_intersection_curves(pnts_line, pPolygons_coords)
  #
  # get nearest
  ints_nearest <- f_pnts_nearest_extending(c(pnt1_x, pnt1_y), ints, Inf)
  #
  # return
  return(ints_nearest)
}
#
#
#
f_pnts_path_lines <- function(list_bnd_sides_right, list_bnd_sides_left, column_anchors, MinStpHrzLen) {
  # #
  # # remove Cz
  # column_anchors <- column_anchors[, ! names(column_anchors) %in% c("Cz"), drop = F]
  #
  # colnames, change Cz to Czp
  colnames(column_anchors) <- c("Cx", "Cy", "Czp", "type", "grp")
  #
  # get pnts_curve
  pnts_curve_profile <- column_anchors[, c("Cx", "Cy")]
  #
  # get count
  if (length(list_bnd_sides_right) != length(list_bnd_sides_left)) { return(NULL) }
  # get
  count_grp <- length(list_bnd_sides_right)
  count_grp_bnd <- (count_grp-2)*2 + 2*1
  # get
  count_anchor <- nrow(column_anchors)
  #
  # initial
  column_anchors_GrpBndInts <- column_anchors[-(1:nrow(column_anchors)), ]
  #
  # initial
  list_GrpBnds <- list()
  #
  # check, if, no group boundary, possibly only one group
  if (count_grp_bnd >= 1) {
    #
    # get intersects
    for (i in 1:count_grp) {
      #
      # get bnd sides for a group
      list_bnd_sides <- list(list_bnd_sides_right[[i]], list_bnd_sides_left[[i]])
      #
      # get bnd sides corners()
      bnd_sides_corners <- f_pnts_bnd_sides_corners(list_bnd_sides)
      #
      # initial
      ints_upper <- c()
      ints_lower <- c()
      #
      # get ints
      if (i != 1 && i != count_grp) {
        #
        # # #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, column_anchors[i*2-1, ])
        #
        # get pnts_curve
        pnts_curve_grp_bnd_upper <- rbind(bnd_sides_corners[4, c("Cx", "Cy")], bnd_sides_corners[1, c("Cx", "Cy")])
        list_GrpBnds[[length(list_GrpBnds)+1]] <- pnts_curve_grp_bnd_upper
        #
        # get intersects
        ints_upper <- f_intersection_curves(pnts_curve_profile, pnts_curve_grp_bnd_upper)
        #
        # update
        ints_upper <- data.frame(ints_upper)
        ints_upper <- cbind(ints_upper, 0, "Intersection: Group boundary upper", i)
        colnames(ints_upper) <- c("Cx", "Cy", "Czp", "type", "grp")
        #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, ints_upper)
        #
        # # #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, column_anchors[i*2, ])
        #
        # get pnts_curve
        pnts_curve_grp_bnd_lower <- rbind(bnd_sides_corners[3, c("Cx", "Cy")], bnd_sides_corners[2, c("Cx", "Cy")])
        list_GrpBnds[[length(list_GrpBnds)+1]] <- pnts_curve_grp_bnd_lower
        #
        # get intersects
        ints_lower <- f_intersection_curves(pnts_curve_profile, pnts_curve_grp_bnd_lower)
        #
        # update
        ints_lower <- data.frame(ints_lower)
        ints_lower <- cbind(ints_lower, 0, "Intersection: Group boundary lower", i)
        colnames(ints_lower) <- c("Cx", "Cy", "Czp", "type", "grp")
        #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, ints_lower)
      }
      else if (i == 1 && i != count_grp) {
        #
        # # #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, column_anchors[1:2, ])
        #
        # get pnts_curve
        pnts_curve_grp_bnd_lower <- rbind(bnd_sides_corners[3, c("Cx", "Cy")], bnd_sides_corners[2, c("Cx", "Cy")])
        list_GrpBnds[[length(list_GrpBnds)+1]] <- pnts_curve_grp_bnd_lower
        #
        # get intersects
        ints_lower <- f_intersection_curves(pnts_curve_profile, pnts_curve_grp_bnd_lower)
        #
        # update
        ints_lower <- data.frame(ints_lower)
        ints_lower <- cbind(ints_lower, 0, "Intersection: Group boundary lower", i)
        colnames(ints_lower) <- c("Cx", "Cy", "Czp", "type", "grp")
        #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, ints_lower)
      }
      else if (i != 1 && i == count_grp) {
        #
        # # #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, column_anchors[i*2-1, ])
        #
        # get pnts_curve
        pnts_curve_grp_bnd_upper <- rbind(bnd_sides_corners[4, c("Cx", "Cy")], bnd_sides_corners[1, c("Cx", "Cy")])
        list_GrpBnds[[length(list_GrpBnds)+1]] <- pnts_curve_grp_bnd_upper
        #
        # get intersects
        ints_upper <- f_intersection_curves(pnts_curve_profile, pnts_curve_grp_bnd_upper)
        #
        # update
        ints_upper <- data.frame(ints_upper)
        ints_upper <- cbind(ints_upper, 0, "Intersection: Group boundary upper", i)
        colnames(ints_upper) <- c("Cx", "Cy", "Czp", "type", "grp")
        #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, ints_upper)
        #
        # # #
        # append
        column_anchors_GrpBndInts <- rbind(column_anchors_GrpBndInts, column_anchors[(i*2):(i*2+1), ])
      }
    }
    #
    # debugging
    if (F) {
      #
      # save
      f_pnts_save_polylines(list_GrpBnds, NA, "list_GrpBnds")
    }
  }
  else {
    #
    # no group boundary
    column_anchors_GrpBndInts <- column_anchors
  }
  #
  # initial column length
  column_anchors_GrpBndInts$Lhrz <- c(0)
  column_anchors_GrpBndInts$Lallp <- c(0)
  #
  # get count
  count_acr_int <- nrow(column_anchors_GrpBndInts)
  #
  # get length 2d
  for (i in 2:count_acr_int) {
    #
    # update
    column_anchors_GrpBndInts[i, "Lhrz"] <- 
      column_anchors_GrpBndInts[i-1, "Lhrz"] + 
      sqrt((column_anchors_GrpBndInts[i, "Cx"] - column_anchors_GrpBndInts[i-1, "Cx"])^2 + 
           (column_anchors_GrpBndInts[i, "Cy"] - column_anchors_GrpBndInts[i-1, "Cy"])^2)
  }
  #
  # get GrpBndInt, index
  index_GrpBndInts <- which(column_anchors_GrpBndInts[, "type"] == "Intersection: Group boundary lower")
  index_GrpBndInts <- c(index_GrpBndInts, which(column_anchors_GrpBndInts[, "type"] == "Intersection: Group boundary upper"))
  #
  # get anchors
  if (length(index_GrpBndInts) == 0) { column_anchors <- column_anchors_GrpBndInts }
  else { column_anchors <- column_anchors_GrpBndInts[-index_GrpBndInts, ] }
  #
  # get GrpBndInt, initial
  column_GrpBndInts <- c()
  #
  # get GrpBndAgls, initial
  column_GrpBndAgls <- c()
  #
  # check, if, no group boundary, possibly only one group
  if (count_grp_bnd >= 1) {
    #
    # get GrpBndInts
    column_GrpBndInts <- column_anchors_GrpBndInts[index_GrpBndInts, ]
    column_GrpBndInts <- column_GrpBndInts[order(column_GrpBndInts$Lhrz), ]
    #
    # get GrpBndAgls
    for (i in 1:count_grp_bnd) {
      #
      # get GrpBnd
      GrpBnd <- list_GrpBnds[[i]]
      #
      # get angle
      GrpBnd_agl <- atan2(GrpBnd[2, "Cy"]-GrpBnd[1, "Cy"], GrpBnd[2, "Cx"]-GrpBnd[1, "Cx"])
      #
      # append
      column_GrpBndAgls <- rbind(column_GrpBndAgls, GrpBnd_agl)
    }
  }
  #
  # get split stations
  list_split_stations <- f_pnts_split_stations(pnts_curve_profile, MinStpHrzLen)
  # get
  strip_station_interval <- list_split_stations[[1]]
  # get
  strip_stations <- list_split_stations[[2]]
  #
  # initial
  strip_stations_angles <- c()
  # initial
  strip_stations_Lhrzs <- c()
  strip_stations_Lallps <- c()
  #
  # initial
  column_anchors_GrpBndInts_stations <- column_anchors_GrpBndInts
  column_anchors_GrpBndInts_stations$station <- c(0)
  #
  # initial
  strip_angle_for_one_group <- NA
  #
  # check, if, no group boundary, possibly only one group
  if (count_grp_bnd == 0) {
    #
    # get 
    dx <- column_anchors[3, "Cx"] - column_anchors[1, "Cx"]
    dy <- column_anchors[3, "Cy"] - column_anchors[1, "Cy"]
    # get angle
    strip_angle_for_one_group <- atan2(dy, dx) - pi/2
    # update
    if (strip_angle_for_one_group > 2*pi) { strip_angle_for_one_group <- strip_angle_for_one_group - 2*pi }
    if (strip_angle_for_one_group < -2*pi) { strip_angle_for_one_group <- strip_angle_for_one_group + 2*pi }
  }
  #
  # get count
  count_station <- nrow(strip_stations)
  #
  # get strip stations
  for (i in 1:count_station) {
    #
    # get length at station
    strip_station_Lhrz <- strip_station_interval * i
    #
    # append
    strip_stations_Lhrzs <- rbind(strip_stations_Lhrzs, strip_station_Lhrz)
    #
    # initial
    strip_angle <- strip_angle_for_one_group
    #
    # check, if, no group boundary, possibly only one group
    if (count_grp_bnd >= 1) {
      #
      # initial
      index_GrpBndInt_upper <- count_grp_bnd
      index_GrpBndInt_lower <- count_grp_bnd+1
      #
      # get index
      for (k in 1:count_grp_bnd) {
        #
        if (column_GrpBndInts[k, "Lhrz"] > strip_station_Lhrz) {
          #
          index_GrpBndInt_lower <- k
          index_GrpBndInt_upper <- k-1
          #
          break
        }
      }
      #
      # if within the initial group
      if (index_GrpBndInt_upper == 0) {
        #
        # get strip angle
        strip_angle <- column_GrpBndAgls[index_GrpBndInt_lower]
      }
      #
      # if within the distal group
      if (index_GrpBndInt_lower == count_grp_bnd+1) {
        #
        # get strip angle
        strip_angle <- column_GrpBndAgls[index_GrpBndInt_upper]
      }
      #
      # if within
      if (index_GrpBndInt_upper != 0 && index_GrpBndInt_lower != count_grp_bnd+1) {
        #
        # get length difference
        Lhrz_diff <- column_GrpBndInts[index_GrpBndInt_lower, "Lhrz"] - column_GrpBndInts[index_GrpBndInt_upper, "Lhrz"]
        Lhrz_diff_station <- strip_station_Lhrz - column_GrpBndInts[index_GrpBndInt_upper, "Lhrz"]
        #
        # get angle
        angle_grp_bnd_upper <- column_GrpBndAgls[index_GrpBndInt_upper]
        angle_grp_bnd_lower <- column_GrpBndAgls[index_GrpBndInt_lower]
        #
        # get angle difference
        angle_diff <- angle_grp_bnd_lower - angle_grp_bnd_upper
        if (angle_diff > 2*pi) { angle_diff <- angle_diff - 2*pi }
        if (angle_diff < -2*pi) { angle_diff <- angle_diff + 2*pi }
        #
        # get angle difference at the station
        angle_diff_station <- angle_diff / Lhrz_diff * Lhrz_diff_station
        #
        # get strip angle
        strip_angle <- angle_grp_bnd_upper + angle_diff_station
        if (strip_angle > pi) { strip_angle <- strip_angle - 2*pi }
        if (strip_angle < -pi) { strip_angle <- strip_angle + 2*pi }
      }
    }
    #
    # append
    strip_stations_angles <- rbind(strip_stations_angles, strip_angle)
    #
    # get station
    station <- cbind(strip_stations[i, ], 0, "Station: Strip boundary", 0, strip_station_Lhrz, 0, i)
    colnames(station) <- c("Cx", "Cy", "Czp", "type", "grp", "Lhrz", "Lallp", "station")
    #
    # get count
    count_acr_int_stn <- nrow(column_anchors_GrpBndInts_stations)
    #
    # initial
    index_anchor_GrpBndInt_station_upper <- c()
    #
    # get index
    for (k in 1:count_acr_int_stn) {
      #
      if (column_anchors_GrpBndInts_stations[k, "Lhrz"] > strip_station_Lhrz) {
        #
        index_anchor_GrpBndInt_station_upper <- k-1
        #
        break
      }
    }
    #
    # update
    if (index_anchor_GrpBndInt_station_upper >= 1) {
      #
      # insert
      column_anchors_GrpBndInts_stations <- 
        rbind(column_anchors_GrpBndInts_stations[1:index_anchor_GrpBndInt_station_upper, ], 
              station, 
              column_anchors_GrpBndInts_stations[(index_anchor_GrpBndInt_station_upper+1):count_acr_int_stn, ])
    }
    else {
      #
      # insert
      column_anchors_GrpBndInts_stations <- 
        rbind(station, 
              column_anchors_GrpBndInts_stations[(index_anchor_GrpBndInt_station_upper+1):count_acr_int_stn, ])
    }
  }
  #
  # get count
  count_acr_int_stn <- nrow(column_anchors_GrpBndInts_stations)
  #
  # get Czp (Cz of plane fitted to pnts)
  for (i in 1:count_acr_int_stn) {
    #
    # check
    if (column_anchors_GrpBndInts_stations[i, "type"] == "Station: Strip boundary" ||
        column_anchors_GrpBndInts_stations[i, "type"] == "Intersection: Group boundary lower" ||
        column_anchors_GrpBndInts_stations[i, "type"] == "Intersection: Group boundary upper") {
      #
      # get length horizontal
      Lhrz <- column_anchors_GrpBndInts_stations[i, "Lhrz"]
      #
      # initial
      index_anchor_upper <- NA
      index_anchor_lower <- NA
      #
      # get index
      for (k in 1:count_anchor) {
        #
        if (column_anchors[k, "Lhrz"] > Lhrz) {
          #
          index_anchor_lower <- k
          index_anchor_upper <- k-1
          #
          break
        }
      }
      #
      # get length ratio
      length_ratio <- 
        (Lhrz - column_anchors[index_anchor_upper, "Lhrz"]) /
        (column_anchors[index_anchor_lower, "Lhrz"] - column_anchors[index_anchor_upper, "Lhrz"])
      #
      # get Czp
      column_anchors_GrpBndInts_stations[i, "Czp"] <- 
        column_anchors[index_anchor_upper, "Czp"] + 
        length_ratio *
        (column_anchors[index_anchor_lower, "Czp"] - column_anchors[index_anchor_upper, "Czp"])
    }
  }
  #
  # get length 3d
  for (i in 2:count_acr_int_stn) {
    #
    # get
    Lallp <- 
      column_anchors_GrpBndInts_stations[i-1, "Lallp"] + 
      sqrt((column_anchors_GrpBndInts_stations[i, "Cx"] - column_anchors_GrpBndInts_stations[i-1, "Cx"])^2 + 
             (column_anchors_GrpBndInts_stations[i, "Cy"] - column_anchors_GrpBndInts_stations[i-1, "Cy"])^2 + 
             (column_anchors_GrpBndInts_stations[i, "Czp"] - column_anchors_GrpBndInts_stations[i-1, "Czp"])^2)
    #
    # update
    column_anchors_GrpBndInts_stations[i, "Lallp"] <- Lallp
    #
    # check
    if (column_anchors_GrpBndInts_stations[i, "type"] == "Station") {
      #
      # append
      strip_stations_Lallps <- rbind(strip_stations_Lallps, Lallp)
    }
  }
  #
  # switch
  column_anchors_GrpBndInts_stations <- column_anchors_GrpBndInts_stations[, c("Cx", "Cy", "Czp", "type", "grp", "station", "Lhrz", "Lallp")]
  #
  # debugging
  if (F) {
    #
    # save
    f_pnts_save_points(column_anchors, NA, "column_anchors")
    f_pnts_save_points(column_anchors_GrpBndInts, NA, "column_anchors_GrpBndInts")
    f_pnts_save_points(column_anchors_GrpBndInts_stations, NA, "column_anchors_GrpBndInts_stations")
    f_pnts_save_points(column_GrpBndInts, NA, "column_GrpBndInts")
  }
  #
  # get total length
  Lhrz_total <- column_anchors_GrpBndInts_stations[count_acr_int_stn, "Lhrz"]
  Lallp_total <- column_anchors_GrpBndInts_stations[count_acr_int_stn, "Lallp"]
  #
  # return
  return(list(list_GrpBnds, strip_stations, strip_stations_angles, strip_stations_Lhrzs, strip_stations_Lallps, Lhrz_total, Lallp_total, column_anchors_GrpBndInts_stations))
}
#
#
#
f_pnts_split_station <- function(pnts, length_cumus, length_split) {
  #
  # check
  stopifnot((length_split >= 0 && length_split <= length_cumus[length(length_cumus)]))
  #
  # get 
  for (i in 1:length(length_cumus)) {
    #
    # check split position
    if (length_split >= length_cumus[i] && length_split <= length_cumus[i+1]) {
      #
      # get split position ratio
      length_diff_ratio <- (length_split - length_cumus[i]) / (length_cumus[i+1] - length_cumus[i])
      #
      # get split position coordinates
      x <- pnts[i, 1] + (pnts[i+1, 1] - pnts[i, 1]) * length_diff_ratio
      y <- pnts[i, 2] + (pnts[i+1, 2] - pnts[i, 2]) * length_diff_ratio
      #
      # return
      return(data.frame(Cx = x, Cy = y))
    }
  }
}
#
#
#
f_pnts_split_stations <- function(pnts, MinStpHrzLen = 0) {
  #
  # check
  stopifnot(MinStpHrzLen > 0)
  #
  # get count
  count_pnt <- nrow(pnts)
  #
  # initial
  length_cumus <- rbind(c(0))
  #
  # get length
  for (i in 2:count_pnt) {
    #
    # get
    length_i <- sqrt((pnts[i, 1] - pnts[i-1, 1])^2 + (pnts[i, 2] - pnts[i-1, 2])^2)
    #
    # append
    length_cumus <- rbind(length_cumus, length_cumus[nrow(length_cumus)] + length_i)
  }
  #
  # get total length
  length_total <- length_cumus[count_pnt]
  #
  # initial
  count_split <- 0
  #
  # get count
  while (T) {
    #
    # update
    if (length_total/(count_split+1) >= MinStpHrzLen) { count_split <- count_split+1 }
    else { break }
  }
  #
  # get stationing of segments
  stationing <- seq(from = 0, to = length_total, length.out = count_split+1)
  #
  # check, do not allow no station
  if (count_split == 1) {
    #
    # update
    count_split <- 2
    #
    # get stationing of segments
    stationing <- seq(from = 0, to = length_total, length.out = count_split+1)
  }
  #
  # initial
  column_split_stations <- c()
  #
  # get
  for (i in 2:count_split) {
    #
    # get
    split_station <- f_pnts_split_station(pnts, length_cumus, stationing[i])
    #
    # append
    column_split_stations <- rbind(column_split_stations, split_station)
  }
  #
  # return
  return(list(length_total/count_split, column_split_stations))
}
#
#
#
f_strips <- function(pRasterDEM, pPolygons, list_path, MinStpHrzLen) {
  #
  # get anchors
  column_anchors <- list_path[[1]]
  # get bnd sides
  list_bnd_sides_right <- list_path[[2]]
  list_bnd_sides_left <- list_path[[3]]
  #
  # get path lines
  list_path_lines <- f_pnts_path_lines(list_bnd_sides_right, list_bnd_sides_left, column_anchors, MinStpHrzLen)
  #
  # get 
  strip_GrpBnds <- list_path_lines[[1]]
  # get
  strip_stations <- list_path_lines[[2]]
  strip_stations_angles <- list_path_lines[[3]]
  # get 
  strip_stations_Lhrzs <- list_path_lines[[4]]
  strip_stations_Lallps <- list_path_lines[[5]]
  # get 
  Lhrz_total <- list_path_lines[[6]]
  Lallp_total <- list_path_lines[[7]]
  # get 
  strip_nodes <- list_path_lines[[8]]
  #
  # initial
  strip_nodes$Czd <- c(0)
  # initial
  strip_nodes$Lalld <- c(0)
  #
  # get count
  count_node <- nrow(strip_nodes)
  #
  # get Czd (Cz of DEM)
  for (i in 1:count_node) {
    #
    # get
    strip_nodes[i, "Czd"] <- raster::extract(pRasterDEM, strip_nodes[i, c("Cx", "Cy")])
  }
  #
  # get length 3d
  for (i in 2:count_node) {
    #
    # get
    Lalld <- 
      strip_nodes[i-1, "Lalld"] + 
      sqrt((strip_nodes[i, "Cx"] - strip_nodes[i-1, "Cx"])^2 + 
             (strip_nodes[i, "Cy"] - strip_nodes[i-1, "Cy"])^2 + 
             (strip_nodes[i, "Czd"] - strip_nodes[i-1, "Czd"])^2)
    #
    # update
    strip_nodes[i, "Lalld"] <- Lalld
  }
  #
  # switch
  strip_nodes <- strip_nodes[, c("Cx", "Cy", "Czd", "Czp", "type", "grp", "station", "Lhrz", "Lalld", "Lallp")]
  #
  # initial
  strip_GrpBnds_extnd <- list()
  #
  # get count
  count_GrpBnd <- length(strip_GrpBnds)
  #
  # check, possibly only one group, no group boundary
  if (count_GrpBnd >= 1) {
    #
    # get GrpBnd extending
    for (i in 1:count_GrpBnd) {
      #
      # get
      GrpBnd <- strip_GrpBnds[[i]]
      #
      # get extnd
      GrpBnd_extnd <- f_pnts_grp_bnd_extending(pPolygons, GrpBnd)
      #
      # append
      strip_GrpBnds_extnd[[i]] <- GrpBnd_extnd
    }
  }
  #
  # get clipping lines
  # debugging, try catch
  if (F) {
    #
    # initial
    out <- NA
    #
    # try catch
    out <- tryCatch (
      #
      ########################################################
      # Try part: define the expression(s) you want to "try" #
      ########################################################
      #
      {
        #
        # Just to highlight: 
        # If you want to use more than one R expression in the "try part" 
        # then you'll have to use curly brackets. 
        # Otherwise, just write the single expression you want to try and 
        #
        message("This is the 'try' part")
        spldf_cls <- f_strips_clipping_polylines(pPolygons, strip_stations, strip_stations_angles)
      },
      #
      ########################################################################
      # Condition handler part: define how you want conditions to be handled #
      ########################################################################
      #
      # Handler when a warning occurs:
      warning = function(cond) {
        #
        message("Get warning:")
        message("Here's the original warning message:")
        message(cond)
        #
        # Choose a return value when such a type of condition occurs
        return(NULL)
      },
      #
      # Handler when an error occurs:
      error = function(cond) {
        #
        # message
        message("Get erros:")
        message("Here's the original error message:")
        message(cond)
        #
        # Choose a return value when such a type of condition occurs
        return(NULL)
      },
      #
      ###############################################
      # Final part: define what should happen AFTER #
      # everything has been tried and/or handled    #
      ###############################################
      #
      finally = {
        #
        # if succeed
        message("Process succeed\n")
      }
    )
    #
    # pause
    if (is.null(out)) {
      #
      # pause at here
      return(NA)
    }
  }
  #
  # get clipping lines
  spldf_cls <- f_strips_clipping_polylines(pPolygons, strip_stations, strip_stations_angles)
  #
  # debugging, output files
  if (F) {
    #
    # write
    rgdal::writeOGR(spldf_cls, dsn = ".", layer = "strip_cls", overwrite = T, driver="ESRI Shapefile")
    #
    # save
    f_pnts_save_polylines(strip_GrpBnds, NA, "strip_GrpBnds")
    #
    # save
    f_pnts_save_polylines(strip_GrpBnds_extnd, NA, "strip_GrpBnds_extnd")
    #
    # save
    f_pnts_save_points(strip_stations, NA, "strip_stations")
  }
  #
  # get clipping lines, sfc
  sfc_spldf_cls <- sf::st_as_sfc(spldf_cls)
  #
  # get SpatialPolygons
  pSpatialPolygons <- sp::SpatialPolygons(list(pPolygons))
  #
  # initial
  sfc_SpatialPolygons0 <- sf::st_as_sfc(pSpatialPolygons)
  sfc_SpatialPolygons <- sfc_SpatialPolygons0
  #
  # initial
  pnt_upper <- column_anchors[1, c("Cx", "Cy")]
  #
  # initial
  list_strip_sp <- list()
  #
  # initial
  strip_Lhrz <- c()
  #
  # initial
  count_strip <- 1
  #
  # initial, there might be invalid clipping polyline
  strip_nodes_for_sp <- strip_nodes
  #
  # get count
  count_cl <- nrow(strip_stations)
  #
  # get SpatialPolygons
  for (i in 1:count_cl) {
    #
    # get clipping polyline
    sfc_cl <- sfc_spldf_cls[[i]]
    #
    # get coords
    sfc_cl_coords <- matrix(sfc_cl, nc = 2, byrow = F)
    #
    # initial
    bool_ints <- FALSE
    #
    # check, possibly only one group, no group boundary
    if (count_GrpBnd >= 1) {
      #
      # check
      for (k in 1:length(strip_GrpBnds_extnd)) {
        #
        # get ints
        ints <- f_intersection_curves(sfc_cl_coords, strip_GrpBnds_extnd[[k]])
        #
        # check
        if (!is.null(ints)) {
          #
          # set
          bool_ints <- TRUE
          break
        }
      }
    }
    #
    # check
    if (bool_ints == TRUE) {
      #
      # append
      strip_nodes_for_sp <- strip_nodes_for_sp[-which(strip_nodes_for_sp[, "station"] == i), ]
    }
    #
    # check
    if (bool_ints == FALSE) {
      #
      # split using line
      sfc_st_split <- lwgeom::st_split(sfc_SpatialPolygons, sfc_cl)
      sfc_st_split_POLYGONs <- st_collection_extract(sfc_st_split, type = c("POLYGON"))
      #
      # get count
      count_POLYGON <- length(sfc_st_split_POLYGONs)
      #
      # debugging
      if (count_POLYGON == 1 || count_POLYGON >= 3) {
        #
        # save
        f_save_list2sp(list_strip_sp, NULL, NA, "list_strip_cp")
        #
        # save spg
        coords <- ((sfc_SpatialPolygons[1])[[1]])[[1]]
        sfc_spg <- f_coords2spg(coords)
        # save
        f_save_list2sp(list(sfc_spg), NULL, NA, "sfc_spg")
        #
        # get spl
        coords <- sfc_cl_coords
        sfc_spl <- f_coords2spl(coords)
        # save
        f_save_list2sp(list(sfc_spl), NULL, NA, "sfc_spl")
        #
        # warning
        warning("Error: count of POLYGON is not 2.")
        #
        # pause at here
        return(NA)
      }
      #
      # get sfc POLYGON
      sfc_POLYGON1 <- sfc_st_split_POLYGONs[1]
      sfc_POLYGON2 <- sfc_st_split_POLYGONs[2]
      #
      # get sfc POLYGON sp
      sfc_POLYGON1_sp <- f_coords2spg((sfc_POLYGON1[[1]])[[1]])
      sfc_POLYGON2_sp <- f_coords2spg((sfc_POLYGON2[[1]])[[1]])
      #
      # get clipping polyline, station
      line_station <- strip_stations[i, 1:2]
      #
      # slightly move forward, or the sp:over will not return right results
      pnt_upper_fwd_x <- pnt_upper[1] + (line_station[1] - pnt_upper[1]) * 10^-6
      pnt_upper_fwd_y <- pnt_upper[2] + (line_station[2] - pnt_upper[2]) * 10^-6
      pnt_upper_fwd <- cbind(pnt_upper_fwd_x, pnt_upper_fwd_y)
      #
      # get shp
      pnt_upper_fwd_shp <- f_pnts2shp_points(pnt_upper_fwd, NA)
      #
      # get over
      sfc_POLYGON1_sp_over <- sp::over(pnt_upper_fwd_shp, sfc_POLYGON1_sp)
      sfc_POLYGON2_sp_over <- sp::over(pnt_upper_fwd_shp, sfc_POLYGON2_sp)
      #
      # check
      if ((!(is.na(sfc_POLYGON1_sp_over)) && !(is.na(sfc_POLYGON2_sp_over))) ||
          (is.na(sfc_POLYGON1_sp_over) && is.na(sfc_POLYGON2_sp_over))) {
        #
        # save
        f_pnts_save_points(pnt_upper, NA, "pnt_upper")
        #
        # save
        f_pnts_save_points(pnt_upper_fwd, NA, "pnt_upper_fwd")
        #
        # save
        f_save_list2sp(list(sfc_POLYGON1_sp, sfc_POLYGON2_sp), NULL, NA, "POLYGONs")
        #
        # warning
        warning("Error: pnt upper is in or not in both POLYGONs.")
        #
        # return
        return(NULL)
      }
      #
      # check
      if (!(is.na(sfc_POLYGON1_sp_over))) {
        #
        # append
        list_strip_sp[[count_strip]] <- sfc_POLYGON1_sp
        count_strip <- count_strip + 1
        #
        # update
        sfc_SpatialPolygons <- sfc_POLYGON2
      }
      #
      # check
      if (!(is.na(sfc_POLYGON2_sp_over))) {
        #
        # append
        list_strip_sp[[count_strip]] <- sfc_POLYGON2_sp
        count_strip <- count_strip + 1
        #
        # update
        sfc_SpatialPolygons <- sfc_POLYGON1
      }
      #
      # update
      pnt_upper <- strip_stations[i, 1:2]
      #
      # append
      strip_Lhrz <- rbind(strip_Lhrz, strip_stations_Lhrzs[i])
    }
  }
  #
  # append
  list_strip_sp[[count_strip]] <- f_coords2spg((sfc_SpatialPolygons[[1]])[[1]])
  #
  # append
  strip_Lhrz <- rbind(strip_Lhrz, Lhrz_total)
  #
  # get length of strip
  for (i in count_strip:2) {
    #
    # cumulative to ...
    strip_Lhrz[i] <- strip_Lhrz[i] - strip_Lhrz[i-1]
  }
  #
  # get areas
  list_areas <- f_strips_areas(pRasterDEM, list_strip_sp, strip_nodes_for_sp)
  # get
  strip_Ahrzs <- list_areas[[1]]
  strip_Aalls <- list_areas[[2]]
  # get
  strip_Lalls <- list_areas[[3]]
  # get
  strip_drops <- list_areas[[4]]
  # get
  strip_nodes2 <- list_areas[[5]]
  #
  # get sp attribute
  strip_df <- cbind(strip_Lhrz, strip_Lalls, strip_Ahrzs, strip_Aalls, strip_drops)
  #
  # colnames
  colnames(strip_df) <- c("Lhrz", "Lall", "Ahrz", "Aall", "Drop")
  #
  # get width
  strip_df$Whrz <- strip_df$Ahrz / strip_df$Lhrz
  strip_df$Wall <- strip_df$Aall / strip_df$Lall
  #
  # switch
  strip_df <- strip_df[, c("Lhrz", "Lall", "Ahrz", "Aall", "Whrz", "Wall", "Drop")]
  #
  # return
  return(list(list_strip_sp, strip_df, strip_nodes2))
}
#
#
#
f_strips_clipping_polylines <- function(pPolygons, strip_stations, strip_stations_angles, CRS_out = NA) {
  #
  # get coords of polygons
  pPolygons_coords <- pPolygons@Polygons[[1]]@coords
  #
  # get limits of polygons
  XMin <- min(pPolygons_coords[, 1])
  XMax <- max(pPolygons_coords[, 1])
  YMin <- min(pPolygons_coords[, 2])
  YMax <- max(pPolygons_coords[, 2])
  #
  # initial
  list_Lines <- list()
  #
  # initial
  data_df <- data.frame()
  #
  # get Lines
  for (i in 1:nrow(strip_stations)) {
    #
    # get station
    line_x0 <- strip_stations[i, 1]
    line_y0 <- strip_stations[i, 2]
    #
    # get angle
    line_agl <- strip_stations_angles[i]
    #
    # initial
    line_x1 <- c()
    line_y1 <- c()
    line_x2 <- c()
    line_y2 <- c()
    #
    # check
    if (line_agl == pi/2 || line_agl == -pi/2) {
      #
      # get x
      line_x1 <- line_x0
      line_x2 <- line_x0
      #
      # get y
      line_y1 <- YMin
      line_y2 <- YMax
    }
    else {
      #
      # get x
      line_x1 <- XMin
      line_x2 <- XMax
      #
      # get y
      line_y1 <- line_y0 + (line_x1-line_x0) * tan(line_agl)
      line_y2 <- line_y0 + (line_x2-line_x0) * tan(line_agl)
    }
    #
    # get pnts, line
    pnts_line1 <- matrix(rbind(c(line_x0, line_y0), c(line_x1, line_y1)), nc = 2, byrow = F)
    #
    # get ints
    pCoords1 <- f_intersection_curves(pnts_line1, pPolygons_coords)
    #
    # get nearest extending
    pCoords1_extnd <- f_pnts_nearest_extending(c(line_x0, line_y0), pCoords1, 10)
    #
    # get pnts, line
    pnts_line2 <- matrix(rbind(c(line_x0, line_y0), c(line_x2, line_y2)), nc = 2, byrow = F)
    #
    # get ints
    pCoords2 <- f_intersection_curves(pnts_line2, pPolygons_coords)
    #
    # get nearest extending
    pCoords2_extnd <- f_pnts_nearest_extending(c(line_x0, line_y0), pCoords2, 10)
    #
    # split polygon with polyline, lwgeom::st_split
    # convert polyline to polygon, sf::st_polygonize
    # these two library do not work perfectly
    # the output polyline can be converted to polygon in ArcGIS
    # so, it is this library having problem
    #
    # make an extending for every polyline
    # do not just stop at its intersection with the polygon
    #
    # get pnts
    pCoords <- matrix(rbind(pCoords1_extnd, pCoords2_extnd), nc = 2, byrow = F)
    #
    # get line
    pLine <- sp::Line(pCoords)
    pLines <- sp::Lines(list(pLine), ID = as.character(i))
    #
    # append
    list_Lines[[i]] <- pLines
    #
    # get data
    data_df <- rbind(data_df, data.frame("ID" = i))
  }
  #
  # get SpatialLines
  pSpatialLines <- sp::SpatialLines(list_Lines)
  #
  # set names
  rownames(data_df) <- c(1:nrow(data_df))
  #
  # get spldf
  spldf <-  sp::SpatialLinesDataFrame(pSpatialLines, data = data_df)
  # plot(spldf, col = c("red"))
  #
  # get prj
  if (!is.na(CRS_out)) { sp::proj4string(spldf) <- CRS_out }
  #
  # return
  return(spldf)
}
#
#
#
f_strips_areas <- function(pRasterDEM, list_strip_sp, strip_nodes) {
  #
  # get count
  count_strip <- length(list_strip_sp)
  #
  # initial
  column_Ahrzs <- data.frame()
  #
  # get area 2d
  for (i in 1:count_strip) {
    #
    # get sp
    strip_sp <- list_strip_sp[[i]]
    #
    # get area
    area2d <- strip_sp@polygons[[1]]@area
    #
    # append
    column_Ahrzs <- rbind(column_Ahrzs, area2d)
  }
  #
  # colnames
  colnames(column_Ahrzs) <- c("Ahrz")
  #
  # get count
  count_node <- nrow(strip_nodes)
  #
  # initial
  list_strip_nodes <- list()
  # initial
  idx_node_start <- 1
  idx_node_end <- 1
  #
  # get nodes for strips
  for (i in 1:count_node) {
    #
    # update
    idx_node_end <- i
    #
    # check
    if (strip_nodes[i, "type"] == "Station: Strip boundary") {
      #
      # append
      list_strip_nodes[[length(list_strip_nodes)+1]] <- strip_nodes[idx_node_start:idx_node_end, ]
      #
      # update
      idx_node_start <- i
    }
  }
  #
  # append, set an extra node for the last strip
  list_strip_nodes[[length(list_strip_nodes)+1]] <- strip_nodes[c(idx_node_start:idx_node_end, idx_node_end), ]
  #
  # initial
  column_Astps <- data.frame()
  column_drops <- data.frame()
  #
  # initial, has additional attributes
  strip_nodes2 <- c()
  #
  # save triangles
  list_triangle_sp <- list()
  #
  # get area 3d and Cz in strip
  for (i in 1:count_strip) {
    #
    # get sp
    strip_sp <- list_strip_sp[[i]]
    #
    # get coords
    coords <- strip_sp@polygons[[1]]@Polygons[[1]]@coords
    #
    # get xy
    Cx <- coords[, 1]
    Cy <- coords[, 2]
    #
    # get differential
    count_coord <- nrow(coords)
    CDx <- Cx[2:count_coord] - Cx[1:(count_coord-1)]
    CDy <- Cy[2:count_coord] - Cy[1:(count_coord-1)]
    CD <- abs(CDx) + abs(CDy)
    #
    # remove duplicate, for robusting polygon triangulation?
    if (length(which(CD==0)) >= 1) {
      #
      # remove
      Cx <- Cx[-(which(CD==0)+1)]
      Cy <- Cy[-(which(CD==0)+1)]
    }
    #
    # get z
    Cz <- raster::extract(pRasterDEM, coords)
    #
    # append
    column_drops <- rbind(column_drops, max(Cz)-min(Cz))
    #
    # polygon triangulation
    # debugging, try catch
    if (F) {
      #
      # initial
      out <- NA
      #
      # try catch
      out <- tryCatch (
        #
        ########################################################
        # Try part: define the expression(s) you want to "try" #
        ########################################################
        #
        {
          #
          # Just to highlight: 
          # If you want to use more than one R expression in the "try part" 
          # then you'll have to use curly brackets. 
          # Otherwise, just write the single expression you want to try and 
          #
          message("This is the 'try' part")
          # tragls <- rgl::triangulate(x = Cx, y = Cy, z = NULL, random = FALSE, plot = FALSE, partial = NA)
          tragls <- decido::earcut(cbind(Cx, Cy))
        },
        #
        ########################################################################
        # Condition handler part: define how you want conditions to be handled #
        ########################################################################
        #
        # Handler when a warning occurs:
        warning = function(cond) {
          #
          # message
          message("Get warning:")
          message("Here's the original warning message:")
          message(cond)
          message("\n")
          #
          # Choose a return value when such a type of condition occurs
          return(NULL)
        },
        #
        # Handler when an error occurs:
        error = function(cond) {
          #
          # message
          message("Get erros:")
          message("Here's the original error message:")
          message(cond)
          message("\n")
          #
          # Choose a return value when such a type of condition occurs
          return(NULL)
        },
        #
        ###############################################
        # Final part: define what should happen AFTER #
        # everything has been tried and/or handled    #
        ###############################################
        #
        finally = {
          #
          # if succeed
          message("Process succeed")
          message("\n")
        }
      )
      #
      # pause
      if (is.null(out)) {
        #
        # save
        f_save_list2sp(list(strip_sp), NULL, NA, "strip_sp")
        #
        # save
        coords <- data.frame(coords)
        colnames(coords) <- c("Cx", "Cy")
        f_pnts_save_points(coords, NA, "strip_coords")
        #
        # pause at here
        return(NA)
      }
    }
    # #
    # # polygon triangulation
    # # this library can not handle some situations
    # tragls <- rgl::triangulate(x = Cx, y = Cy, z = NULL, random = FALSE, plot = FALSE, partial = NA)
    #
    # polygon triangulation
    tragls <- decido::earcut(cbind(Cx, Cy))
    tragls <- matrix(tragls, 3)
    #
    # save triangles, set a switch
    if (F) {
      #
      # append
      for (k in 1:ncol(tragls)) {
        #
        # get index
        idx1 <- tragls[1, k]
        idx2 <- tragls[2, k]
        idx3 <- tragls[3, k]
        #
        # get coords
        coords <- c(Cx[idx1], Cy[idx1])
        coords <- rbind(coords, c(Cx[idx2], Cy[idx2]))
        coords <- rbind(coords, c(Cx[idx3], Cy[idx3]))
        coords <- rbind(coords, c(Cx[idx1], Cy[idx1]))
        #
        # append
        list_triangle_sp[[length(list_triangle_sp)+1]] <- f_coords2spg(coords)
      }
    }
    #
    # initial
    area3d <- 0
    #
    # get area
    for (k in 1:ncol(tragls)) {
      #
      # get index
      idx1 <- tragls[1, k]
      idx2 <- tragls[2, k]
      idx3 <- tragls[3, k]
      #
      # get edge
      a <- sqrt((Cx[idx2]-Cx[idx1])^2 + (Cy[idx2]-Cy[idx1])^2 + (Cz[idx2]-Cz[idx1])^2)
      b <- sqrt((Cx[idx3]-Cx[idx2])^2 + (Cy[idx3]-Cy[idx2])^2 + (Cz[idx3]-Cz[idx2])^2)
      c <- sqrt((Cx[idx1]-Cx[idx3])^2 + (Cy[idx1]-Cy[idx3])^2 + (Cz[idx1]-Cz[idx3])^2)
      #
      # update
      s <- (a+b+c)/2
      area <- sqrt(s*(s-a)*(s-b)*(s-c))
      area3d <- area3d + area
    }
    #
    # append
    column_Astps <- rbind(column_Astps, area3d)
    #
    # debugging
    if (F) {
      #
      # initial
      area2d <- 0
      #
      # get area
      for (k in 1:ncol(tragls)) {
        #
        # get index
        idx1 <- tragls[1, k]
        idx2 <- tragls[2, k]
        idx3 <- tragls[3, k]
        #
        # get edge
        a <- sqrt((Cx[idx2]-Cx[idx1])^2 + (Cy[idx2]-Cy[idx1])^2)
        b <- sqrt((Cx[idx3]-Cx[idx2])^2 + (Cy[idx3]-Cy[idx2])^2)
        c <- sqrt((Cx[idx1]-Cx[idx3])^2 + (Cy[idx1]-Cy[idx3])^2)
        #
        # update
        s <- (a+b+c)/2
        area <- sqrt(s*(s-a)*(s-b)*(s-c))
        area2d <- area2d + area
      }
      #
      # check
      if (abs(area2d - column_Ahrzs[i, 1]) >= 10^-6) {
        #
        # warning
        warning("Error: horizontal areas are inconsistent.")
        #
        # pause at here
        return(NULL)
      }
    }
    #
    # initial
    strip_nodes2_i <- list_strip_nodes[[i]]
    strip_nodes2_i[, "Lstp"] <- c(0)
    strip_nodes2_i[, "Czs"] <- c(NA)
    #
    # get Cz surface
    for (j in 1:(nrow(strip_nodes2_i)-1)) {
      #
      # get node
      node <- strip_nodes2_i[j, ]
      #
      # node xy extnd
      node_x_extnd <- node[1, "Cx"]
      node_y_extnd <- node[1, "Cy"]
      #
      # if Station or Initial group anchor
      if (node[1, "type"] == "Station: Strip boundary" ||
          node[1, "type"] == "Anchor: Initial group anchor") {
        #
        # get
        node_x_extnd <- node_x_extnd + signif(strip_nodes2_i[j+1, "Cx"] - node_x_extnd) * 10^-6
        node_y_extnd <- node_y_extnd + signif(strip_nodes2_i[j+1, "Cy"] - node_y_extnd) * 10^-6
      }
      #
      # if Distal group anchor
      else if (node[1, "type"] == "Anchor: Distal group anchor") {
        #
        # get
        node_x_extnd <- node_x_extnd + signif(strip_nodes2_i[j-1, "Cx"] - node_x_extnd) * 10^-6
        node_y_extnd <- node_y_extnd + signif(strip_nodes2_i[j-1, "Cy"] - node_y_extnd) * 10^-6
      }
      #
      # if identical to the first station
      else if (node[1, "Cx"] == strip_nodes2_i[1, "Cx"] &&
          node[1, "Cy"] == strip_nodes2_i[1, "Cy"]) {
        #
        # get
        node_x_extnd <- node_x_extnd + signif(strip_nodes2_i[1+2, "Cx"] - node_x_extnd) * 10^-6
        node_y_extnd <- node_y_extnd + signif(strip_nodes2_i[1+2, "Cy"] - node_y_extnd) * 10^-6
      }
      #
      # if identical to the second station
      else if (node[1, "Cx"] == strip_nodes2_i[nrow(strip_nodes2_i), "Cx"] &&
          node[1, "Cy"] == strip_nodes2_i[nrow(strip_nodes2_i), "Cy"]) {
        #
        # get
        node_x_extnd <- node_x_extnd + signif(strip_nodes2_i[nrow(strip_nodes2_i)-2, "Cx"] - node_x_extnd) * 10^-6
        node_y_extnd <- node_y_extnd + signif(strip_nodes2_i[nrow(strip_nodes2_i)-2, "Cy"] - node_y_extnd) * 10^-6
      }
      #
      # get triangle
      for (k in 1:ncol(tragls)) {
        #
        # get index
        idx1 <- tragls[1, k]
        idx2 <- tragls[2, k]
        idx3 <- tragls[3, k]
        #
        # get xyz
        tx <- rbind(Cx[idx1], Cx[idx2], Cx[idx3], Cx[idx1])
        ty <- rbind(Cy[idx1], Cy[idx2], Cy[idx3], Cy[idx1])
        tz <- rbind(Cz[idx1], Cz[idx2], Cz[idx3], Cz[idx1])
        #
        # point.in.polygon
        pip <- sp::point.in.polygon(node_x_extnd, node_y_extnd, tx, ty, mode.checked = FALSE)
        #
        # check
        if (pip != 0) {
          #
          # get normal
          nx <- (Cy[idx2]-Cy[idx1])*(Cz[idx3]-Cz[idx1]) - (Cz[idx2]-Cz[idx1])*(Cy[idx3]-Cy[idx1])
          ny <- (Cz[idx2]-Cz[idx1])*(Cx[idx3]-Cx[idx1]) - (Cx[idx2]-Cx[idx1])*(Cz[idx3]-Cz[idx1])
          nz <- (Cx[idx2]-Cx[idx1])*(Cy[idx3]-Cy[idx1]) - (Cy[idx2]-Cy[idx1])*(Cx[idx3]-Cx[idx1])
          #
          # get Czs, nz won't be zero?
          Czs <- -(nx*(node[1, "Cx"]-Cx[idx1]) + ny*(node[1, "Cy"]-Cy[idx1])) / nz + Cz[idx1]
          #
          # update
          strip_nodes2_i[j, "Czs"] <- Czs
          #
          # break
          break
        }
      }
      #
      # check
      if (is.na(strip_nodes2_i[j, "Czs"])) {
        #
        # save
        f_save_list2sp(list_strip_sp, NULL, NA, "strips_sp")
        #
        # save
        f_pnts_save_points(strip_nodes, NA, "strip_nodes")
        #
        # save
        f_save_list2sp(list(strip_sp), NULL, NA, "strip_sp")
        #
        # save
        f_pnts_save_points(strip_nodes2_i, NA, "strip_nodes2_i")
        #
        # save
        f_pnts_save_points(node, NA, "strip_node")
        #
        # save
        f_pnts_save_points(data.frame(Cx = node_x_extnd, Cy = node_y_extnd), NA, "strip_node_extnd")
        #
        # warning
        warning("Error: Czs (Cz surface) is zero.")
        #
        # pause at here
        return(NULL)
      }
    }
    #
    # append
    strip_nodes2 <- rbind(strip_nodes2, strip_nodes2_i[1:(nrow(strip_nodes2_i)-1), ])
  }
  #
  # save triangles, set a switch
  if (F) {
    #
    # save
    f_save_list2sp(list_triangle_sp, NULL, NA, "triangles")
  }
  #
  # colnames
  colnames(column_Astps) <- c("Astp")
  colnames(column_drops) <- c("Drop")
  #
  # initial
  column_Lstps <- data.frame()
  #
  # get length strip
  for (i in 2:count_node) {
    #
    # get
    strip_nodes2[i, "Lstp"] <- strip_nodes2[i-1, "Lstp"] + 
      sqrt((strip_nodes2[i, "Cx"]-strip_nodes2[i-1, "Cx"])^2 +
           (strip_nodes2[i, "Cy"]-strip_nodes2[i-1, "Cy"])^2 +
           (strip_nodes2[i, "Czs"]-strip_nodes2[i-1, "Czs"])^2)
    #
    # check
    if (strip_nodes2[i, "type"] == "Station: Strip boundary") {
      #
      # append
      column_Lstps <- rbind(column_Lstps, strip_nodes2[i, "Lstp"])
    }
  }
  #
  # switch
  strip_nodes2 <- strip_nodes2[, c("Cx", "Cy", "Czd", "Czp", "Czs", "type", "grp", "station", "Lhrz", "Lalld", "Lallp", "Lstp")]
  #
  # append
  column_Lstps <- rbind(column_Lstps, strip_nodes2[count_node, "Lstp"])
  #
  # get length strip
  for (i in count_strip:2) {
    #
    # cumulative to ...
    column_Lstps[i, 1] <- column_Lstps[i, 1] - column_Lstps[i-1, 1]
  }
  #
  # colnames, change Lstp to Lalls
  colnames(strip_nodes2) <- c("Cx", "Cy", "Czd", "Czp", "Czs", "type", "grp", "station", "Lhrz", "Lalld", "Lallp", "Lalls")
  #
  # return 
  return(list(column_Ahrzs, column_Astps, column_Lstps, column_drops, strip_nodes2))
}
#
#
#
f_pnts_grp_bnd_extending <- function(pPolygons, GrpBnd) {
  #
  # get coords of polygons
  pPolygons_coords <- pPolygons@Polygons[[1]]@coords
  #
  # get limits of polygons
  XMin <- min(pPolygons_coords[, 1])
  XMax <- max(pPolygons_coords[, 1])
  YMin <- min(pPolygons_coords[, 2])
  YMax <- max(pPolygons_coords[, 2])
  #
  # get
  GrpBnd_x1 <- GrpBnd[1, 1]
  GrpBnd_y1 <- GrpBnd[1, 2]
  GrpBnd_x2 <- GrpBnd[2, 1]
  GrpBnd_y2 <- GrpBnd[2, 2]
  GrpBnd_x0 <- (GrpBnd_x1 + GrpBnd_x2) / 2
  GrpBnd_y0 <- (GrpBnd_y1 + GrpBnd_y2) / 2
  #
  # initial
  line_x1 <- c()
  line_y1 <- c()
  line_x2 <- c()
  line_y2 <- c()
  #
  # check
  if (GrpBnd_x1 == GrpBnd_x2) {
    #
    # get x
    line_x1 <- GrpBnd_x1
    line_x2 <- GrpBnd_x1
    #
    # get y
    line_y1 <- YMin
    line_y2 <- YMax
  }
  else {
    #
    # get x
    line_x1 <- XMin
    line_x2 <- XMax
    #
    # get y
    line_y1 <- GrpBnd_y1 + (line_x1 - GrpBnd_x1) * (GrpBnd_y2 - GrpBnd_y1) / (GrpBnd_x2 - GrpBnd_x1)
    line_y2 <- GrpBnd_y1 + (line_x2 - GrpBnd_x1) * (GrpBnd_y2 - GrpBnd_y1) / (GrpBnd_x2 - GrpBnd_x1)
  }
  #
  # get pnts, line
  pnts_line1 <- matrix(rbind(c(GrpBnd_x0, GrpBnd_y0), c(line_x1, line_y1)), nc = 2, byrow = F)
  #
  # get ints
  ints1 <- f_intersection_curves(pnts_line1, pPolygons_coords)
  #
  # get nearest
  ints1_nearest <- f_pnts_nearest_extending(c(GrpBnd_x0, GrpBnd_y0), ints1, Inf)
  #
  # get pnts, line
  pnts_line2 <- matrix(rbind(c(GrpBnd_x0, GrpBnd_y0), c(line_x2, line_y2)), nc = 2, byrow = F)
  #
  # get ints
  ints2 <- f_intersection_curves(pnts_line2, pPolygons_coords)
  #
  # get nearest
  ints2_nearest <- f_pnts_nearest_extending(c(GrpBnd_x0, GrpBnd_y0), ints2, Inf)
  #
  # get extedning
  GrpBnd_extnd <- rbind(ints1_nearest, ints2_nearest)
  colnames(GrpBnd_extnd) <- c("Cx", "Cy")
  #
  # return
  return(GrpBnd_extnd)
}
#
#
#
f_pnts_nearest_extending <- function(pnt0, pnts, extnd = Inf) {
  #
  # get count
  count_pnt <- nrow(pnts)
  #
  # check
  if (count_pnt == 1) {
    #
    # get pnt
    pnt_extending_x <- pnts[1, 1] + (pnts[1, 1] - pnt0[1])/extnd
    pnt_extending_y <- pnts[1, 2] + (pnts[1, 2] - pnt0[2])/extnd
    #
    # return
    return(c(pnt_extending_x, pnt_extending_y)) 
  }
  #
  # initial
  column_dist <- c()
  #
  # get dist
  for (i in 1:count_pnt) {
    #
    # get dist
    dist <- sqrt((pnts[i, 1]- pnt0[1])^2 + (pnts[i, 2] - pnt0[2])^2)
    #
    # append
    column_dist <- rbind(column_dist, dist)
  }
  #
  # order
  ID_order <- order(column_dist)
  # get ID
  ID_nearest <- ID_order[1]
  ID_nearest2 <- ID_order[2]
  #
  # get pnt
  pnt_extending_x <- pnts[ID_nearest, 1] + (pnts[ID_nearest2, 1] - pnts[ID_nearest, 1])/extnd
  pnt_extending_y <- pnts[ID_nearest, 2] + (pnts[ID_nearest2, 2] - pnts[ID_nearest, 2])/extnd
  #
  # return
  return(c(pnt_extending_x, pnt_extending_y))
}
#
#
#
f_coords2spg <- function(coords) {
  #
  # get SpatialPolygons
  pPolygon <- sp::Polygon(coords)
  pPolygons <- sp::Polygons(list(pPolygon), ID = "1")
  pSpatialPolygons <- sp::SpatialPolygons(list(pPolygons))
  #
  # return
  return(pSpatialPolygons)
}
#
#
#
f_coords2spl <- function(coords) {
  #
  # get SpatialLines
  pLine <- sp::Line(coords)
  pLines <- sp::Lines(list(pLine), ID = "1")
  pSpatialLines <- sp::SpatialLines(list(pLines))
  #
  # return
  return(pSpatialLines)
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
    # get coords
    int.coords <- int.coords.line
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
    # get coords
    int.coords <- rbind(int.coords, int.coords.line)
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
  # colnames
  colnames(int.coords) <- c("Cx", "Cy")
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
    return(NA)
  }
  #
  # return
  return(CoordZ)
}
#
#
#
f_pnts_distances <- function(column_anchors) {
  #
  # get count
  count_anchor <- nrow(column_anchors)
  #
  # if only one anchor, return Inf
  if (count_anchor == 1) {
    #
    # return
    return(Inf)
  }
  #
  # initial
  column_distances <- c() 
  #
  # get distances
  for (i in 1:(count_anchor-1)) {
    #
    # get distance
    distance <- 0
    distance <- distance + (column_anchors[i+1, 1] - column_anchors[i, 1])^2
    distance <- distance + (column_anchors[i+1, 2] - column_anchors[i, 2])^2
    distance <- distance + (column_anchors[i+1, 3] - column_anchors[i, 3])^2
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
f_pnts_save_polylines <- function(list_pnts, CRS_out, file_out) {
  #
  # initial
  list_spldf <- list()
  #
  # get spldf
  for (i in 1:length(list_pnts)) {
    #
    # get pnts
    pnts <- list_pnts[[i]]
    # #
    # # get data_df
    # data_df <- 
    #
    # get spldf
    list_spldf[[i]] <- f_pnts2shp_polyline(pnts, CRS_out, data_df = NULL)
  }
  #
  # save
  f_save_list2spldf(list_spldf, CRS_out, file_out)
}
#
#
#
f_pnts_save_proflle <- function(pnts_prfl, paras_strip, paras_input, CRS_out, file_out, LineID = 1) {
  #
  # get parameters for landslide profile
  paras <- f_pnts2paras_profile(pnts_prfl, paras_strip)
  #
  # get data.frame
  paras_input <- data.frame(iMGPC = paras_input[1], iMGAD = paras_input[2], iMEAR = paras_input[3], iMEDA = paras_input[4], iMSHL = paras_input[5])
  #
  # get data_df, inputs leed paras
  data_df <- cbind(paras_input, data.frame(ID = LineID), paras)
  #
  # get spldf
  spldf <- f_pnts2shp_polyline(pnts_prfl, CRS_out, data_df)
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
  if (!is.na(CRS_out)) { sp::proj4string(PntsCoords_sptdf) = CRS_out }
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
    # colnames
    colnames(data_df) <- c("c")
    #
    # rownames
    rownames(data_df) <- c("a")
  }
  #
  # get spldf
  spldf <- f_data2spldf(data_spl, data_df)
  #
  # get prj
  if (!is.na(CRS_out)) { sp::proj4string(spldf) = CRS_out }
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
  # get counts
  count_col <- (x_max - x_min) / raster_cellsize
  count_row <- (y_max - y_min) / raster_cellsize
  #
  # get counts
  raster::ncol(pnts_raster0) <- count_col
  raster::nrow(pnts_raster0) <- count_row
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
f_pnts2paras_profile <- function(pnts, paras_strip) {
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
    x1 = pnts[i-1, "Cx"]
    y1 = pnts[i-1, "Cy"]
    z1 = pnts[i-1, "Czs"]
    #
    x2 = pnts[i, "Cx"]
    y2 = pnts[i, "Cy"]
    z2 = pnts[i, "Czs"]
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
  d_overall <- sqrt((pnts[count_pnts, "Cx"] - pnts[1, "Cx"])^2 + 
                    (pnts[count_pnts, "Cy"] - pnts[1, "Cy"])^2 +
                    (pnts[count_pnts, "Czs"] - pnts[1, "Czs"])^2)
  T_overall <- s_profile / d_overall
  #
  # get tortuosity horizontal
  d_horizontal <- sqrt((pnts[count_pnts, "Cx"] - pnts[1, "Cx"])^2 + 
                       (pnts[count_pnts, "Cy"] - pnts[1, "Cy"])^2)
  T_horizontal <- l_profile / d_horizontal
  #
  # get tortuosity longitudinal
  d_longitudinal <- sqrt((l_profile)^2 +
                         (h_profile)^2)
  T_longitudinal <- s_profile / d_longitudinal
  #
  # initial
  paras <- c()
  #
  # get paras
  if (is.null(paras_strip)) {
    #
    # width no
    paras <- data.frame(Dall = d_overall, Lall = s_profile,
                        Dhrz = d_horizontal, Lhrz = l_profile,
                        Dlng = d_longitudinal,
                        Hfnl = -h_profile,
                        Hmax = -min(h_series), Hint = H_integral,
                        Gmax = max(g_series), Gint = G_integral,
                        Tall = T_overall, Thrz = T_horizontal, Tlng = T_longitudinal,
                        miuapp = miuapp_profile)
  }
  else {
    #
    # width yes
    paras <- data.frame(Dall = d_overall, Lall = s_profile,
                        Dhrz = d_horizontal, Lhrz = l_profile,
                        Dlng = d_longitudinal,
                        Hfnl = -h_profile,
                        Hmax = -min(h_series), Hint = H_integral,
                        Gmax = max(g_series), Gint = G_integral,
                        Tall = T_overall, Thrz = T_horizontal, Tlng = T_longitudinal,
                        miuapp = miuapp_profile, 
                        # Lall2 = paras_strip[3], 
                        Aall = paras_strip[4], Wall = paras_strip[4] / paras_strip[3], epsall = paras_strip[3]^2 / paras_strip[4],
                        # Lhrz2 = paras_strip[1], 
                        Ahrz = paras_strip[2], Whrz = paras_strip[2] / paras_strip[1], epshrz = paras_strip[1]^2 / paras_strip[2])
  }
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
f_save_list2sp <- function(list_sp, data_df, CRS_out, file_out) {
  #
  # initial
  spdf <- NULL
  #
  # get class
  class_sp <- class(list_sp[[1]])[1]
  #
  # SpatialPolygons
  if (class_sp == "SpatialPolygons") {
    #
    # initial
    list_spdf <- list()
    #
    # get spdf
    for (i in 1:length(list_sp)) { 
      #
      # get data
      data_df_i <- data.frame(i)
      colnames(data_df_i) <- "ID"
      rownames(data_df_i) <- i
      #
      # append
      if (!is.null(data_df)) { data_df_i <- cbind(data_df_i, data_df[i, ]) }
      #
      # append
      list_spdf[[i]] <- sp::SpatialPolygonsDataFrame(list_sp[[i]], data = data_df_i, match.ID = F)
    }
    #
    # initial
    spdf <- list_spdf[[1]]
    #
    # append spdf
    if (length(list_spdf) >= 2) {
      #
      # append
      for (i in 2:length(list_spdf)) {
        #
        # append
        spdf <- rbind(spdf, list_spdf[[i]], makeUniqueIDs = TRUE)
      }
    }
  }
  #
  # SpatialLines
  else if (class_sp == "SpatialLines") {
    #
    # initial
    list_Lines <- list()
    #
    # initial
    data_df_ID <- data.frame()
    #
    # get stack
    for (i in 1:length(list_sp)) {
      #
      # get SpatialLines
      pSpatialLines <- list_sp[[i]]
      #
      # get Lines
      pLines <- pSpatialLines@lines
      #
      # change ID
      pLines[[1]]@ID <- toString(i)
      #
      # update list
      list_Lines[[i]] <- pLines[[1]]
      #
      # get data
      data_df_i <- data.frame(i)
      colnames(data_df_i) <- "ID"
      rownames(data_df_i) <- i
      #
      # get data
      data_df_ID <- rbind(data_df_ID, data_df_i)
    }
    #
    # get SpatialLines
    pSpatialLines = sp::SpatialLines(list_Lines)
    #
    # get data_df
    if (!is.null(data_df)) { data_df <- cbind(data_df_ID, data_df) }
    else { data_df <- data_df_ID }
    #
    # get spdf
    spdf <-  sp::SpatialLinesDataFrame(pSpatialLines, data = data_df)
  }
  #
  # check
  else {
    #
    # Warning
    warning("Warning: not right type of sp (spatial).")
  }
  #
  # get prj
  if (!is.na(CRS_out)) { sp::proj4string(spdf) <- CRS_out }
  #
  # write
  rgdal::writeOGR(spdf, dsn = ".", layer = file_out, overwrite = T, driver="ESRI Shapefile")
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
  list_Lines <- list()
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
    list_Lines[[i]] <- pLines
    #
    # get data
    data_df <- rbind(data_df, spldf_i@data)
  }
  #
  # get SpatialLines
  pSpatialLines = sp::SpatialLines(list_Lines)
  #
  # set names
  rownames(data_df) <- c(1:nrow(data_df))
  #
  # get spldf
  spldf <-  sp::SpatialLinesDataFrame(pSpatialLines, data = data_df)
  # plot(spldf, col = c("red"))
  #
  # get prj
  if (!is.na(CRS_out)) { sp::proj4string(spldf) <- CRS_out }
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