#
# ###################################################################### #
# ###################################################################### #
#                                                                        #
#                                ALPA                                    #
#                                                                        #
#      An R-Script for for automatic analysis of landslide profile       #
#                                                                        #
#                             Version 2.7                                #
#                                                                        #
#                          November 21th, 2022                           #
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
ALPA <- function(FileDEM, FileLandslides, MinGrpPntCnt = 3, MinGrpAcrDst = 0, MinEndAptRto = 0, MaxEndDvtAgl = 180, EndHdlAgrCfg = 1, MinStpHrzLen = 30, FilePrefix = "", OutputTemp = F) {
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
  file_out_profile <- paste(FilePrefix, "_lsd_profile2d", sep = "", collapse = NULL)
  file_out_xlsx <- paste(FilePrefix, "_lsd_profile2d.xlsx", sep = "", collapse = NULL)
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
    spldf_profile <- f_lasld_split(raster_dem, pPolygons, MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, EndHdlAgrCfg, MinStpHrzLen, FilePrefix, i, OutputTemp)
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
f_lasld_split <- function(pRasterDEM, pPolygons, MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, EndHdlAgrCfg, MinStpHrzLen, FilePrefix, FileID, OutputTemp) {
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
  file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_lsd_pln2d", sep = "", collapse = NULL)
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
    file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts", sep = "", collapse = NULL)
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
    file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts", sep = "", collapse = NULL)
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
    file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grp.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts_split, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grp", sep = "", collapse = NULL)
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
    # get finalization
    Finalization <- min(data.frame(list_grpm)) == 1
    #
    # check
    if (Finalization) {
      #
      # get path (centers and bnd sides)
      if (EndHdlAgrCfg <= 6) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, EndHdlAgrCfg, T)
      }
      else if (EndHdlAgrCfg == 7) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 3, T)
        #
        # if fail, do not change algorithm for generating anchors in finalization
        if (is.null(list_pnts_path)) { list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 1, T) }
      }
      else if (EndHdlAgrCfg == 8) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 5, T)
        #
        # if fail, do not change algorithm for generating anchors in finalization
        if (is.null(list_pnts_path)) { list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 1, T) }
      }
      else if (EndHdlAgrCfg == 9) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 5, T)
        #
        # if fail, do not change algorithm for generating anchors in finalization
        if (is.null(list_pnts_path)) { list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 3, T) }
      }
      else if (EndHdlAgrCfg == 10) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 4, T)
        #
        # if fail, do not change algorithm for generating anchors in finalization
        if (is.null(list_pnts_path)) { list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 2, T) }
      }
      else if (EndHdlAgrCfg == 11) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 6, T)
        #
        # if fail, do not change algorithm for generating anchors in finalization
        if (is.null(list_pnts_path)) { list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 2, T) }
      }
      else if (EndHdlAgrCfg == 12) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 6, T)
        #
        # if fail, do not change algorithm for generating anchors in finalization
        if (is.null(list_pnts_path)) { list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 4, T) }
      }
    }
    else {
      #
      # get path (centers and bnd sides)
      if (EndHdlAgrCfg <= 6) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, EndHdlAgrCfg, F)
      }
      else if (EndHdlAgrCfg == 7 || EndHdlAgrCfg == 8) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 1, F)
      }
      else if (EndHdlAgrCfg == 9) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 3, F)
      }
      else if (EndHdlAgrCfg == 10 || EndHdlAgrCfg == 11) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 2, F)
      }
      else if (EndHdlAgrCfg == 12) {
        #
        # get
        list_pnts_path <- f_pnts_path(list_pnts, IDB_max, pPolygons, 4, F)
      }
    }
    #
    # get centers
    pnts_split_anchors <- list_pnts_path[[1]]
    #
    # if save
    if (OutputTemp || Finalization) {
      #
      # get length
      length_prfl <- sum(f_pnts_distances(pnts_split_anchors[, c("Cx", "Cy")]))
      #
      # initial
      count_strip <- 0
      #
      # get count
      while (T) {
        #
        # update
        if (length_prfl/(count_strip+1) >= MinStpHrzLen) { count_strip <- count_strip+1 }
        else { break }
      }
      #
      # check, at least two strip, do not allow no station
      count_strip <- max(count_strip, 2)
      #
      # get strips, strip might be merged
      list_strips <- f_strips(pPolygons, list_pnts_path, count_strip)
      #
      # get
      list_strip_sp <- list_strips[[1]]
      strip_nodes_for_sp <- list_strips[[2]]
      #
      # get areas
      list_areas <- f_strips_areas(pRasterDEM, list_strip_sp, strip_nodes_for_sp)
      #
      # get 
      strip_sp <- list_strip_sp
      # get 
      strip_df <- list_areas[[1]]
      # get 
      strip_nodes <- list_areas[[2]]
      # get 
      strip_triangle_sp <- list_areas[[3]]
      strip_triangle_df <- list_areas[[4]]
      #
      # get paras
      paras_strip <- c(sum(strip_df[, "Lhrz"]), sum(strip_df[, "Ahrz"]), sum(strip_df[, "Lall"]), sum(strip_df[, "Aall"]))
    }
    #
    # save, for every split
    if (OutputTemp) {
      #
      # save pnts
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), ".xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split, raster_CRS, file_out)
      #
      # save anchors
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_anchors.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split_anchors, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_anchors", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split_anchors, raster_CRS, file_out)
      #
      # save strips
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_df, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips", sep = "", collapse = NULL)
      f_save_list2sp(strip_sp, strip_df, raster_CRS, file_out)
      #
      # save strips
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_triangles.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_triangle_df, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_triangles", sep = "", collapse = NULL)
      f_save_list2sp(strip_triangle_sp, strip_triangle_df, raster_CRS, file_out)
      #
      # save nodes
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_nodes.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_nodes, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_nodes", sep = "", collapse = NULL)
      f_pnts_save_points(strip_nodes, raster_CRS, file_out)
      #
      # save shp pln
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_nodes_pln2d", sep = "", collapse = NULL)
      spldf_profile <- f_pnts_save_proflle(strip_nodes, paras_strip, c(MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, EndHdlAgrCfg, MinStpHrzLen), raster_CRS, file_out, FileID)
      #
      # save paras
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sprintf("%03d", pnts_split[nrow(pnts_split), "grp"]), "_strips_nodes_pln2d.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(spldf_profile@data, file_out)
    }
    #
    # save, and break, if no group need split
    if (Finalization) {
      # 
      # set SMS, to be "SMG", stop because all subgroup stop split
      pnts_split[, "SMS"] <- "SMG"
      #
      # save pnts
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split, raster_CRS, file_out)
      #
      # save anchors
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_anchors.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split_anchors, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_anchors", sep = "", collapse = NULL)
      f_pnts_save_points(pnts_split_anchors, raster_CRS, file_out)
      #
      # save strips
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_strips.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_df, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_strips", sep = "", collapse = NULL)
      f_save_list2sp(strip_sp, strip_df, raster_CRS, file_out)
      #
      # save strips
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_strips_triangles.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_triangle_df, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_strips_triangles", sep = "", collapse = NULL)
      f_save_list2sp(strip_triangle_sp, strip_triangle_df, raster_CRS, file_out)
      #
      # save nodes
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_strips_nodes.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(strip_nodes, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_strips_nodes", sep = "", collapse = NULL)
      f_pnts_save_points(strip_nodes, raster_CRS, file_out)
      #
      # save shp pln
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_strips_nodes_pln2d", sep = "", collapse = NULL)
      spldf_profile <- f_pnts_save_proflle(strip_nodes, paras_strip, c(MinGrpPntCnt, MinGrpAcrDst, MinEndAptRto, MaxEndDvtAgl, EndHdlAgrCfg, MinStpHrzLen), raster_CRS, file_out, FileID)
      #
      # save paras
      file_out <- paste(FilePrefix, "_lsd", sprintf("%06d", FileID), "_pnts_grps_strips_nodes_pln2d.xlsx", sep = "", collapse = NULL)
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
          pnts_i[, "SMG"] <- "SGL"
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
          # get pnts_i_bnd
          pnts_i_upper_bnd <- pnts_i_upper[which(pnts_i_upper[, "IDB"] != 0), ]
          #
          # get sorted
          pnts_i_upper_bnd_ascending <- f_pnts_bnd_ascending(pnts_i_upper_bnd)
          #
          # get pnts0
          pnts0_i_upper <- rbind(pnts_i_upper_bnd_ascending[1, c("Cx", "Cy")], pnts_i_upper_bnd_ascending[nrow(pnts_i_upper_bnd_ascending), c("Cx", "Cy")])
          #
          # get bnd sides
          list_end_pnts_i_upper <- f_pnts_bnd_end_even(pnts_i_upper, pnts0_i_upper)
          #
          # check
          if (!is.list(list_end_pnts_i_upper) && is.na(list_end_pnts_i_upper)) {
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
          # get pnts_i_bnd
          pnts_i_lower_bnd <- pnts_i_lower[which(pnts_i_lower[, "IDB"] != 0), ]
          #
          # get sorted
          pnts_i_lower_bnd_ascending <- f_pnts_bnd_ascending(pnts_i_lower_bnd)
          #
          # get pnts0
          pnts0_i_lower <- rbind(pnts_i_lower_bnd_ascending[nrow(pnts_i_lower_bnd_ascending), c("Cx", "Cy")], pnts_i_lower_bnd_ascending[1, c("Cx", "Cy")])
          #
          # get bnd sides
          list_end_pnts_i_lower <- f_pnts_bnd_end_even(pnts_i_lower, pnts0_i_lower)
          #
          # check
          if (!is.list(list_end_pnts_i_lower) && is.na(list_end_pnts_i_lower)) {
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
        # #
        # # show step time
        # f_show_time(paste("Check end group prolong done.", "\n", sep = ""), time_step_start)
        #
        # make a copy
        list_pnts_new_for_check <- list_pnts_new
        #
        # update list_pnts_new
        list_pnts_new_for_check[[grpu]] <- pnts_i_upper
        list_pnts_new_for_check[[grpu + 1]] <- pnts_i_lower
        #
        # if still have pnts, append
        if ((i+1) <= length(list_pnts)) {
          #
          # append
          for (k in (i+1):length(list_pnts)) {
            #
            # get pnts_i
            pnts_k <- list_pnts[[k]]
            list_pnts_new_for_check[[length(list_pnts_new_for_check) + 1]] <- pnts_k
          }
        }
        #
        # get anchors and distances (before split)
        pnts_split_anchors_before_split <- pnts_split_anchors
        #
        # get path (centers and bnd sides)
        if (EndHdlAgrCfg <= 6) {
          #
          # get
          list_pnts_path <- f_pnts_path(list_pnts_new_for_check, IDB_max, pPolygons, EndHdlAgrCfg, F)
        }
        else if (EndHdlAgrCfg == 7 || EndHdlAgrCfg == 8) {
          #
          # get
          list_pnts_path <- f_pnts_path(list_pnts_new_for_check, IDB_max, pPolygons, 1, F)
        }
        else if (EndHdlAgrCfg == 9) {
          #
          # get
          list_pnts_path <- f_pnts_path(list_pnts_new_for_check, IDB_max, pPolygons, 3, F)
        }
        else if (EndHdlAgrCfg == 10 || EndHdlAgrCfg == 11) {
          #
          # get
          list_pnts_path <- f_pnts_path(list_pnts_new_for_check, IDB_max, pPolygons, 2, F)
        }
        else if (EndHdlAgrCfg == 12) {
          #
          # get
          list_pnts_path <- f_pnts_path(list_pnts_new_for_check, IDB_max, pPolygons, 4, F)
        }
        #
        # check
        if (is.null(list_pnts_path)) {
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
        # #
        # # show step time
        # f_show_time(paste("Check group corners done.", "\n", sep = ""), time_step_start)
        #
        # get centers
        pnts_split_anchors <- list_pnts_path[[1]]
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
        # #
        # # show step time
        # f_show_time(paste("Check end group deviation done.", "\n", sep = ""), time_step_start)
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
        # #
        # # show step time
        # f_show_time(paste("Check minimum distance done.", "\n", sep = ""), time_step_start)
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
  # get deviation angle
  vectors_agl <- acos(vectors_cos) * 180 / pi
  #
  # numerical error, necessary?
  if (vectors_cos > 1) { vectors_agl <- 0 }
  if (vectors_cos < -1) { vectors_agl <- 180 }
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
f_pnts_bnd_sides <- function(pnts_bnd, IDB_max) {
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
f_pnts_bnd_sides_corners_check <- function(list_bnd_sides) {
  #
  # get corners
  bnd_sides_corners <- f_pnts_bnd_sides_corners(list_bnd_sides)
  #
  # initial
  diff41 <- T
  #
  # check
  if ((bnd_sides_corners[4, "Cx"] == bnd_sides_corners[1, "Cx"]) &&
      (bnd_sides_corners[4, "Cy"] == bnd_sides_corners[1, "Cy"]) ) {
    #
    #
    diff41 <- F
  }
  #
  # initial
  diff32 <- T
  #
  # check
  if ((bnd_sides_corners[3, "Cx"] == bnd_sides_corners[2, "Cx"]) &&
      (bnd_sides_corners[3, "Cy"] == bnd_sides_corners[2, "Cy"]) ) {
    #
    #
    diff32 <- F
  }
  #
  # return
  return(c(diff41, diff32))
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
f_pnts_bnd_end_even <- function(pnts, pnts0) {
  #
  # do not use circular statistics
  # just loop all the bnd pnts, to see which bnd pnt or bnd pair splits most equally
  #
  # also, check if the most equally bnd pnt is blocked by other bnd pnt
  # straight flow path that crosses the bnd is not accepted
  #
  # it might be slightly slower than using circular statistics
  # but much more robust and concise
  #
  # situation is very complex
  # so, check unblocked for all bnd pnts
  # although might cost more time, but more robust and concise
  #
  # get pnt0
  pnt0 <- data.frame((pnts0[1, "Cx"] + pnts0[2, "Cx"])/2, (pnts0[1, "Cy"] + pnts0[2, "Cy"])/2)
  colnames(pnt0) <- c("Cx", "Cy")
  #
  # initial, remove pnt0 from pnts, might be redundant
  pnts_pnt0_rm <- c()
  #
  # initial
  pnts_angles <- c()
  #
  # get angles
  for (i in 1:nrow(pnts)) {
    #
    # get xy
    lx <- as.numeric(pnts[i, "Cx"] - pnt0[1])
    ly <- as.numeric(pnts[i, "Cy"] - pnt0[2])
    #
    # if not the same
    if (lx != 0 || ly != 0) {
      #
      # append
      pnts_pnt0_rm <- rbind(pnts_pnt0_rm, pnts[i, ])
      #
      # get angle
      pnts_angles <- c(pnts_angles, atan2(ly, lx))
    }
  }
  #
  # get bnd pnts
  pnts_pnt0_rm_bnd <- pnts_pnt0_rm[which(pnts_pnt0_rm[, "IDB"] != 0), ]
  #
  # get sorted
  pnts_pnt0_rm_bnd_ascending <- f_pnts_bnd_ascending(pnts_pnt0_rm_bnd)
  #
  # get pnts_curve, for bnd pnts
  pnts_curve_bnd <- pnts_pnt0_rm_bnd_ascending[, c("Cx", "Cy")]
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
    lx <- as.numeric(pnts_curve_bnd[i, "Cx"] - pnt0[1])
    ly <- as.numeric(pnts_curve_bnd[i, "Cy"] - pnt0[2])
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
    ints <- f_pnts_bnd_end_even_ints(pnts_curve_bnd, pnt0, pnt_bnd)
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
      # for example, pnt0 is circled by pnts
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
      pnts_curve_segment <- rbind(as.numeric(pnt0))
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
  bnd_sides_right <- pnts_pnt0_rm_bnd_ascending[1:ID_split_right, ]
  bnd_sides_left <- pnts_pnt0_rm_bnd_ascending[ID_split_left:nrow(pnts_pnt0_rm_bnd_ascending), ]
  #
  # return
  return(list(bnd_sides_right, bnd_sides_left))
}
#
#
#
f_pnts_bnd_end_even_ints <- function(pnts_curve_bnd, pnt0, pnt_bnd) {
  #
  # get pnts_curve
  pnts_curve_segment <- rbind(as.numeric(pnt0), as.numeric(pnt_bnd[1, c("Cx", "Cy")]))
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
  # landing_x <- pnt0[1]
  # landing_y <- pnt0[2]
  # #
  # # get bnd pnt, xy
  # pnt_bnd_x <- pnt_bnd[1, c("Cx")]
  # pnt_bnd_y <- pnt_bnd[1, c("Cy")]
  # #
  # # get, landing pnt
  # if (pnt_bnd_x == pnt0[1]) {
  #   #
  #   # update
  #   landing_x <- pnt0[1]
  #   #
  #   # update
  #   if (pnt_bnd_y < pnt0[2]) { landing_y <- Ymin }
  #   else { landing_y <- Ymax }
  #   #
  # }
  # else {
  #   #
  #   # update
  #   if (pnt_bnd_x < pnt0[1]) { landing_x <- Xmin }
  #   else { landing_x <- Xmax }
  #   #
  #   # update
  #   landing_y <- (pnt_bnd_y - pnt0[2]) / (pnt_bnd_x - pnt0[1]) * (landing_x - pnt0[1]) + pnt0[2]
  # }
  # #
  # # get pnts_curve
  # pnts_curve_segment <- rbind(as.numeric(pnt0), as.numeric(c(landing_x, landing_y)))
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
f_pnts_bnd_end_mbb <- function(pnts_bnd, pnts0) {
  #
  # convex hull, to reduce pnts
  pnts_chull <- pnts_bnd[chull(pnts_bnd[, c("Cx", "Cy")]), ]
  pnts_chull <- pnts_chull[, c("Cx", "Cy")]
  #
  # get mbb
  colnames(pnts_chull) <- c("x", "y")
  mbb <- shotGroups::getMinBBox(pnts_chull)
  #
  # get edges
  mbb_edges <- mbb$pts
  mbb_edges <- rbind(mbb_edges, mbb_edges[1, ])
  colnames(mbb_edges) <- c("Cx", "Cy")
  #
  # debugging
  if (F) {
    #
    # save
    f_pnts_save_polylines(list(mbb_edges), NA, "mbb")
  }
  #
  # get Agl0
  Agl0 <- atan2(pnts0[2, "Cy"]-pnts0[1, "Cy"], pnts0[2, "Cx"]-pnts0[1, "Cx"])
  #
  # initial
  column_agl <- c()
  column_agl_icld <- c()
  # 
  # get, angle of edge
  for (i in 1:4) {
    #
    # get angle
    agl <- atan2(mbb_edges[i+1, "Cy"]-mbb_edges[i, "Cy"], mbb_edges[i+1, "Cx"]-mbb_edges[i, "Cx"])
    #
    # get angle, included
    agl_icld <- f_angle_included(Agl0, agl)
    #
    # check
    if (agl_icld > pi/2) { 
      #
      # reverse
      agl <- agl - pi
      agl_icld <- pi - agl_icld
    }
    #
    # check
    if (agl > pi) { agl <- agl - pi*2 }
    if (agl < -pi) { agl <- agl + pi*2 }
    #
    # append
    column_agl <- rbind(column_agl, agl)
    column_agl_icld <- rbind(column_agl_icld, agl_icld)
  }
  #
  # initial, middle pnts
  column_mp <- c()
  # 
  # get, middle pnts
  for (i in 1:4) {
    #
    # get coords
    x <- (mbb_edges[i, "Cx"] + mbb_edges[i+1, "Cx"])/2
    y <- (mbb_edges[i, "Cy"] + mbb_edges[i+1, "Cy"])/2
    #
    # append
    column_mp <- rbind(column_mp, c(x, y))
  }
  #
  # colnames
  colnames(column_mp) <- c("Cx", "Cy")
  column_mp <- data.frame(column_mp)
  #
  # get pnt0
  pnt0 <- data.frame((pnts0[1, "Cx"] + pnts0[2, "Cx"])/2, (pnts0[1, "Cy"] + pnts0[2, "Cy"])/2)
  colnames(pnt0) <- c("Cx", "Cy")
  #
  # get coords
  x0 <- pnt0[1, "Cx"]
  y0 <- pnt0[1, "Cy"]
  #
  # initial, distance
  column_dist <- c()
  # 
  # get, distance
  for (i in 1:4) {
    #
    # get coords
    x <- column_mp[i, "Cx"]
    y <- column_mp[i, "Cy"]
    #
    # get dist
    dist <- sqrt((x - x0)^2 + (y - y0)^2)
    #
    # append
    column_dist <- rbind(column_dist, dist)
  }
  #
  # get idx
  idxs <- which(column_agl_icld == min(column_agl_icld))
  #
  # owing to numerical errors
  # possibly, two subtense have different included angles
  if (length(idxs) == 1) {
    #
    idx1 <- idxs[1]
    idx2 <- which(column_agl_icld == min(column_agl_icld[-idx1]))
  }
  # 
  # owing to numerical errors
  # possibly, two subtense have the same included angle
  else {
    #
    idx1 <- idxs[1]
    idx2 <- idxs[2]
  }
  #
  # get idx
  if (column_dist[idx1] >= column_dist[idx2]) { idx <- idx1 }
  else { idx <- idx2 }
  #
  # get agl and end_Anchor0
  GrpBndAgl_end <- column_agl[idx]
  end_Anchor0 <- data.frame(column_mp[idx, ])
  #
  # get split
  list_split <- f_pnts_group_end_split(pnts_bnd, pnt0, end_Anchor0)
  bnd_sides_right <- list_split[[1]]
  bnd_sides_left <- list_split[[2]]
  #
  # return
  return(list(bnd_sides_right, bnd_sides_left, GrpBndAgl_end))
}
#
#
#
f_pnts_bnd_end_quad <- function(pnts_bnd, pnts0) {
  #
  # convex hull, to reduce pnts
  # and, also, to guarantee quadrilateral be convex
  pnts_chull <- pnts_bnd[chull(pnts_bnd[, c("Cx", "Cy")]), ]
  pnts_chull <- pnts_chull[, c("Cx", "Cy")]
  #
  # get count
  count_pnt <- nrow(pnts_chull)
  #
  # get pnt0
  pnt0 <- data.frame((pnts0[1, "Cx"] + pnts0[2, "Cx"])/2, (pnts0[1, "Cy"] + pnts0[2, "Cy"])/2)
  colnames(pnt0) <- c("Cx", "Cy")
  #
  # get Agl0
  Agl0 <- atan2(pnts0[2, "Cy"]-pnts0[1, "Cy"], pnts0[2, "Cx"]-pnts0[1, "Cx"])
  #
  # keep original
  Agl0_original <- Agl0
  #
  # get perpendicular of Agl0
  Agl0_perp <- Agl0 + pi/2
  Agl0_perp_x <- cos(Agl0_perp)
  Agl0_perp_y <- sin(Agl0_perp)
  #
  # initial
  dots_sum <- 0
  #
  # get dots
  for (i in 1:count_pnt) {
    #
    # append
    dots_sum <- dots_sum + (pnts_chull[i, "Cx"] - pnt0[1]) * Agl0_perp_x + (pnts_chull[i, "Cy"] - pnt0[2]) * Agl0_perp_y
  }
  #
  # check
  if (dots_sum < 0 ) { Agl0 <- Agl0 + pi }
  #
  # check
  if (Agl0 > pi) { Agl0 <- Agl0 - pi*2 }
  #
  # get initials
  par_init <- c(Agl0+pi, 0.5, 0.5)
  #
  # get optim
  if (F) {
    #
    # get bounds
    bound_lower <- c(Agl0 + pi/2, 0.01, 0.01)
    bound_upper <- c(Agl0 + pi/2 + pi, 0.99, 0.99)
    #
    # optim
    res <- optim(par = par_init, f_pnts_bqa, pnts = pnts_chull, Agl0 = Agl0, pnt0 = pnt0,
                 lower = bound_lower, upper = bound_upper, method = "L-BFGS-B")
    # # optim
    # # using SANN is too inefficient
    # # and, because no bounds, get wrong results?
    # res <- optim(par = par_init, f_pnts_bqa, pnts = pnts_chull, Agl0 = Agl0, pnt0 = pnt0,
    #              method = "SANN")
  }
  else {
    #
    #
    # using TEN intervals
    # to avoid getting local minimum?
    #
    #
    # get bounds
    column_bounds <- Agl0 + pi/2 + pi*(0:10)/10
    # get bounds, count
    count_bounds <- length(column_bounds) - 1
    # get bounds
    column_bound_lower <- column_bounds[1:count_bounds]
    column_bound_upper <- column_bounds[2:(count_bounds+1)]
    #
    # initial
    list_res <- list()
    column_value <- c()
    # 
    # get optim
    for (i in 1:count_bounds) {
      #
      # get bounds
      bound_lower <- c(column_bound_lower[i], 0.01, 0.01)
      bound_upper <- c(column_bound_upper[i], 0.99, 0.99)
      #
      # optim
      res <- optim(par = par_init, f_pnts_bqa, pnts = pnts_chull, Agl0 = Agl0, pnt0 = pnt0,
                   lower = bound_lower, upper = bound_upper, method = "L-BFGS-B")
      #
      # append
      list_res[[(length(list_res)+1)]] <- res
      column_value <- rbind(column_value, res$value)
    }
    #
    # get optim
    idx <- which(column_value == min(column_value))
    res <- list_res[[idx[1]]]
  }
  #
  # get optim
  Aglx_optim <- res$par
  #
  # get Agls
  Agl2 <- Aglx_optim[1]
  # get Agls
  Agl1_ratio <- Aglx_optim[2]
  Agl3_ratio <- Aglx_optim[3]
  # get Agls
  Agl1_min <- max(Agl2-pi, Agl0)
  Agl3_min <- max(Agl0+pi*2-pi, Agl2)
  # get Agls
  Agl1_max <- min(Agl0+pi, Agl2)
  Agl3_max <- min(Agl2+pi, Agl0+pi*2)
  # get Agls
  Agl1 <- Agl1_min + (Agl1_max-Agl1_min) * Agl1_ratio
  Agl3 <- Agl3_min + (Agl3_max-Agl3_min) * Agl3_ratio
  #
  # get agls
  Agls <- c(Agl0, Agl1, Agl2, Agl3)
  #
  # get ints
  list_ints <- f_pnts_bqc(pnts_chull, Agls, pnt0)
  # get ints
  int01 <- list_ints[[1]]
  int12 <- list_ints[[2]]
  int23 <- list_ints[[3]]
  int30 <- list_ints[[4]]
  #
  # debugging
  if (F) {
    # #
    # # save
    # f_pnts_save_points(rbind(int01, int12, int23, int30), NA, "int")
    #
    # get
    quadrilateral_edges <- rbind(int01, int12, int23, int30, int01)
    colnames(quadrilateral_edges) <- c("Cx", "Cy")
    #
    # save
    f_pnts_save_polylines(list(quadrilateral_edges), NA, "quadrilateral")
  }
  #
  # get
  GrpBndAgl_end <- Agl2
  #
  # check
  agl_icld <- f_angle_included(Agl0_original, GrpBndAgl_end)
  # reverse
  if (agl_icld > pi/2) { GrpBndAgl_end <- GrpBndAgl_end - pi }
  #
  # check
  if (GrpBndAgl_end > pi) { GrpBndAgl_end <- GrpBndAgl_end - pi*2 }
  if (GrpBndAgl_end < -pi) { GrpBndAgl_end <- GrpBndAgl_end + pi*2 }
  #
  # get end anchor, only direction
  end_Anchor0 <- data.frame((int12[1, "Cx"] + int23[1, "Cx"])/2, (int12[1, "Cy"] + int23[1, "Cy"])/2)
  colnames(end_Anchor0) <- c("Cx", "Cy")
  #
  # get split
  list_split <- f_pnts_group_end_split(pnts_bnd, pnt0, end_Anchor0)
  bnd_sides_right <- list_split[[1]]
  bnd_sides_left <- list_split[[2]]
  #
  # return
  return(list(bnd_sides_right, bnd_sides_left, GrpBndAgl_end))
}
#
#
#
f_pnts_bqa <- function(pnts, Agl0, pnt0, Aglx) {
  #
  #
  # bqa, bounding quadrilateral area
  #
  #
  # get Agls
  Agl2 <- Aglx[1]
  # get Agls
  Agl1_ratio <- Aglx[2]
  Agl3_ratio <- Aglx[3]
  # get Agls
  Agl1_min <- max(Agl2-pi, Agl0)
  Agl3_min <- max(Agl0+pi*2-pi, Agl2)
  # get Agls
  Agl1_max <- min(Agl0+pi, Agl2)
  Agl3_max <- min(Agl2+pi, Agl0+pi*2)
  # get Agls
  Agl1 <- Agl1_min + (Agl1_max-Agl1_min) * Agl1_ratio
  Agl3 <- Agl3_min + (Agl3_max-Agl3_min) * Agl3_ratio
  #
  # get agls
  Agls <- c(Agl0, Agl1, Agl2, Agl3)
  #
  #
  # get ints
  list_ints <- f_pnts_bqc(pnts, Agls, pnt0)
  #
  # get area
  area <- f_pnts_bqca(list_ints)
  #
  #
  # return
  return(area)
}
#
#
#
f_pnts_bqc <- function(pnts, Agls, pnt0) {
  #
  #
  # bqc, bounding quadrilateral corners
  #
  #
  # get agls
  Agl0 <- Agls[1]
  Agl1 <- Agls[2]
  Agl2 <- Agls[3]
  Agl3 <- Agls[4]
  #
  #
  # get tangency point
  pnt1 <- f_pnts_tangency(pnts, Agl0, Agl1)
  # check
  if (is.null(pnt1)) { return(Inf) }
  #
  # get tangency point
  pnt2 <- f_pnts_tangency(pnts, Agl1, Agl2)
  # check
  if (is.null(pnt2)) { return(Inf) }
  #
  # get tangency point
  pnt3 <- f_pnts_tangency(pnts, Agl2, Agl3)
  # check
  if (is.null(pnt3)) { return(Inf) }
  #
  #
  # get pntb
  pnt0b <- f_pnts_pntb(pnt0, Agl0, 100)
  # get pntb
  pnt1b <- f_pnts_pntb(pnt1, Agl1, 100)
  # get pntb
  pnt2b <- f_pnts_pntb(pnt2, Agl2, 100)
  # get pntb
  pnt3b <- f_pnts_pntb(pnt3, Agl3, 100)
  #
  #
  # get int
  int01 <- f_intersection_curves2(rbind(pnt0, pnt0b), rbind(pnt1, pnt1b))
  # check
  if (is.null(int01)) { return(Inf) }
  #
  # get int
  int12 <- f_intersection_curves2(rbind(pnt1, pnt1b), rbind(pnt2, pnt2b))
  # check
  if (is.null(int12)) { return(Inf) }
  #
  # get int
  int23 <- f_intersection_curves2(rbind(pnt2, pnt2b), rbind(pnt3, pnt3b))
  # check
  if (is.null(int23)) { return(Inf) }
  #
  # get int
  int30 <- f_intersection_curves2(rbind(pnt3, pnt3b), rbind(pnt0, pnt0b))
  # check
  if (is.null(int30)) { return(Inf) }
  #
  # debugging
  if (F) {
    #
    if ((abs(int01[1] - int12[1]) <= 10^-3 && abs(int01[2] == int12[2]) <= 10^-3) ||
        (abs(int01[1] - int23[1]) <= 10^-3 && abs(int01[2] == int23[2]) <= 10^-3) ||
        (abs(int01[1] - int30[1]) <= 10^-3 && abs(int01[2] == int30[2]) <= 10^-3)) {
      #
      # save
      f_pnts_save_points(rbind(pnt0, pnt1, pnt2, pnt3), NA, "pnt")
      f_pnts_save_points(rbind(pnt0b, pnt1b, pnt2b, pnt3b), NA, "pntb")
      f_pnts_save_points(rbind(int01, int12, int23, int30), NA, "int")
      #
      # return
      return(NULL)
    }
  }
  #
  # return
  return(list(int01, int12, int23, int30))
}
#
#
#
f_pnts_bqca <- function(list_ints) {
  #
  # get ints
  int01 <- list_ints[[1]]
  int12 <- list_ints[[2]]
  int23 <- list_ints[[3]]
  int30 <- list_ints[[4]]
  #
  # get sides
  a <- sqrt((int12[1]-int01[1])^2 + (int12[2]-int01[2])^2)
  b <- sqrt((int23[1]-int12[1])^2 + (int23[2]-int12[2])^2)
  c <- sqrt((int30[1]-int23[1])^2 + (int30[2]-int23[2])^2)
  d <- sqrt((int01[1]-int30[1])^2 + (int01[2]-int30[2])^2)
  # get diagonal
  p <- sqrt((int23[1]-int01[1])^2 + (int23[2]-int01[2])^2)
  #
  # get area triangle 1
  s <- (a + b + p) / 2
  area1 <- sqrt(s * (s-a) * (s-b) * (s-p))
  # get area triangle 2
  s <- (c + d + p) / 2
  area2 <- sqrt(s * (s-c) * (s-d) * (s-p))
  #
  # return
  return(area1 + area2)
}
#
#
#
f_pnts_tangency <- function(pnts, Agl0, Agl) {
  #
  # get count
  count_pnt <- nrow(pnts)
  #
  # get perpendicular
  perp_x <- cos(Agl - pi/2)
  perp_y <- sin(Agl - pi/2)
  #
  # get direction of the former edge
  dir_x <- cos(Agl0)
  dir_y <- sin(Agl0)
  #
  # check direction
  if ((perp_x*dir_x + perp_y*dir_y) == 0) {
    #
    # return
    return(NULL)
  }
  else if ((perp_x*dir_x + perp_y*dir_y) < 0) {
    #
    # reverse
    perp_x <- -perp_x
    perp_y <- -perp_y
  }
  #
  # initial, dot production
  dots <- c()
  #
  # get dots
  for (i in 1:count_pnt) {
    #
    # append
    dots <- rbind(dots, pnts[i, "Cx"] * perp_x + pnts[i, "Cy"] * perp_y)
  }
  #
  # get idx, possibly more than 1
  idx <- which(dots == max(dots))
  #
  # return
  return(pnts[idx[1], ])
}
#
#
#
f_pnts_pntb <- function(pnt, Agl, step) {
  #
  # get pnts
  pnt1_x <- pnt[1]
  pnt1_y <- pnt[2]
  #
  # get dirs
  dir_x <- cos(Agl)
  dir_y <- sin(Agl)
  #
  # get pnt1
  if (dir_x == 0) {
    #
    pnt2_x <- pnt1_x
    # pnt2_y <- pnt1_y * 2
    pnt2_y <- pnt1_y + step
  }
  else {
    #
    # pnt2_x <- pnt1_x * 2
    pnt2_x <- pnt1_x + step
    pnt2_y <- dir_y / dir_x * (pnt2_x - pnt1_x) + pnt1_y
  }
  #
  # get pnt2b
  pnt2 <- data.frame(pnt2_x, pnt2_y)
  colnames(pnt2) <- c("Cx", "Cy")
  #
  # return
  return(pnt2)
}
#
#
#
f_angle_included <- function(agl1, agl2) {
  #
  # get
  agl_icld <- agl2 - agl1
  #
  # check
  if (agl_icld > pi) { agl_icld <- agl_icld - pi*2 }
  if (agl_icld < -pi) { agl_icld <- agl_icld + pi*2 }
  #
  # return
  return(abs(agl_icld))
}
#
#
#
f_pnts_group_end_split <- function(pnts_bnd, pnt0, end_Anchor0) {
  #
  # get x0 and y0
  x0 <- pnt0[1, "Cx"]
  y0 <- pnt0[1, "Cy"]
  #
  # get profile line
  lx_prfl <- end_Anchor0[1, "Cx"] - x0
  ly_prfl <- end_Anchor0[1, "Cy"] - y0
  ln_prfl <- sqrt(lx_prfl^2 + ly_prfl^2)
  #
  # get sorted
  pnts_bnd_ascending <- f_pnts_bnd_ascending(pnts_bnd)
  #
  # initial, remove pnt0 from pnts, might be redundant
  pnts_bnd_ascending_pnt0_rm <- c()
  #
  # initial, included angle
  column_agl_icld <- c()
  #
  # get angles
  for (i in 1:nrow(pnts_bnd_ascending)) {
    #
    # get xy
    lx <- as.numeric(pnts_bnd_ascending[i, "Cx"] - x0)
    ly <- as.numeric(pnts_bnd_ascending[i, "Cy"] - y0)
    ln <- sqrt(lx^2 + ly^2)
    #
    # if not the same
    if (lx != 0 || ly != 0) {
      #
      # append
      pnts_bnd_ascending_pnt0_rm <- rbind(pnts_bnd_ascending_pnt0_rm, pnts_bnd_ascending[i, ])
      #
      # dot product
      dp <- lx * lx_prfl + ly * ly_prfl
      #
      # included angle
      agl_icld_cos <- dp / ln / ln_prfl
      agl_icld <- acos(agl_icld_cos)
      #
      # numerical error, necessary?
      if (agl_icld_cos > 1) { agl_icld <- 0 }
      if (agl_icld_cos < -1) { agl_icld <- pi }
      #
      # append
      column_agl_icld <- rbind(column_agl_icld, agl_icld)
    }
  }
  #
  # get index(s)
  idxs <- which(column_agl_icld == min(column_agl_icld))
  #
  # check
  if (length(idxs) <= 1) {
    #
    # get index
    idx <- idxs
  }
  else {
    #
    # initial
    column_dist <- c()
    #
    for (i in 1:length(idxs)) {
      #
      # get dist
      x <- pnts_bnd_ascending_pnt0_rm[idxs[i], "Cx"]
      y <- pnts_bnd_ascending_pnt0_rm[idxs[i], "Cy"]
      dist <- sqrt((x - x0)^2 + (y - y0)^2)
      column_dist <- rbind(column_dist, dist)
    }
    #
    # get minimum
    idx <- idxs[which(column_dist == min(column_dist))]
  }
  #
  # get pnts bnd sides
  # both sides keep split pnt?
  # not necessarily real left and right, just for a mark here
  pnts_bnd_sides_right <- pnts_bnd_ascending_pnt0_rm[1:idx, ]
  pnts_bnd_sides_left <- pnts_bnd_ascending_pnt0_rm[idx:nrow(pnts_bnd_ascending_pnt0_rm), ]
  # #
  # # get end anchor
  # end_Anchor <- pnts_bnd_ascending_pnt0_rm[idx]
  #
  # return
  return(list(pnts_bnd_sides_right, pnts_bnd_sides_left))
}
#
#
#
f_pnts_path <- function(list_pnts, IDB_max, pPolygons, EndAlgorithm, Finalization) {
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
    # get launch anchor, ID
    ID_LAcr <- LaunchAnchor[, "ID"]
    # get launch anchor, IDB
    IDB_LAcr <- LaunchAnchor[, "IDB"]
    # get launch anchor, IDB
    IDB_LAcr1m <- IDB_LAcr - 1
    if (IDB_LAcr1m < 1) { IDB_LAcr1m <- IDB_max }
    # get launch anchor, IDB
    IDB_LAcr1p <- IDB_LAcr + 1
    if (IDB_LAcr1p > IDB_max) { IDB_LAcr1p <- 1 }
    # get launch anchor, ID
    ID_LAcr1m <- which(pnts_bnd[, "IDB"] == IDB_LAcr1m)
    ID_LAcr1p <- which(pnts_bnd[, "IDB"] == IDB_LAcr1p)
    #
    # remove launch anchor
    pnts_bnd_rm_LAcr <- pnts_bnd[-which(pnts_bnd[, "ID"] == ID_LAcr), ]
    #
    # get, for end group
    list_end <- f_pnts_bnd_end_mbb(pnts_bnd_rm_LAcr, pnts_bnd[c(ID_LAcr1m, ID_LAcr1p), c("Cx", "Cy")])
    # get, for end group
    bnd_sides_right <- list_end[[1]]
    bnd_sides_left <- list_end[[2]]
    # get, for end group
    GrpBndAgl_initial <- list_end[[3]]
    GrpBndAgl_distal <- list_end[[3]]
    #
    # get bnd sides
    list_bnd_sides <- list(bnd_sides_right, bnd_sides_left)
    #
    # switch, do not need here?
    list_bnd_sides <- f_pnts_bnd_sides_switch2(list_bnd_sides, LaunchAnchor, IDB_max)
    #
    # get bnd sides
    bnd_sides_right <- list_bnd_sides[[1]]
    bnd_sides_left <- list_bnd_sides[[2]]
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
    return(list(column_anchors, list(bnd_sides_right), list(bnd_sides_left), c(GrpBndAgl_initial, GrpBndAgl_distal)))
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
    # get pnts_i_bnd
    pnts_i_bnd <- pnts_i[which(pnts_i[, "IDB"] != 0), ]
    #
    # get bnd sides
    if ((i == 1) || (i == count_pnts)) {
      #
      # get sorted
      pnts_i_bnd_ascending <- f_pnts_bnd_ascending(pnts_i_bnd)
      #
      # check
      if (i == 1) { 
        #
        # get pnts0
        pnts0 <- rbind(pnts_i_bnd_ascending[1, c("Cx", "Cy")], pnts_i_bnd_ascending[nrow(pnts_i_bnd_ascending), c("Cx", "Cy")])
      }
      else {
        #
        # get pnts0
        pnts0 <- rbind(pnts_i_bnd_ascending[nrow(pnts_i_bnd_ascending), c("Cx", "Cy")], pnts_i_bnd_ascending[1, c("Cx", "Cy")])
      }
      #
      # Algorithm 1
      if (EndAlgorithm == 1) {
        #
        # even
        list_end <- f_pnts_bnd_end_even(pnts_i, pnts0)
        #
        # get, for end group
        bnd_sides_right <- list_end[[1]]
        bnd_sides_left <- list_end[[2]]
        #
        # check
        if (Finalization) {
          #
          # mbb
          list_end <- f_pnts_bnd_end_mbb(pnts_i_bnd, pnts0)
          #
          # get, for end group
          GrpBndAgl <- list_end[[3]]
        }
      }
      #
      # Algorithm 2
      if (EndAlgorithm == 2) {
        #
        # even
        list_end <- f_pnts_bnd_end_even(pnts_i, pnts0)
        #
        # get, for end group
        bnd_sides_right <- list_end[[1]]
        bnd_sides_left <- list_end[[2]]
        #
        # check
        if (Finalization) {
          #
          # quad
          list_end <- f_pnts_bnd_end_quad(pnts_i_bnd, pnts0)
          #
          # get, for end group
          GrpBndAgl <- list_end[[3]]
        }
      }
      #
      # Algorithm 3
      if (EndAlgorithm == 3) {
        #
        # mbb
        list_end <- f_pnts_bnd_end_mbb(pnts_i_bnd, pnts0)
        #
        # get, for end group
        bnd_sides_right <- list_end[[1]]
        bnd_sides_left <- list_end[[2]]
        #
        # check
        if (Finalization) {
          # #
          # # mbb
          # list_end <- f_pnts_bnd_end_mbb(pnts_i_bnd, pnts0)
          #
          # get, for end group
          GrpBndAgl <- list_end[[3]]
        }
      }
      #
      # Algorithm 4
      if (EndAlgorithm == 4) {
        #
        # mbb
        list_end <- f_pnts_bnd_end_mbb(pnts_i_bnd, pnts0)
        #
        # get, for end group
        bnd_sides_right <- list_end[[1]]
        bnd_sides_left <- list_end[[2]]
        #
        # check
        if (Finalization) {
          #
          # quad
          list_end <- f_pnts_bnd_end_quad(pnts_i_bnd, pnts0)
          #
          # get, for end group
          GrpBndAgl <- list_end[[3]]
        }
      }
      #
      # Algorithm 5
      if (EndAlgorithm == 5) {
        #
        # quad
        list_end <- f_pnts_bnd_end_quad(pnts_i_bnd, pnts0)
        #
        # get, for end group
        bnd_sides_right <- list_end[[1]]
        bnd_sides_left <- list_end[[2]]
        #
        # check
        if (Finalization) {
          #
          # mbb
          list_end <- f_pnts_bnd_end_mbb(pnts_i_bnd, pnts0)
          #
          # get, for end group
          GrpBndAgl <- list_end[[3]]
        }
      }
      #
      # Algorithm 6
      if (EndAlgorithm == 6) {
        #
        # quad
        list_end <- f_pnts_bnd_end_quad(pnts_i_bnd, pnts0)
        #
        # get, for end group
        bnd_sides_right <- list_end[[1]]
        bnd_sides_left <- list_end[[2]]
        #
        # check
        if (Finalization) {
          # #
          # # quad
          # list_end <- f_pnts_bnd_end_quad(pnts_i_bnd, pnts0)
          #
          # get, for end group
          GrpBndAgl <- list_end[[3]]
        }
      }
      #
      # get bnd sides
      list_pnts_i_bnd_sides <- list(bnd_sides_right, bnd_sides_left)
      #
      # check
      if (i == 1) { 
        #
        # get, for end group
        if (Finalization) { GrpBndAgl_initial <- GrpBndAgl }
        else { GrpBndAgl_initial <- NULL }
        #
        # get launch anchor
        # before switch, bnd_sides_right is split at idx
        LaunchAnchor <- bnd_sides_right[nrow(bnd_sides_right), ]
      }
      else {
        #
        # get, for end group
        if (Finalization) { GrpBndAgl_distal <- GrpBndAgl }
        else { GrpBndAgl_distal <- NULL }
      }
    }
    else {
      #
      # get bnd sides
      list_pnts_i_bnd_sides <- f_pnts_bnd_sides(pnts_i_bnd, IDB_max)
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
    # check corners
    pnts_i_bnd_sides_corners_check <- f_pnts_bnd_sides_corners_check(list_pnts_i_bnd_sides)
    #
    # check corners
    if (i == 1) {
      #
      # check
      if ((pnts_i_bnd_sides_corners_check[2] == F)) {
        #
        # return
        return(NULL)
      }
    }
    else if (i == count_pnts) {
      #
      # check
      if ((pnts_i_bnd_sides_corners_check[1] == F)) {
        #
        # return
        return(NULL)
      }
    }
    else {
      #
      # check
      if ((pnts_i_bnd_sides_corners_check[1] == F) ||
          (pnts_i_bnd_sides_corners_check[2] == F)) {
        #
        # return
        return(NULL)
      }
    }
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
  return(list(data.frame(column_anchors), list_bnd_sides_right, list_bnd_sides_left, c(GrpBndAgl_initial, GrpBndAgl_distal)))
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
f_pnts_path_lines <- function(list_bnd_sides_right, list_bnd_sides_left, column_anchors, GrpBndAgl_ends, count_strip) {
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
  if (count_grp >= 2) {
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
        # debugging
        if (is.null(ints_upper)) {
          #
          # warning
          warning("Error: get intersection failed.")
          #
          # return
          return(NULL)
        }
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
        # debugging
        if (is.null(ints_lower)) {
          #
          # warning
          warning("Error: get intersection failed.")
          #
          # return
          return(NULL)
        }
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
        # debugging
        if (is.null(ints_lower)) {
          #
          # warning
          warning("Error: get intersection failed.")
          #
          # return
          return(NULL)
        }
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
        # debugging
        if (is.null(ints_upper)) {
          #
          # warning
          warning("Error: get intersection failed.")
          #
          # return
          return(NULL)
        }
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
  if (count_grp >= 2) {
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
      GrpBndAgl <- atan2(GrpBnd[2, "Cy"]-GrpBnd[1, "Cy"], GrpBnd[2, "Cx"]-GrpBnd[1, "Cx"])
      #
      # append
      column_GrpBndAgls <- rbind(column_GrpBndAgls, GrpBndAgl)
    }
  }
  else {
    #
    # Get a single GrpBndAgl for one group situation
    # also, direct from left side to right side
    #
    # get
    dx <- column_anchors[3, "Cx"] - column_anchors[1, "Cx"]
    dy <- column_anchors[3, "Cy"] - column_anchors[1, "Cy"]
    #
    # get angle
    GrpBndAgl_for_one_group <- atan2(dy, dx) - pi/2
  }
  #
  # get GrpBndAgl_ends
  GrpBndAgl_end_initial <- GrpBndAgl_ends[1]
  GrpBndAgl_end_distal <- GrpBndAgl_ends[2]
  #
  # get split stations
  list_split_stations <- f_pnts_split_stations(pnts_curve_profile, count_strip)
  # get
  strip_station_interval <- list_split_stations[[1]]
  # get
  strip_stations <- list_split_stations[[2]]
  #
  # initial
  strip_stations_angles <- c()
  #
  # initial
  column_anchors_GrpBndInts_stations <- column_anchors_GrpBndInts
  column_anchors_GrpBndInts_stations$station <- c(0)
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
    # check, if, no group boundary, possibly only one group
    if (count_grp_bnd >= 2) {
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
        # get length difference
        Lhrz_diff <- column_GrpBndInts[index_GrpBndInt_lower, "Lhrz"] - 0
        Lhrz_diff_station <- strip_station_Lhrz - 0
        #
        # get angle
        angle_grp_bnd_upper <- GrpBndAgl_end_initial
        angle_grp_bnd_lower <- column_GrpBndAgls[index_GrpBndInt_lower]
      }
      #
      # if within the distal group
      if (index_GrpBndInt_lower == count_grp_bnd+1) {
        #
        # get length difference
        Lhrz_diff <- column_anchors[count_anchor, "Lhrz"] - column_GrpBndInts[index_GrpBndInt_upper, "Lhrz"]
        Lhrz_diff_station <- strip_station_Lhrz - column_GrpBndInts[index_GrpBndInt_upper, "Lhrz"]
        #
        # get angle
        angle_grp_bnd_upper <- column_GrpBndAgls[index_GrpBndInt_upper]
        angle_grp_bnd_lower <- GrpBndAgl_end_distal
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
      }
    }
    else {
      #
      # get
      if (strip_station_Lhrz <= column_anchors[2, "Lhrz"]) {
        #
        # get length difference
        Lhrz_diff <- column_anchors[2, "Lhrz"] - 0
        Lhrz_diff_station <- strip_station_Lhrz - 0
        #
        # get angle
        angle_grp_bnd_upper <- GrpBndAgl_end_initial
        angle_grp_bnd_lower <- GrpBndAgl_for_one_group
      }
      else {
        #
        # get length difference
        Lhrz_diff <- column_anchors[3, "Lhrz"] - column_anchors[2, "Lhrz"]
        Lhrz_diff_station <- strip_station_Lhrz - column_anchors[2, "Lhrz"]
        #
        # get angle
        angle_grp_bnd_upper <- GrpBndAgl_for_one_group
        angle_grp_bnd_lower <- GrpBndAgl_end_distal
      }
    }
    #
    # get angle difference
    angle_diff <- angle_grp_bnd_lower - angle_grp_bnd_upper
    #
    # update
    if (angle_diff > pi*2) { angle_diff <- angle_diff - pi*2 }
    if (angle_diff < -pi*2) { angle_diff <- angle_diff + pi*2 }
    # update
    if (angle_diff > pi) { angle_diff <- angle_diff - pi*2 }
    if (angle_diff < -pi) { angle_diff <- angle_diff + pi*2 }
    #
    # get angle difference at the station
    angle_diff_station <- angle_diff / Lhrz_diff * Lhrz_diff_station
    #
    # get strip angle
    strip_angle <- angle_grp_bnd_upper + angle_diff_station
    #
    # update
    if (strip_angle > pi*2) { strip_angle <- strip_angle - 2*pi }
    if (strip_angle < -pi*2) { strip_angle <- strip_angle + 2*pi }
    # update
    if (strip_angle > pi) { strip_angle <- strip_angle - 2*pi }
    if (strip_angle < -pi) { strip_angle <- strip_angle + 2*pi }
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
  # return
  return(list(list_GrpBnds, strip_stations, strip_stations_angles, column_anchors_GrpBndInts_stations))
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
f_pnts_split_stations <- function(pnts, count_split = 1) {
  #
  # check
  stopifnot(count_split > 0)
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
  # get stationing of segments
  stationing <- seq(from = 0, to = length_total, length.out = count_split+1)
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
f_strips <- function(pPolygons, list_path, count_strip) {
  #
  # get anchors
  column_anchors <- list_path[[1]]
  # get bnd sides
  list_bnd_sides_right <- list_path[[2]]
  list_bnd_sides_left <- list_path[[3]]
  # get group boundary angle
  GrpBndAgl_ends <- list_path[[4]]
  #
  # get path lines
  list_path_lines <- f_pnts_path_lines(list_bnd_sides_right, list_bnd_sides_left, column_anchors, GrpBndAgl_ends, count_strip)
  #
  # get 
  strip_GrpBnds <- list_path_lines[[1]]
  # get
  strip_stations <- list_path_lines[[2]]
  strip_stations_angles <- list_path_lines[[3]]
  # get 
  strip_nodes <- list_path_lines[[4]]
  #
  # initial
  strip_nodes$Czd <- c(0)
  strip_nodes$Lalld <- c(0)
  #
  # switch
  strip_nodes <- strip_nodes[, c("Cx", "Cy", "Czd", "Czp", "type", "grp", "station", "Lhrz", "Lalld", "Lallp")]
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
    f_pnts_save_points(strip_stations, NA, "strip_stations")
  }
  #
  # get clipping lines, sfc
  sfc_spldf_cls <- sf::st_as_sfc(spldf_cls)
  #
  # get count, there might be stations outside of landslide polygon
  count_cl <- length(sfc_spldf_cls)
  #
  # initial
  strip_stations_IDs_inside <- c()
  strip_nodes_IDs_outside <- c()
  #
  # get SpatialPolygons
  pSpatialPolygons <- sp::SpatialPolygons(list(pPolygons))
  #
  # check
  for (i in 1:count_cl) {
    #
    # get coords
    spnt <- strip_stations[i, c("Cx", "Cy")]
    #
    # get shp
    sshp <- f_pnts2shp_points(spnt, NA)
    #
    # point in polygon or not
    sshp_over <- sp::over(sshp, pSpatialPolygons)
    #
    # check
    if (!(is.na(sshp_over))) {
      #
      # append
      strip_stations_IDs_inside <- rbind(strip_stations_IDs_inside, i)
    }
    else {
      #
      # append
      strip_nodes_IDs_outside <- rbind(strip_nodes_IDs_outside, which(strip_nodes[, "station"] == i))
    }
  }
  #
  # get
  sfc_spldf_cls_for_sp <- sfc_spldf_cls[strip_stations_IDs_inside]
  # get
  strip_stations_for_sp <- strip_stations[strip_stations_IDs_inside, ]
  strip_stations_angles_for_sp <- strip_stations_angles[strip_stations_IDs_inside]
  # get
  if (length(strip_nodes_IDs_outside) >= 1) { strip_nodes_for_sp <- strip_nodes[-strip_nodes_IDs_outside, ] }
  else { strip_nodes_for_sp <- strip_nodes }
  #
  # get count, there might be invalid clipping polyline
  count_cl <- length(sfc_spldf_cls_for_sp)
  #
  # initial
  strip_stations_IDs_for_sp <- strip_stations_IDs_inside
  #
  # initial
  bool_ints <- TRUE
  #
  # check clipping lines
  while (bool_ints && count_cl >= 2) {
    #
    # initial
    bool_ints <- FALSE
    #
    # initial
    column_ints_count <- c(1:count_cl) - c(1:count_cl)
    #
    # get
    for (i in 1:(count_cl-1)) {
      #
      # get clipping polyline
      sfc_cl_i <- sfc_spldf_cls_for_sp[[i]]
      #
      # get coords
      sfc_cl_i_coords <- matrix(sfc_cl_i, nc = 2, byrow = F)
      #
      # check
      for (k in (i+1):count_cl) {
        #
        # get clipping polyline
        sfc_cl_k <- sfc_spldf_cls_for_sp[[k]]
        #
        # get coords
        sfc_cl_k_coords <- matrix(sfc_cl_k, nc = 2, byrow = F)
        #
        # get ints
        ints <- f_intersection_curves(sfc_cl_i_coords, sfc_cl_k_coords)
        #
        # check
        if (!is.null(ints)) {
          #
          # set
          bool_ints <- TRUE
          #
          # update
          column_ints_count[i] <- column_ints_count[i] + 1
          column_ints_count[k] <- column_ints_count[k] + 1
        }
      }
    }
    #
    # check
    if (max(column_ints_count) >= 1) {
      #
      # set
      bool_ints <- TRUE
      #
      # get
      idxs <- which(column_ints_count == max(column_ints_count))
      #
      # initial
      agl_diffs <- c()
      #
      # choose
      for (i in 1:length(idxs)) {
        #
        # get idx
        idxs_i <- idxs[i]
        #
        # initial
        agl_diff_upper <- 0
        agl_diff_lower <- 0
        #
        # get upper
        if (idxs_i-1 >= 1) {
          #
          # get
          agl_diff_upper <- strip_stations_angles_for_sp[idxs_i] - strip_stations_angles_for_sp[idxs_i-1]
          if (agl_diff_upper >= 360) { agl_diff_upper <- agl_diff_upper - 360}
          if (agl_diff_upper <= -360) { agl_diff_upper <- agl_diff_upper + 360}
        }
        #
        # get lower
        if (idxs_i+1 <= count_cl) {
          #
          # get
          agl_diff_lower <- strip_stations_angles_for_sp[idxs_i] - strip_stations_angles_for_sp[idxs_i+1]
          if (agl_diff_lower >= 360) { agl_diff_lower <- agl_diff_lower - 360}
          if (agl_diff_lower <= -360) { agl_diff_lower <- agl_diff_lower + 360}
        }
        #
        # get
        agl_diff_i <- abs(agl_diff_upper) + abs(agl_diff_lower)
        #
        # append
        agl_diffs <- rbind(agl_diffs, agl_diff_i)
      }
      #
      # get idx
      idx <- idxs[which(agl_diffs == max(agl_diffs))]
      #
      # check, possibly still more than 1
      if (length(idx) > 1) { idx <- idx[1] }
      #
      # append
      sfc_spldf_cls_for_sp <- sfc_spldf_cls_for_sp[-idx]
      #
      # append
      strip_stations_for_sp <- strip_stations_for_sp[-idx, ]
      #
      # append
      strip_stations_angles_for_sp <- strip_stations_angles_for_sp[-idx]
      #
      # get ID
      ID <- strip_stations_IDs_for_sp[idx, 1]
      #
      # append
      strip_stations_IDs_for_sp <- data.frame(strip_stations_IDs_for_sp[-idx, ])
      #
      # get index
      idx_station_in_nodes <- which(strip_nodes_for_sp[, "station"] == ID)
      #
      # append
      strip_nodes_for_sp <- strip_nodes_for_sp[-idx_station_in_nodes, ]
    }
    #
    # update count
    count_cl <- length(sfc_spldf_cls_for_sp)
  }
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
  count_strip <- 1
  #
  # get count
  count_cl <- length(sfc_spldf_cls_for_sp)
  #
  # get SpatialPolygons
  # and, update station Cx and Cy, it will slightly change
  # when consider station as the intersection of profile and strip
  # also update initial and distal group anchor
  for (i in 1:count_cl) {
    #
    # get clipping polyline
    sfc_cl <- sfc_spldf_cls_for_sp[[i]]
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
      coords <- matrix(sfc_cl, nc = 2, byrow = F)
      sfc_spl <- f_coords2spl(coords)
      # save
      f_save_list2sp(list(sfc_spl), NULL, NA, "sfc_spl")
      #
      # save
      f_pnts_save_points(strip_nodes_for_sp, NA, "strip_nodes_for_sp")
      #
      # save
      f_pnts_save_points(strip_stations_for_sp, NA, "strip_stations_for_sp")
      #
      # warning
      warning("Error: count of POLYGON is not 2.")
      #
      # pause at here
      return(NULL)
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
    line_station <- strip_stations_for_sp[i, 1:2]
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
    pnt_upper <- strip_stations_for_sp[i, 1:2]
    #
    # get coords, strip
    strip_sp <- list_strip_sp[[length(list_strip_sp)]]
    pnts_strip <- strip_sp@polygons[[1]]@Polygons[[1]]@coords
    pnts_strip <- data.frame(pnts_strip)
    colnames(pnts_strip) <- c("Cx", "Cy")
    #
    # get ints, profile with strip
    int.coords <- f_intersection_curves(strip_nodes_for_sp[, c("Cx", "Cy")], pnts_strip)
    #
    # debugging
    if (is.null(int.coords)) {
      #
      # save
      f_save_list2sp(list(strip_sp), NULL, NA, "strip")
      #
      # save
      f_pnts_save_points(strip_nodes, NA, "strip_nodes")
      #
      # save
      f_pnts_save_points(strip_nodes_for_sp, NA, "strip_nodes_for_sp")
      #
      # save
      f_pnts_save_points(int.coords, NA, "int.coords")
      #
      # warning
      warning("Error: profile and strip do not have intersections.")
      #
      # return
      return(NULL)
    }
    #
    # colnames
    int.coords <- data.frame(int.coords)
    colnames(int.coords) <- c("Cx", "Cy")
    #
    # debugging
    if (nrow(int.coords) >= 3) {
      # #
      # # save
      # f_save_list2sp(list(strip_sp), NULL, NA, "strip")
      # #
      # # save
      # f_pnts_save_points(strip_nodes, NA, "strip_nodes")
      # #
      # # save
      # f_pnts_save_points(strip_nodes_for_sp, NA, "strip_nodes_for_sp")
      # #
      # # save
      # f_pnts_save_points(int.coords, NA, "int.coords")
      #
      # warning
      warning("Error: profile and strip have more than 2 intersections.")
      # #
      # # return
      # return(NULL)
      #
      # keep the first and last intersections?
      # sometimes, profile will go through outside landslide polygon?
      int.coords <- int.coords[c(1, nrow(int.coords)), ]
    }
    #
    # get index
    idx_station_in_nodes <- which(strip_nodes_for_sp[, "station"] == strip_stations_IDs_for_sp[i, 1])
    #
    # check
    if (nrow(int.coords) == 1) {
      #
      # update Cx and Cy of station
      strip_nodes_for_sp[idx_station_in_nodes, "Cx"] <- int.coords[1, 1]
      strip_nodes_for_sp[idx_station_in_nodes, "Cy"] <- int.coords[1, 2]
    }
    #
    # check
    if (nrow(int.coords) == 2) {
      #
      # get dists
      dist1 <- sqrt((strip_nodes_for_sp[idx_station_in_nodes, "Cx"] - int.coords[1, 1])^2 + (strip_nodes_for_sp[idx_station_in_nodes, "Cy"] - int.coords[1, 2])^2)
      dist2 <- sqrt((strip_nodes_for_sp[idx_station_in_nodes, "Cx"] - int.coords[2, 1])^2 + (strip_nodes_for_sp[idx_station_in_nodes, "Cy"] - int.coords[2, 2])^2)
      #
      # update Cx and Cy of station
      if (dist1 <= dist2) {
        #
        # replace
        strip_nodes_for_sp[idx_station_in_nodes, "Cx"] <- int.coords[1, 1]
        strip_nodes_for_sp[idx_station_in_nodes, "Cy"] <- int.coords[1, 2]
      }
      else {
        #
        # replace
        strip_nodes_for_sp[idx_station_in_nodes, "Cx"] <- int.coords[2, 1]
        strip_nodes_for_sp[idx_station_in_nodes, "Cy"] <- int.coords[2, 2]
      }
      #
      # check, if we have the first strip
      if (length(list_strip_sp) == 1) {
        #
        # update Cx and Cy of the initial group anchor
        if (dist1 <= dist2) {
          #
          # replace
          strip_nodes_for_sp[1, "Cx"] <- int.coords[2, 1]
          strip_nodes_for_sp[1, "Cy"] <- int.coords[2, 2]
        }
        else {
          #
          # replace
          strip_nodes_for_sp[1, "Cx"] <- int.coords[1, 1]
          strip_nodes_for_sp[1, "Cy"] <- int.coords[1, 2]
        }
      }
    }
  }
  #
  # append
  list_strip_sp[[count_strip]] <- f_coords2spg((sfc_SpatialPolygons[[1]])[[1]])
  #
  # get count
  count_node <- nrow(strip_nodes_for_sp)
  #
  # update Cx and Cy for the distal group anchor
  if (T) {
    #
    # initial
    Cx <- strip_nodes_for_sp[count_node, "Cx"]
    Cy <- strip_nodes_for_sp[count_node, "Cy"]
    #
    # until Cx and Cy do not update
    while (T) {
      #
      # get coords, strip
      strip_sp <- list_strip_sp[[count_strip]]
      pnts_strip <- strip_sp@polygons[[1]]@Polygons[[1]]@coords
      pnts_strip <- data.frame(pnts_strip)
      colnames(pnts_strip) <- c("Cx", "Cy")
      #
      # coords <- pPolygons@Polygons[[1]]@coords
      #
      # get ints, profile with strip
      int.coords <- f_intersection_curves(strip_nodes_for_sp[, c("Cx", "Cy")], pnts_strip)
      #
      # debugging
      if (is.null(int.coords) || nrow(int.coords) >= 3) {
        #
        # save
        f_save_list2sp(list(strip_sp), NULL, NA, "strip")
        #
        # save
        f_pnts_save_points(strip_nodes, NA, "strip_nodes")
        #
        # save
        f_pnts_save_points(strip_nodes_for_sp, NA, "strip_nodes_for_sp")
        #
        # save
        f_pnts_save_points(int.coords, NA, "int.coords")
        #
        # warning
        warning("Error: profile and strip do not have 2 intersections.")
        #
        # return
        return(NULL)
      }
      #
      # colnames
      int.coords <- data.frame(int.coords)
      colnames(int.coords) <- c("Cx", "Cy")
      #
      # check
      if (nrow(int.coords) == 1) {
        # #
        # # update Cx and Cy of station
        # strip_nodes_for_sp[count_node, "Cx"] <- int.coords[1, 1]
        # strip_nodes_for_sp[count_node, "Cy"] <- int.coords[1, 2]
      }
      #
      # check
      if (nrow(int.coords) == 2) {
        #
        # get dists
        dist1 <- sqrt((strip_nodes_for_sp[count_node, "Cx"] - int.coords[1, 1])^2 + (strip_nodes_for_sp[count_node, "Cy"] - int.coords[1, 2])^2)
        dist2 <- sqrt((strip_nodes_for_sp[count_node, "Cx"] - int.coords[2, 1])^2 + (strip_nodes_for_sp[count_node, "Cy"] - int.coords[2, 2])^2)
        #
        # update Cx and Cy of the distal group anchor
        if (dist1 <= dist2) {
          #
          # replace
          strip_nodes_for_sp[count_node, "Cx"] <- int.coords[1, 1]
          strip_nodes_for_sp[count_node, "Cy"] <- int.coords[1, 2]
        }
        else {
          #
          # replace
          strip_nodes_for_sp[count_node, "Cx"] <- int.coords[2, 1]
          strip_nodes_for_sp[count_node, "Cy"] <- int.coords[2, 2]
        }
      }
      #
      # check
      if ((Cx - strip_nodes_for_sp[count_node, "Cx"]) == 0 && 
          (Cy - strip_nodes_for_sp[count_node, "Cy"]) == 0) {
        #
        # break
        break
      }
      else {
        #
        # update
        Cx <- strip_nodes_for_sp[count_node, "Cx"]
        Cy <- strip_nodes_for_sp[count_node, "Cy"]
      }
    }
  }
  #
  # return
  return(list(list_strip_sp, strip_nodes_for_sp))
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
    pCoords1_extnd <- f_pnts_nearest_extending(c(line_x0, line_y0), pCoords1, 10^3)
    #
    # get pnts, line
    pnts_line2 <- matrix(rbind(c(line_x0, line_y0), c(line_x2, line_y2)), nc = 2, byrow = F)
    #
    # get ints
    pCoords2 <- f_intersection_curves(pnts_line2, pPolygons_coords)
    #
    # get nearest extending
    pCoords2_extnd <- f_pnts_nearest_extending(c(line_x0, line_y0), pCoords2, 10^3)
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
  # get area2d
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
  # initial
  column_coords <- c()
  #
  # initial
  column_triangles <- c()
  #
  # initial
  column_Asrfs <- data.frame()
  column_drops <- data.frame()
  #
  # initial
  list_triangle_sp <- list()
  # initial
  column_triangles_df <- data.frame()
  #
  # get triangles
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
        return(NULL)
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
    tragls <- t(tragls)
    #
    # append
    column_triangles <- rbind(column_triangles, tragls+length(column_coords)/3)
    #
    # append
    column_coords <- rbind(column_coords, cbind(Cx, Cy, Cz))
    #
    # save triangles, set a switch
    if (T) {
      #
      # append
      for (k in 1:nrow(tragls)) {
        #
        # get index
        idx1 <- tragls[k, 1]
        idx2 <- tragls[k, 2]
        idx3 <- tragls[k, 3]
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
    area2d <- 0
    area3d <- 0
    #
    # get area
    for (k in 1:nrow(tragls)) {
      #
      # get index
      idx1 <- tragls[k, 1]
      idx2 <- tragls[k, 2]
      idx3 <- tragls[k, 3]
      #
      # get edge
      a <- sqrt((Cx[idx2]-Cx[idx1])^2 + (Cy[idx2]-Cy[idx1])^2)
      b <- sqrt((Cx[idx3]-Cx[idx2])^2 + (Cy[idx3]-Cy[idx2])^2)
      c <- sqrt((Cx[idx1]-Cx[idx3])^2 + (Cy[idx1]-Cy[idx3])^2)
      #
      # get area2d
      s <- (a+b+c)/2
      area2d_k <- sqrt(s*(s-a)*(s-b)*(s-c))
      #
      # get edge
      a <- sqrt((Cx[idx2]-Cx[idx1])^2 + (Cy[idx2]-Cy[idx1])^2 + (Cz[idx2]-Cz[idx1])^2)
      b <- sqrt((Cx[idx3]-Cx[idx2])^2 + (Cy[idx3]-Cy[idx2])^2 + (Cz[idx3]-Cz[idx2])^2)
      c <- sqrt((Cx[idx1]-Cx[idx3])^2 + (Cy[idx1]-Cy[idx3])^2 + (Cz[idx1]-Cz[idx3])^2)
      #
      # get area3d
      s <- (a+b+c)/2
      area3d_k <- sqrt(s*(s-a)*(s-b)*(s-c))
      #
      # update
      area2d <- area2d + area2d_k
      area3d <- area3d + area3d_k
      #
      # append
      column_triangles_df <- rbind(column_triangles_df, c(area2d_k, area3d_k))
    }
    #
    # append
    column_Asrfs <- rbind(column_Asrfs, area3d)
    #
    # debugging, checking area2d
    if (F) {
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
  }
  #
  # colnames
  colnames(column_Asrfs) <- c("Asrf")
  colnames(column_drops) <- c("Drop")
  # colnames
  colnames(column_triangles_df) <- c("Area2D", "Area3D")
  #
  # initial
  column_ints <- data.frame()
  #
  # get intersections with triangles
  for (i in 1:nrow(column_triangles)) {
    #
    # get ints, edge1
    idxo <- column_triangles[i, 1]
    idxp <- column_triangles[i, 2]
    pnts_edge <- rbind(column_coords[idxo, c("Cx", "Cy")], column_coords[idxp, c("Cx", "Cy")])
    int.coords <- f_intersection_curves(strip_nodes[, c("Cx", "Cy")], pnts_edge)
    #
    if (is.null(int.coords) == FALSE) { column_ints <- rbind(column_ints, cbind(int.coords, i)) }
    #
    # get ints, edge2
    idxo <- column_triangles[i, 2]
    idxp <- column_triangles[i, 3]
    pnts_edge <- rbind(column_coords[idxo, c("Cx", "Cy")], column_coords[idxp, c("Cx", "Cy")])
    int.coords <- f_intersection_curves(strip_nodes[, c("Cx", "Cy")], pnts_edge)
    #
    if (is.null(int.coords) == FALSE) { column_ints <- rbind(column_ints, cbind(int.coords, i)) }
    #
    # get ints, edge3
    idxo <- column_triangles[i, 3]
    idxp <- column_triangles[i, 1]
    pnts_edge <- rbind(column_coords[idxo, c("Cx", "Cy")], column_coords[idxp, c("Cx", "Cy")])
    int.coords <- f_intersection_curves(strip_nodes[, c("Cx", "Cy")], pnts_edge)
    #
    if (is.null(int.coords) == FALSE) { column_ints <- rbind(column_ints, cbind(int.coords, i)) }
  }
  #
  # colnames
  colnames(column_ints) <- c("Cx", "Cy", "triangle")
  #
  # remove duplicated
  idx_duplicated <- duplicated(rbind(strip_nodes[, c("Cx", "Cy")], column_ints[, c("Cx", "Cy")]))
  idx_duplicated <- idx_duplicated[(nrow(strip_nodes)+1):length(idx_duplicated)]
  # remove
  column_ints_rmdp <- column_ints[!idx_duplicated, ]
  #
  # initial, has additional attributes
  column_ints_rmdp2 <- c()
  #
  # initial
  column_positions <- c()
  #
  # get count
  count_node <- nrow(strip_nodes)
  #
  # get dists, must use Lhrz
  column_dists_nebr <- strip_nodes[2:count_node, "Lhrz"] - strip_nodes[1:(count_node-1), "Lhrz"]
  #
  # get positions
  for (i in 1:nrow(column_ints_rmdp)) {
    #
    # get
    int <- column_ints_rmdp[i, c("Cx", "Cy")]
    int$Czd <- 0
    int$Czp <- 0
    int$type <- "Intersection: Triangle boundary"
    int$grp <- 0
    int$station <- 0
    int$Lhrz <- 0
    int$Lalld <- 0
    int$Lallp <- 0
    #
    # append
    column_ints_rmdp2 <- rbind(column_ints_rmdp2, int)
    #
    # get Cx and Cy
    Cx <- int$Cx
    Cy <- int$Cy
    #
    # initial
    column_agls <- c()
    #
    # get included angle
    for (k in 1:(count_node-1)) {
      #
      # get vector
      vector_upper_x <- strip_nodes[k, "Cx"] - Cx
      vector_upper_y <- strip_nodes[k, "Cy"] - Cy
      #
      # get vector
      vector_lower_x <- strip_nodes[k+1, "Cx"] - Cx
      vector_lower_y <- strip_nodes[k+1, "Cy"] - Cy
      #
      # get length
      vector_upper_length <- sqrt(vector_upper_x^2+ vector_upper_y^2)
      vector_lower_length <- sqrt(vector_lower_x^2+ vector_lower_y^2)
      #
      # get included angle
      if (vector_upper_length == 0 || vector_lower_length == 0) {
        #
        # append
        column_agls <- rbind(column_agls, pi)
      }
      else {
        #
        # get included angle
        agl_cos <- vector_upper_x*vector_lower_x + vector_upper_y*vector_lower_y
        agl_cos <- agl_cos / vector_upper_length
        agl_cos <- agl_cos / vector_lower_length
        # get
        agl <- acos(agl_cos)
        #
        # numerical error, necessary?
        if (agl_cos > 1) { agl <- 0 }
        if (agl_cos < -1) { agl <- pi }
        #
        # append
        column_agls <- rbind(column_agls, agl)
      }
    }
    #
    # check position, do not use equal, because of numerical error?
    positions <- data.frame(which(column_agls == max(column_agls)))
    #
    # check count, possibly two positions have pi?
    if (nrow(positions) >= 2) {
      #
      # choose one
      positions <- positions[1, 1]
    }
    #
    # get dists
    column_dists_upper <- (Cx - strip_nodes[1:(count_node-1), "Cx"])^2
    column_dists_upper <- column_dists_upper + (Cy - strip_nodes[1:(count_node-1), "Cy"])^2
    column_dists_upper <- sqrt(column_dists_upper)
    #
    # get dists
    column_dists_lower <- (Cx - strip_nodes[2:count_node, "Cx"])^2
    column_dists_lower <- column_dists_lower + (Cy - strip_nodes[2:count_node, "Cy"])^2
    column_dists_lower <- sqrt(column_dists_lower)
    #
    # initial, must use relative difference
    column_dists_relative_diff <- c()
    #
    # check position, do not use equal, because of numerical error?
    for (k in 1:(count_node-1)) {
      #
      # append
      column_dists_relative_diff <- rbind(column_dists_relative_diff, abs(column_dists_upper[k] + column_dists_lower[k] - column_dists_nebr[k]) / column_dists_nebr[k])
    }
    #
    # check position, again
    if(column_dists_relative_diff[positions[1, 1]] != min(column_dists_relative_diff)) {
      # #
      # # use relative difference, still not robust
      # # drop this unnecessary check?
      # #
      # # save
      # f_pnts_save_points(strip_nodes, NA, "strip_nodes")
      # #
      # # save
      # f_pnts_save_points(int, NA, "int")
      # #
      # warning
      warning("warning: position for ints found not consistent.")
      # #
      # # pause at here
      # return(NULL)
    }
    #
    # append
    column_positions <- rbind(column_positions, positions[1, 1])
  }
  #
  # sort, for the first time
  column_ints_rmdp2 <- column_ints_rmdp2[order(column_positions), ]
  column_positions <- column_positions[order(column_positions)]
  # do not need to save triangle?
  # other types of nodes do not have triangle marked
  # anyway, you need to search one by one
  # column_ints_rmdp2_triangle <- column_ints_rmdp[order(column_positions), "triangle"]
  #
  # unique
  column_positions_unique <- unique(column_positions)
  #
  # sort, for the second time
  column_ints_rmdp2_ordered <- c()
  #
  # sort, for the second time
  for (i in 1:length(column_positions_unique)) {
    #
    # get index
    idxs <- which(column_positions == column_positions_unique[i])
    #
    # get count
    count_idxs <- length(idxs)
    #
    # check
    if (count_idxs >= 2) {
      #
      # initial
      column_dists_upper <- c()
      #
      # get dists from the upper node
      for (k in 1:count_idxs) {
        #
        # get Cx and Cy
        Cx <- column_ints_rmdp2[idxs[k], "Cx"]
        Cy <- column_ints_rmdp2[idxs[k], "Cy"]
        #
        # get dist
        dists_upper <- (Cx - strip_nodes[column_positions_unique[i], "Cx"])^2
        dists_upper <- dists_upper + (Cy - strip_nodes[column_positions_unique[i], "Cy"])^2
        dists_upper <- sqrt(dists_upper)
        #
        # append
        column_dists_upper <- rbind(column_dists_upper, dists_upper)
      }
      #
      # order idxs
      idxs_ordered <- idxs[order(column_dists_upper)]
      #
      # append
      for (k in 1:count_idxs) {
        #
        # append
        column_ints_rmdp2_ordered <- rbind(column_ints_rmdp2_ordered, column_ints_rmdp2[idxs_ordered[k], ])
      }
    }
    else {
      #
      # append
      column_ints_rmdp2_ordered <- rbind(column_ints_rmdp2_ordered, column_ints_rmdp2[idxs, ])
    }
  }
  #
  # initial
  strip_nodes2 <- c()
  #
  # initial
  idx_start <- 0
  idx_end <- 0
  #
  # insert
  for (i in 1:length(column_positions_unique)) {
    #
    # get index
    idx_start <- idx_end + 1
    idx_end <- column_positions_unique[i]
    #
    # get for position
    column_ints_rmdp2_ordered_i <- column_ints_rmdp2_ordered[which(column_positions==idx_end), ]
    #
    # insert
    strip_nodes2 <- rbind(strip_nodes2, strip_nodes[idx_start:idx_end, ], column_ints_rmdp2_ordered_i)
  }
  #
  # insert
  strip_nodes2 <- rbind(strip_nodes2, strip_nodes[(idx_end+1):count_node, ])
  #
  # debugging, check
  if (T) {
    #
    # check
    if (nrow(strip_nodes2) != (nrow(strip_nodes) + nrow(column_ints_rmdp2_ordered))) {
      #
      # warning
      warning("Error: add intersections of triangles in nodes failed.")
      #
      # pause at here
      return(NULL)
    }
  }
  #
  # get count
  count_node <- nrow(strip_nodes2)
  #
  # get length 2d
  for (i in 2:count_node) {
    #
    # get
    Lhrz <- 
      strip_nodes2[i-1, "Lhrz"] + 
      sqrt((strip_nodes2[i, "Cx"] - strip_nodes2[i-1, "Cx"])^2 + 
             (strip_nodes2[i, "Cy"] - strip_nodes2[i-1, "Cy"])^2)
    #
    # update
    strip_nodes2[i, "Lhrz"] <- Lhrz
  }
  #
  # get Czd (Cz of DEM)
  for (i in 1:count_node) {
    #
    # get
    strip_nodes2[i, "Czd"] <- raster::extract(pRasterDEM, strip_nodes2[i, c("Cx", "Cy")])
  }
  #
  # get length 3d (for Czd)
  for (i in 2:count_node) {
    #
    # get
    Lalld <- 
      strip_nodes2[i-1, "Lalld"] + 
      sqrt((strip_nodes2[i, "Cx"] - strip_nodes2[i-1, "Cx"])^2 + 
             (strip_nodes2[i, "Cy"] - strip_nodes2[i-1, "Cy"])^2 + 
             (strip_nodes2[i, "Czd"] - strip_nodes2[i-1, "Czd"])^2)
    #
    # update
    strip_nodes2[i, "Lalld"] <- Lalld
  }
  #
  # get Czp (Cz of plane fitted to pnts), only for non- anchors
  for (i in 1:count_node) {
    #
    # check
    if (strip_nodes2[i, "type"] != "Anchor: Initial group anchor" &&
        strip_nodes2[i, "type"] != "Anchor: Group center" &&
        strip_nodes2[i, "type"] != "Anchor: Inter-group center" &&
        strip_nodes2[i, "type"] != "Anchor: Distal group anchor") {
      #
      # get length horizontal
      Lhrz <- strip_nodes2[i, "Lhrz"]
      #
      # initial
      index_anchor_upper <- NA
      index_anchor_lower <- NA
      #
      # get index
      for (k in 1:count_node) {
        #
        # check
        if (strip_nodes2[k, "type"] == "Anchor: Initial group anchor" ||
            strip_nodes2[k, "type"] == "Anchor: Group center" ||
            strip_nodes2[k, "type"] == "Anchor: Inter-group center" ||
            strip_nodes2[k, "type"] == "Anchor: Distal group anchor") {
          #
          # check
          if (strip_nodes2[k, "Lhrz"] > Lhrz) {
            #
            index_anchor_lower <- k
            #
            break
          }
          #
          # get
          index_anchor_upper <- k
        }
      }
      #
      # get length ratio
      length_ratio <- 
        (Lhrz - strip_nodes2[index_anchor_upper, "Lhrz"]) /
        (strip_nodes2[index_anchor_lower, "Lhrz"] - strip_nodes2[index_anchor_upper, "Lhrz"])
      #
      # get Czp
      strip_nodes2[i, "Czp"] <- 
        strip_nodes2[index_anchor_upper, "Czp"] + 
        length_ratio *
        (strip_nodes2[index_anchor_lower, "Czp"] - strip_nodes2[index_anchor_upper, "Czp"])
    }
  }
  #
  # get length 3d (for Czp), for all nodes
  for (i in 2:count_node) {
    #
    # get
    Lallp <- 
      strip_nodes2[i-1, "Lallp"] + 
      sqrt((strip_nodes2[i, "Cx"] - strip_nodes2[i-1, "Cx"])^2 + 
             (strip_nodes2[i, "Cy"] - strip_nodes2[i-1, "Cy"])^2 + 
             (strip_nodes2[i, "Czp"] - strip_nodes2[i-1, "Czp"])^2)
    #
    # update
    strip_nodes2[i, "Lallp"] <- Lallp
  }
  #
  # owing to numerical error
  # ints are not exactly the same
  # in "remove duplicated", duplicated ints might be not removed?
  # do "remove duplicated" for end anchors
  if ((strip_nodes2[count_node, "Lhrz"] - strip_nodes2[count_node-1, "Lhrz"]) <= 10^-6) {
    #
    # remove
    strip_nodes2 <- rbind(strip_nodes2[1:(count_node-2), ], strip_nodes2[count_node, ])
  }
  if ((strip_nodes2[2, "Lhrz"] - strip_nodes2[1, "Lhrz"]) <= 10^-6) {
    #
    # remove
    strip_nodes2 <- rbind(strip_nodes2[1, ], strip_nodes2[3:count_node, ])
  }
  #
  # get count
  count_node <- nrow(strip_nodes2)
  #
  # debugging, save
  if (F) {
    #
    # save
    f_save_list2sp(list_triangle_sp, NULL, NA, "triangles")
    #
    # save
    f_save_list2sp(list_strip_sp, NULL, NA, "strips")
    #
    # save
    f_pnts_save_points(strip_nodes, NA, "strip_nodes")
    #
    # save
    f_pnts_save_points(column_ints, NA, "strip_ints")
    #
    # save
    f_pnts_save_points(column_ints_rmdp, NA, "strip_ints_rmdp")
    #
    # save
    f_pnts_save_points(strip_nodes2, NA, "strip_nodes2")
  }
  #
  # initial
  strip_nodes2[, "Czs"] <- c(NA)
  strip_nodes2[, "Lalls"] <- c(0)
  #
  # get Czs (Cz of surface)
  for (i in 1:count_node) {
    #
    # get node
    node <- strip_nodes2[i, ]
    #
    # node xy
    node_x <- node[1, "Cx"]
    node_y <- node[1, "Cy"]
    # node xy extnd
    node_x_extnd <- node_x
    node_y_extnd <- node_y
    #
    # still need to extend initial and distal group anchor?
    # although you also update their Cx and Cy
    # to the intersection with profile
    # and, do not extend only a little amount
    # as you already get intersection with profile
    # just use the middle point to locate triangle
    # it will not go to the next triangle
    if (T) {
      #
      # if Initial group anchor
      if (node[1, "type"] == "Anchor: Initial group anchor") {
        #
        # get
        node_x_extnd <- (strip_nodes2[i+1, "Cx"] + node_x_extnd) / 2
        node_y_extnd <- (strip_nodes2[i+1, "Cy"] + node_y_extnd) / 2
      }
      #
      # if Distal group anchor
      if (node[1, "type"] == "Anchor: Distal group anchor") {
        #
        # get
        node_x_extnd <- (strip_nodes2[i-1, "Cx"] + node_x_extnd) / 2
        node_y_extnd <- (strip_nodes2[i-1, "Cy"] + node_y_extnd) / 2
      }
    }
    #
    # initial
    pip <- 0
    #
    # get Czs
    for (k in 1:nrow(column_triangles)) {
      #
      # get index
      idx1 <- column_triangles[k, 1]
      idx2 <- column_triangles[k, 2]
      idx3 <- column_triangles[k, 3]
      #
      # get coords
      Cx1 <- column_coords[idx1, 1]
      Cy1 <- column_coords[idx1, 2]
      Cz1 <- column_coords[idx1, 3]
      # get coords
      Cx2 <- column_coords[idx2, 1]
      Cy2 <- column_coords[idx2, 2]
      Cz2 <- column_coords[idx2, 3]
      # get coords
      Cx3 <- column_coords[idx3, 1]
      Cy3 <- column_coords[idx3, 2]
      Cz3 <- column_coords[idx3, 3]
      #
      # get xyz
      tx <- rbind(Cx1, Cx2, Cx3, Cx1)
      ty <- rbind(Cy1, Cy2, Cy3, Cy1)
      tz <- rbind(Cz1, Cz2, Cz3, Cz1)
      #
      # point.in.polygon
      pip <- sp::point.in.polygon(node_x_extnd, node_y_extnd, tx, ty, mode.checked = FALSE)
      #
      # check
      if (pip != 0) {
        #
        # get normal
        nx <- (Cy2-Cy1)*(Cz3-Cz1) - (Cz2-Cz1)*(Cy3-Cy1)
        ny <- (Cz2-Cz1)*(Cx3-Cx1) - (Cx2-Cx1)*(Cz3-Cz1)
        nz <- (Cx2-Cx1)*(Cy3-Cy1) - (Cy2-Cy1)*(Cx3-Cx1)
        #
        # get Czs, nz won't be zero?
        Czs <- -(nx*(node_x-Cx1) + ny*(node_y-Cy1)) / nz + Cz1
        #
        # update
        strip_nodes2[i, "Czs"] <- Czs
        #
        # break
        break
      }
    }
    #
    # check
    if (pip == 0) {
      #
      # sometimes, profile will go through outside landslide polygon?
      # lead to intersections with triangles on the boundary of landslide polygon?
      #
      # get
      node_x_extnd1 <- (strip_nodes2[i+1, "Cx"] + node_x_extnd) / 2
      node_y_extnd1 <- (strip_nodes2[i+1, "Cy"] + node_y_extnd) / 2
      #
      # get
      node_x_extnd2 <- (strip_nodes2[i-1, "Cx"] + node_x_extnd) / 2
      node_y_extnd2 <- (strip_nodes2[i-1, "Cy"] + node_y_extnd) / 2
      #
      # get Czs
      for (k in 1:nrow(column_triangles)) {
        #
        # get index
        idx1 <- column_triangles[k, 1]
        idx2 <- column_triangles[k, 2]
        idx3 <- column_triangles[k, 3]
        #
        # get coords
        Cx1 <- column_coords[idx1, 1]
        Cy1 <- column_coords[idx1, 2]
        Cz1 <- column_coords[idx1, 3]
        # get coords
        Cx2 <- column_coords[idx2, 1]
        Cy2 <- column_coords[idx2, 2]
        Cz2 <- column_coords[idx2, 3]
        # get coords
        Cx3 <- column_coords[idx3, 1]
        Cy3 <- column_coords[idx3, 2]
        Cz3 <- column_coords[idx3, 3]
        #
        # get xyz
        tx <- rbind(Cx1, Cx2, Cx3, Cx1)
        ty <- rbind(Cy1, Cy2, Cy3, Cy1)
        tz <- rbind(Cz1, Cz2, Cz3, Cz1)
        #
        # point.in.polygon
        pip1 <- sp::point.in.polygon(node_x_extnd1, node_y_extnd1, tx, ty, mode.checked = FALSE)
        pip2 <- sp::point.in.polygon(node_x_extnd2, node_y_extnd2, tx, ty, mode.checked = FALSE)
        #
        # check
        if (pip1 != 0 || pip2 != 0) {
          #
          # get normal
          nx <- (Cy2-Cy1)*(Cz3-Cz1) - (Cz2-Cz1)*(Cy3-Cy1)
          ny <- (Cz2-Cz1)*(Cx3-Cx1) - (Cx2-Cx1)*(Cz3-Cz1)
          nz <- (Cx2-Cx1)*(Cy3-Cy1) - (Cy2-Cy1)*(Cx3-Cx1)
          #
          # get Czs, nz won't be zero?
          Czs <- -(nx*(node_x-Cx1) + ny*(node_y-Cy1)) / nz + Cz1
          #
          # update
          strip_nodes2[i, "Czs"] <- Czs
          #
          # break
          break
        }
      }
    }
    #
    # check
    if (is.na(strip_nodes2[i, "Czs"])) {
      #
      # save
      f_save_list2sp(list_strip_sp, NULL, NA, "strips")
      #
      # save
      f_save_list2sp(list_triangle_sp, NULL, NA, "strips_triangles")
      #
      # save
      f_pnts_save_points(strip_nodes2, NA, "strip_nodes2")
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
  # get length 3d (for Czs), update for all nodes
  for (i in 2:count_node) {
    #
    # get
    Lalls <- 
      strip_nodes2[i-1, "Lalls"] + 
      sqrt((strip_nodes2[i, "Cx"] - strip_nodes2[i-1, "Cx"])^2 + 
             (strip_nodes2[i, "Cy"] - strip_nodes2[i-1, "Cy"])^2 + 
             (strip_nodes2[i, "Czs"] - strip_nodes2[i-1, "Czs"])^2)
    #
    # update
    strip_nodes2[i, "Lalls"] <- Lalls
  }
  #
  # initial
  strip_nodes2$Cz <- c(0)
  strip_nodes2$Lall <- c(0)
  #
  # switch
  strip_nodes2 <- strip_nodes2[, c("Cx", "Cy", "Cz", "Czd", "Czp", "Czs", "type", "grp", "station", "Lhrz", "Lall", "Lalld", "Lallp", "Lalls")]
  #
  # set a choice here
  iElevation = 1
  #
  # use elevation of DEM
  if (iElevation == 1) {
    #
    # get, use d
    strip_nodes2$Cz <- strip_nodes2$Czd
    strip_nodes2$Lall <- strip_nodes2$Lalld
  }
  # use elevation of fitted plane
  else if (iElevation == 2) {
    #
    # get, use p
    strip_nodes2$Cz <- strip_nodes2$Czp
    strip_nodes2$Lall <- strip_nodes2$Lallp
  }
  # use elevation of TIN surface
  else if (iElevation == 3) {
    #
    # get, use s
    strip_nodes2$Cz <- strip_nodes2$Czs
    strip_nodes2$Lall <- strip_nodes2$Lalls
  }
  # use elevation of DEM
  else {
    #
    # get, use d
    strip_nodes2$Cz <- strip_nodes2$Czd
    strip_nodes2$Lall <- strip_nodes2$Lalld
  }
  #
  # initial length strip
  column_Lhrzs <- data.frame()
  column_Lalls <- data.frame()
  column_Lallss <- data.frame()
  #
  # get length strip
  for (i in 2:count_node) {
    #
    # check
    if (strip_nodes2[i, "type"] == "Station: Strip boundary") {
      #
      # append
      column_Lhrzs <- rbind(column_Lhrzs, strip_nodes2[i, "Lhrz"])
      column_Lalls <- rbind(column_Lalls, strip_nodes2[i, "Lall"])
      column_Lallss <- rbind(column_Lallss, strip_nodes2[i, "Lalls"])
    }
  }
  #
  # append
  column_Lhrzs <- rbind(column_Lhrzs, strip_nodes2[count_node, "Lhrz"])
  column_Lalls <- rbind(column_Lalls, strip_nodes2[count_node, "Lall"])
  column_Lallss <- rbind(column_Lallss, strip_nodes2[count_node, "Lalls"])
  #
  # get length strip
  for (i in count_strip:2) {
    #
    # cumulative to ...
    column_Lhrzs[i, 1] <- column_Lhrzs[i, 1] - column_Lhrzs[i-1, 1]
    column_Lalls[i, 1] <- column_Lalls[i, 1] - column_Lalls[i-1, 1]
    column_Lallss[i, 1] <- column_Lallss[i, 1] - column_Lallss[i-1, 1]
  }
  #
  # get width
  column_Whrzs <- column_Ahrzs / column_Lhrzs
  column_Walls <- column_Asrfs / column_Lallss
  # get Aall
  column_Aalls <- column_Lalls * column_Walls
  #
  # get sp attribute
  strip_df <- cbind(column_Lhrzs, column_Lalls, column_Ahrzs, column_Aalls, column_Whrzs, column_Walls, column_drops)
  #
  # colnames
  colnames(strip_df) <- c("Lhrz", "Lall", "Ahrz", "Aall", "Whrz", "Wall", "Drop")
  #
  # get
  triangle_df <- column_triangles_df
  #
  # return
  return(list(strip_df, strip_nodes2, list_triangle_sp, triangle_df))
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
        int.pts <- rgeos::gIntersection(sl1, sl2, byid = TRUE)
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
      return(NULL)
    }
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
f_intersection_curves2 <- function(pnts1, pnts2) {
  #
  # get coords
  x1 <- pnts1[1, 1]
  y1 <- pnts1[1, 2]
  # get coords
  x2 <- pnts1[2, 1]
  y2 <- pnts1[2, 2]
  #
  # get coords
  x3 <- pnts2[1, 1]
  y3 <- pnts2[1, 2]
  # get coords
  x4 <- pnts2[2, 1]
  y4 <- pnts2[2, 2]
  #
  # get denominator
  D = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
  #
  # check
  if (D == 0) {
    #
    # return
    return(NULL)
  }
  #
  # get ints
  intx <- 0
  intx <- intx + (x3 - x4) * (x1*y2 - y1*x2)
  intx <- intx - (x1 - x2) * (x3*y4 - y3*x4)
  intx <- intx / D
  # get ints
  inty <- 0
  inty <- inty + (y3 - y4) * (x1*y2 - y1*x2)
  inty <- inty - (y1 - y2) * (x3*y4 - y3*x4)
  inty <- inty / D
  # get ints
  int <- data.frame(intx, inty)
  colnames(int) <- c("Cx", "Cy")
  #
  # return
  return(int)
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
  # check
  if (ncol(column_anchors) == 3) {
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
  }
  #
  # check
  if (ncol(column_anchors) == 2) {
    #
    # get distances
    for (i in 1:(count_anchor-1)) {
      #
      # get distance
      distance <- 0
      distance <- distance + (column_anchors[i+1, 1] - column_anchors[i, 1])^2
      distance <- distance + (column_anchors[i+1, 2] - column_anchors[i, 2])^2
      distance <- sqrt(distance)
      #
      # update distances
      column_distances <- rbind(column_distances, distance)
    }
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
  paras_input <- data.frame(iMGPC = paras_input[1], iMGAD = paras_input[2], iMEAR = paras_input[3], iMEDA = paras_input[4], iEHAC = paras_input[5], iMSHL = paras_input[6])
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
    z1 = pnts[i-1, "Cz"]
    #
    x2 = pnts[i, "Cx"]
    y2 = pnts[i, "Cy"]
    z2 = pnts[i, "Cz"]
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
                    (pnts[count_pnts, "Cz"] - pnts[1, "Cz"])^2)
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