#
# ###################################################################### #
# ###################################################################### #
#                                                                        #
#                                ALPG                                    #
#                                                                        #
#     An R-Script for for automatic generation of landslide profile      #
#                                                                        #
#                             Version 1.0                                #
#                                                                        #
#                           July 7th, 2019                               #
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
# ######### commonly it is NOT SUGGESTED to change this script  #########
#
#
#
ALPG <- function(FileDEM, FileLandslides, MinGrpPnt = 3, MinGrpDist = 0, FilePrefix = "", OutputTemp = F) {
  #
  # check input
  library(dendextend)
  if (is.natural.number(MinGrpPnt) == F) {
    #
    # disp
    disp("Error: MinGrpPnt is not a natural number.")
    #
    # return
    return(NA)
  }
  #
  # check input
  if (MinGrpDist < 0) {
    #
    # disp
    disp("Error: MinGrpDist is negative.")
    #
    # return
    return(NA)
  }
  #
  # import vector layer
  library(rgdal)
  shp_lasld <- readOGR(FileLandslides)
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
  raster_dem <- raster(FileDEM)
  raster_CRS <- proj4string(raster_dem)
  #
  # get file
  file_out <- paste(FilePrefix, "_lasd_profile2d", sep = "", collapse = NULL)
  #
  # split
  for (i in 1:count_polygons) {
    #
    # get Polygons
    pPolygons <- list_polygons[[i]]
    #
    # split using on Polygons
    spldf_profile <- f_lasld_split(FileDEM, pPolygons, MinGrpPnt, MinGrpDist, FilePrefix, i, OutputTemp)
    #
    # update list
    if (is.null(spldf_profile) == F) {
      #
      # append list
      list_spldf_profile[[length(list_spldf_profile)+1]] <- spldf_profile
      #
      # list2shp
      f_save_list2spldf(list_spldf_profile, raster_CRS, file_out)
    }
  }
}
#
#
#
f_lasld_split <- function(FileDEM, pPolygons, MinGrpPnt, MinGrpDist, FilePrefix, FileID, OutputTemp) {
  #
  # library
  library(raster)
  library(rgdal)
  #
  # initial
  spldf_profile <- NULL
  #
  # open raster layer
  raster_dem <- raster(FileDEM)
  #
  # get CRS
  raster_CRS <- proj4string(raster_dem)
  #
  # get cell size
  raster_cellsize <- (res(raster_dem))[1]
  #
  # get SpatialPolygons
  pSpatialPolygons <- SpatialPolygons(list(pPolygons), proj4string = CRS(raster_CRS))
  #
  # rasterize polygon
  shp_lasld <- f_polygon_rasterize(raster_dem, pSpatialPolygons)
  # plot(shp_lasld, main = "polygon", axes = TRUE, border = "blue")
  #
  # set possible file_out (for the line of rasterized landslide plolygon)
  file_out <- paste(FilePrefix, as.character(FileID), "_lasld_pln2d", sep = "", collapse = NULL)
  #
  # get pnts
  pnts <- f_pnts_read(raster_dem, shp_lasld, file_out)
  #
  # if read pnts fail
  if (max(pnts[, "IDB"]) == 0) {
    #
    # error
    str_failed <- paste("Error: Read pnts failed for the ", as.character(FileID), "th landslide.", sep = "", collapse = NULL)
    disp(str_failed)
    #
    # save pnts
    file_out <- paste(FilePrefix, as.character(FileID), "_pnts.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, as.character(FileID), "_pnts", sep = "", collapse = NULL)
    f_pnts_save_shp(pnts, raster_CRS, file_out)
    #
    # return
    return(NA)
  }
  #
  # save the original pnts
  # if (OutputTemp) {
  if (T) {
    #
    # save pnts
    file_out <- paste(FilePrefix, as.character(FileID), "_pnts.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, as.character(FileID), "_pnts", sep = "", collapse = NULL)
    f_pnts_save_shp(pnts, raster_CRS, file_out)
  }
  #
  # initial pnts_split
  pnts_split <- pnts
  #
  # if normal is vertical (horizontal plane)
  if (pnts_split[1, "Dx"] == 0 && pnts_split[1, "Dy"] == 0 && pnts_split[1, "Dz"] == 0) { 
    #
    # error, if when reading pnts
    disp("Error: The first fitted plane for the landslide is horizontal.")
    #
    # set SMS, to be "LSH" (horizontal landslide)
    pnts_split[, "SMS"] <- "LSH"
    #
    # save pnts
    file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp.xlsx", sep = "", collapse = NULL)
    f_pnts_save_xls(pnts_split, file_out)
    #
    # save shp
    file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp", sep = "", collapse = NULL)
    f_pnts_save_shp(pnts, raster_CRS, file_out)
    #
    # return
    return(pnts_split)
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
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp", as.character(pnts_split[nrow(pnts_split), "grp"]), ".xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp", as.character(pnts_split[nrow(pnts_split), "grp"]), "", sep = "", collapse = NULL)
      f_pnts_save_shp(pnts_split, raster_CRS, file_out)
    }
    #
    # get centers
    pnts_split_anchors <- f_pnts_anchors(list_pnts, IDB_max)
    #
    # save for every split
    if (OutputTemp) {
      #
      # save centers
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp", as.character(pnts_split[nrow(pnts_split), "grp"]), "_anchors.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split_anchors, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp", as.character(pnts_split[nrow(pnts_split), "grp"]), "_anchors", sep = "", collapse = NULL)
      f_pnts_save_shp(pnts_split_anchors, raster_CRS, file_out)
      #
      # save shp pln
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp", as.character(pnts_split[nrow(pnts_split), "grp"]), "_anchors_pln2d", sep = "", collapse = NULL)
      f_pnts_save_shp_pln(pnts_split_anchors, raster_CRS, file_out, FileID)
    }
    #
    # break, if no group need split
    if (min(data.frame(list_grpm)) == 1) {
      # 
      # set SMS, to be "SMG", stop because all subgroup stop split
      pnts_split[, "SMS"] <- "SMG"
      #
      # save pnts
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp", sep = "", collapse = NULL)
      f_pnts_save_shp(pnts_split, raster_CRS, file_out)
      #
      # save centers
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp_anchors.xlsx", sep = "", collapse = NULL)
      f_pnts_save_xls(pnts_split_anchors, file_out)
      #
      # save shp
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp_anchors", sep = "", collapse = NULL)
      f_pnts_save_shp(pnts_split_anchors, raster_CRS, file_out)
      #
      # save shp pln
      file_out <- paste(FilePrefix, as.character(FileID), "_pnts_grp_anchors_pln2d", sep = "", collapse = NULL)
      spldf_profile <- f_pnts_save_shp_pln(pnts_split_anchors, raster_CRS, file_out, FileID)
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
    # for start pnts
    pnts_start <- list_pnts[[1]]
    pnts_start_bnd <- pnts_start[which(pnts_start[, "IDB"] != 0), ]
    #
    # for end pnts
    pnts_end <- list_pnts[[length(list_pnts)]]
    pnts_end_bnd <- pnts_end[which(pnts_end[, "IDB"] != 0), ]
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
        if ( isempty( which(pnts1[, "ID"] == pntu_i[1, "ID"]) ) ) {
          #
          # switch upper and lower
          pnts_i_upper <- pnts2
          pnts_i_lower <- pnts1
        }
        #
        # get clusters, check spatial contiguities
        list_clusters_upper <- f_pnts_clustering(pnts_i_upper, raster_cellsize * sqrt(2))
        list_clusters_lower <- f_pnts_clustering(pnts_i_lower, raster_cellsize * sqrt(2))
        #
        # save original
        pnts_i_upper_original <- pnts_i_upper
        pnts_i_lower_original <- pnts_i_lower
        #
        # allow for the upper subset
        if (length(list_clusters_upper) >= 2) {
          # #
          # # debugging
          # f_pnts_save_shp(pnts_i_upper, raster_CRS, "pnts_i_upper")
          # f_pnts_save_shp(pnts_i_lower, raster_CRS, "pnts_i_lower")
          #
          # update the upper and lower subsets
          for (k in 1:length(list_clusters_upper)) {
            #
            # get pnts
            pnts_cluster <- list_clusters_upper[[k]]
            #
            # assign
            if ( isempty( which(pnts_cluster[, "ID"] == pntu_i[1, "ID"]) ) ) {
              #
              # to lower
              pnts_i_lower <- rbind(pnts_i_lower, pnts_cluster)
            }
            else {
              #
              # to upper
              pnts_i_upper <- pnts_cluster
            }
          }
        }
        #
        # allow for the lower subset, also???
        # possibly both list_clusters_upper and list_clusters_lower >= 2, ???
        else if (length(list_clusters_lower) >= 2) {
          # #
          # # debugging
          # f_pnts_save_shp(pnts_i_upper, raster_CRS, "pnts_i_upper")
          # f_pnts_save_shp(pnts_i_lower, raster_CRS, "pnts_i_lower")
          #
          # initial
          ID_cluster_largest <- 1
          count_cluster_largest <- nrow(list_clusters_lower[[ID_cluster_largest]])
          #
          # get largest cluster, assign to lower
          for (k in 2:length(list_clusters_lower)) {
            #
            # get count
            count_cluster <- nrow(list_clusters_lower[[k]])
            #
            # update
            if (count_cluster > count_cluster_largest) {
              #
              # update
              ID_cluster_largest <- k
              count_cluster_largest <- count_cluster
            }
          }
          #
          # update the upper and lower subsets
          for (k in 1:length(list_clusters_lower)) {
            #
            # get pnts
            pnts_cluster <- list_clusters_lower[[k]]
            #
            # assign
            if (k == ID_cluster_largest) {
              #
              # to lower
              pnts_i_lower <- pnts_cluster
            }
            else {
              #
              # to upper
              pnts_i_upper <- rbind(pnts_i_upper, pnts_cluster)
            }
          }
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
        }
        #
        # debugging
        if ( (nrow(pnts_i_upper) + nrow(pnts_i_lower) ) != nrow(pnts_i)) {
          #
          # disp
          disp("Error: count of pnts in clustering is not consistent.")
          #
          # return
          return(NA)
        }
        #
        # get clusters, check spatial contiguities, again
        list_clusters_upper <- f_pnts_clustering(pnts_i_upper, raster_cellsize * sqrt(2))
        list_clusters_lower <- f_pnts_clustering(pnts_i_lower, raster_cellsize * sqrt(2))
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
          # f_pnts_save_shp(pnts_i_upper_original, raster_CRS, "pnts_i_upper_original")
          # f_pnts_save_shp(pnts_i_lower_original, raster_CRS, "pnts_i_lower_original")
          # #
          # # debugging
          # f_pnts_save_shp(pnts_i_upper, raster_CRS, "pnts_i_upper")
          # f_pnts_save_shp(pnts_i_lower, raster_CRS, "pnts_i_lower")
          # #
          # # disp
          # disp("Error: subsets can be clustered after combination.")
          # #
          # # return
          # return(NA)
        }
        #
        # count of pnt in sub- pnts is smaller than a threshold
        if (nrow(pnts_i_upper) < MinGrpPnt || nrow(pnts_i_lower) < MinGrpPnt) {
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
        #
        # check if sub- pnts have boundary points on two sides
        grpm_i_upper <- f_pnts_bndcorners_check(pnts_i_upper, IDB_max, pnts_start_bnd, pnts_end_bnd, boolEND_upper)
        grpm_i_lower <- f_pnts_bndcorners_check(pnts_i_lower, IDB_max, pnts_start_bnd, pnts_end_bnd, boolEND_lower)
        #
        # if sub- pnts do not have boundary points on two sides
        if (grpm_i_upper == 1 || grpm_i_lower == 1) {
          # #
          # # debugging
          # f_pnts_save_shp(pnts_i_upper, raster_CRS, "pnts_i_upper")
          # f_pnts_save_shp(pnts_i_lower, raster_CRS, "pnts_i_lower")
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
          f_pnts_save_shp(pnts_i_upper_original, raster_CRS, "pnts_i_upper_original")
          f_pnts_save_shp(pnts_i_lower_original, raster_CRS, "pnts_i_lower_original")
          #
          # debugging
          f_pnts_save_shp(pnts_i_upper, raster_CRS, "pnts_i_upper")
          f_pnts_save_shp(pnts_i_lower, raster_CRS, "pnts_i_lower")
          #
          # disp
          disp("Error: spatial in-contiguities are not expected here.")
          #
          # return
          return(NA)
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
        # get centers
        pnts_split_anchors <- f_pnts_anchors(list_pnts_new_anchor, IDB_max)
        #
        # get distances
        pnts_split_distances <- f_pnts_distances(pnts_split_anchors)
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
        if (pnts_split_distances_min <= MinGrpDist) {
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
  # library
  library(raster)
  library(rgdal)
  #
  # # plot raster_module
  # plot(raster_module, col = rev(terrain.colors(50)))
  #
  # # plot shapefile
  # # notice that you use add = T to add a layer on top of an existing plot in R.
  # plot(shp_pgn, main = "polygon", axes = TRUE, border = "blue")
  #
  # mask the raster using the polygon
  raster_module_mask <- mask(raster_module, shp_pgn)
  #
  # # plot
  # plot(raster_module_mask, main = "masked raster")
  # plot(shp_pgn, add = TRUE)
  #
  # set one
  raster_module_mask_set1 <- calc(raster_module_mask, function(x) x / x)
  #
  # # plot
  # plot(raster_module_mask_set1, main = "masked raster set1")
  # plot(shp_pgn, add = TRUE)
  #
  # raster to polygon
  raster_module_mask_set1_pgn <- rasterToPolygons(raster_module_mask_set1, fun = NULL, na.rm = TRUE, dissolve = TRUE)
  #
  # # plot shapefile
  # # notice that you use add = T to add a layer on top of an existing plot in R.
  # plot(raster_module_mask_set1_pgn, main = "polygon set1", axes = TRUE, border = "blue")
  # #
  # # write shapefile
  # writeOGR(raster_module_mask_set1_pgn, ".", "raster_module_mask_set1_pgn", driver = "ESRI Shapefile")
  #
  # return
  return(raster_module_mask_set1_pgn)
}
#
#
#
f_pnts_read <- function(raster_dem, shp_lasld, file_out = "lasld_pln") {
  #
  # library
  library(raster)
  library(rgdal)
  #
  # mask the raster using the vector extent
  raster_dem_mask <- mask(raster_dem, shp_lasld)
  #
  # get pnts from raster
  pnts_from_raster <- rasterToPoints(raster_dem_mask)
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
  # writeOGR(shp_lasld_pln, ".", "shp_lasld_pln", driver = "ESRI Shapefile")
  #
  # initial
  column_IDB <- c()
  #
  # maybe more than 1 (polygon) lines, difficult to handle
  if (length(shp_lasld_pln@lines) > 1) {
    #
    # disp
    disp("Error: shp_lasld_pln has more than 1 lines, read pnts failed.")
    #
    # save
    writeOGR(shp_lasld_pln, ".", file_out, overwrite_layer = T, driver = "ESRI Shapefile")
    #
    # get column_IDB
    column_IDB <- matrix(0, count_pnt, 1)
  }
  else if (length(shp_lasld_pln@lines[[1]]@Lines) > 1) {
    #
    # disp
    disp("Error: lines in shp_lasld_pln has more than 1 Lines, read pnts failed.")
    #
    # save
    writeOGR(shp_lasld_pln, ".", file_out, overwrite_layer = T, driver = "ESRI Shapefile")
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
    raster_cellsize <- (res(raster_dem))[1]
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
        # disp
        disp("Error: pnt_IDB is NA, read pnts failed.")
        #
        # save
        writeOGR(shp_lasld_pln, ".", file_out, overwrite_layer = T, driver = "ESRI Shapefile")
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
  library(pracma)
  pnts_strike <- cross(pnts_normal, vector_plumb)
  pnts_strike <- f_vector_unit(pnts_strike)
  #
  # vertical normal (horizontal plane)
  if (pnts_strike[1] ==0 && pnts_strike[2] == 0 && pnts_strike[3] ==0) { 
    #
    # waring
    disp("Warning: A fitted plane is horizontal.")
    #
    # use the dip of mother plane
    return(c(pnts[1, "Dx"], pnts[1, "Dy"], pnts[1, "Dz"])) 
  }
  #
  # get dip
  pnts_dip <- cross(pnts_strike, pnts_normal)
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
f_pnts_clustering <- function(pnts, distance) {
  #
  # initial list
  list_clusters <- list()
  #
  # get pnts
  pnts_remain <- pnts
  #
  # get count
  count_pnt_remain <- nrow(pnts_remain)
  #
  # clustering
  while (count_pnt_remain > 0) {
    #
    # if only 1 pnt
    if (count_pnt_remain == 1) {
      #
      # update
      list_clusters[[length(list_clusters) + 1]] <- pnts_remain
      #
      # break
      break
    }
    #
    # initial
    pnts_cluster <- pnts_remain[1, ]
    #
    # initial
    pnts_others <- pnts_remain[2:count_pnt_remain, ]
    #
    # initial
    boolCluster = TRUE
    #
    # clustering
    while (boolCluster) {
      #
      # get count
      count_pnt_others <- nrow(pnts_others)
      #
      # if others are empty
      if (count_pnt_others == 0) { break }
      #
      # clustering
      for (i in 1:count_pnt_others) {
        #
        # get count
        count_pnt_cluster <- nrow(pnts_cluster)
        #
        # get new member
        for (k in 1:count_pnt_cluster) {
          #
          # get distance
          distance_pp2 <- (pnts_others[i, "Cx"] - pnts_cluster[k, "Cx"])^2
          distance_pp2 <- (pnts_others[i, "Cy"] - pnts_cluster[k, "Cy"])^2 + distance_pp2
          distance_pp <- sqrt(distance_pp2)
          #
          # debugging
          if (is.na(distance_pp) || is.nan(distance_pp) || is.null(distance_pp) || isempty(distance_pp)) {
            #
            # return
            retur(NA)
          }
          #
          # replace
          if (distance_pp <= distance) {
            #
            # update cluster
            pnts_cluster <- rbind(pnts_cluster, pnts_others[i, ])
            #
            # delete
            ID_delete <- which(pnts_others[, "ID"] == pnts_others[i, "ID"])
            pnts_others <- pnts_others[-ID_delete, ]
            #
            # break
            break
          }
        }
        #
        # if have new member
        if (nrow(pnts_cluster) > count_pnt_cluster) {
          #
          # set true
          boolCluster = TRUE
          #
          # break
          break
        }
        #
        # if do not have new member, set FALSE
        boolCluster = FALSE
      }
    }
    #
    # update
    list_clusters[[length(list_clusters) + 1]] <- pnts_cluster
    #
    # update
    pnts_remain <- pnts_others
    #
    # get count
    count_pnt_remain <- nrow(pnts_remain)
  }
  #
  # initial
  count_pnt_list <- 0
  #
  # get count
  for (i in 1:length(list_clusters)) {
    #
    count_pnt_list <- count_pnt_list + nrow(list_clusters[[i]])
  }
  #
  # debugging
  if (count_pnt_list != nrow(pnts)) {
    #
    # disp
    disp("Error: count of pnts in clustering is not consistent.")
    #
    # return
    return(NA)
  }
  #
  # return
  return(list_clusters)
}
#
#
#
f_pnts_bndcorners_check <- function(pnts, IDB_max, pnts_start_bnd, pnts_end_bnd, boolEND) {
  #
  # get pnts_bnd
  pnts_bnd <- pnts[which(pnts[, "IDB"] != 0), ]
  #
  # get bndcorners
  bc_pnts_bnd <- f_pnts_bndcorners(pnts_bnd, IDB_max)
  #
  # check spatial impartiality
  if (is.na(bc_pnts_bnd)) {
    #
    # spatial impartiality
    return(1)
  }
  #
  # can not be 1???
  if (nrow(bc_pnts_bnd) == 1) {
    #
    # disp
    disp("Error: bndcorners count of initial or distal pnts is 1.")
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
      # disp
      disp("Error: bndcorners count of initial or distal pnts is 1.")
      return(NA)
    }
  }
  #
  # for debugging
  pnts_start_bnd_input <- pnts_start_bnd
  pnts_end_bnd_input <- pnts_end_bnd
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
    ID_start <- which(pnts_start_bnd[, "ID"] == pnts_bnd[i, "ID"])
    if (length(ID_start) != 0) {
      #
      pnts_start_bnd <- pnts_start_bnd[-ID_start, ]
    }
    #
    # remove duplicates
    ID_end <- which(pnts_end_bnd[, "ID"] == pnts_bnd[i, "ID"])
    if (length(ID_end) != 0) {
      #
      pnts_end_bnd <- pnts_end_bnd[-ID_end, ]
    }
  }
  #
  # get bndcorners of initial and distal pnts
  bc_pnts_start_bnd <- f_pnts_bndcorners(pnts_start_bnd, IDB_max)
  bc_pnts_end_bnd <- f_pnts_bndcorners(pnts_end_bnd, IDB_max)
  #
  # check end pnts after splitting
  if (nrow(bc_pnts_start_bnd) > 2 || nrow(bc_pnts_end_bnd) > 2) {
    #
    # not ok for end pnts
    return(3)
  }
  else if (nrow(bc_pnts_start_bnd) < 2 || nrow(bc_pnts_end_bnd) < 2) {
    #
    # disp
    disp("Error: bndcorners count of initial or distal pnts is 1.")
    return(NA)
  }
  #
  # get IDBs of initial and distal pnts
  IDB_start1 <- bc_pnts_start_bnd[1, "IDB"]
  IDB_start2 <- bc_pnts_start_bnd[2, "IDB"]
  IDB_end1 <- bc_pnts_end_bnd[1, "IDB"]
  IDB_end2 <- bc_pnts_end_bnd[2, "IDB"]
  #
  # circle from start1
  IDB_start1_start2 <- IDB_start2 - IDB_start1
  if (IDB_start1_start2 < 0) { IDB_start1_start2 <- IDB_start1_start2 + IDB_max }
  IDB_start1_end1 <- IDB_end1 - IDB_start1
  if (IDB_start1_end1 < 0) { IDB_start1_end1 <- IDB_start1_end1 + IDB_max }
  #
  # if get start2 first
  if (IDB_start1_start2 < IDB_start1_end1) {
    #
    # replace
    IDB_tmp <- IDB_start1
    IDB_start1 <- IDB_start2
    IDB_start2 <- IDB_tmp
  }
  #
  # circle from start1
  IDB_start1_end1 <- IDB_end1 - IDB_start1
  if (IDB_start1_end1 < 0) { IDB_start1_end1 <- IDB_start1_end1 + IDB_max }
  IDB_start1_end2 <- IDB_end2 - IDB_start1
  if (IDB_start1_end2 < 0) { IDB_start1_end2 <- IDB_start1_end2 + IDB_max }
  #
  # if get end2 first
  if (IDB_start1_end2 < IDB_start1_end1) {
    #
    # replace
    IDB_tmp <- IDB_end1
    IDB_end1 <- IDB_end2
    IDB_end2 <- IDB_tmp
  }
  #
  # circle from start1
  IDB_start1_start2 <- IDB_start2 - IDB_start1
  if (IDB_start1_start2 < 0) { IDB_start1_start2 <- IDB_start1_start2 + IDB_max }
  IDB_start1_end1 <- IDB_end1 - IDB_start1
  if (IDB_start1_end1 < 0) { IDB_start1_end1 <- IDB_start1_end1 + IDB_max }
  IDB_start1_end2 <- IDB_end2 - IDB_start1
  if (IDB_start1_end2 < 0) { IDB_start1_end2 <- IDB_start1_end2 + IDB_max }
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
    # circle from start1
    IDB_start1_bc <- IDB_bc - IDB_start1
    if (IDB_start1_bc < 0) { IDB_start1_bc <- IDB_start1_bc + IDB_max }
    #
    # check sides
    if (IDB_start1_bc < IDB_start1_end1) {
      #
      # 1st side
      count_bc_pnts_bnd_1 <- count_bc_pnts_bnd_1 + 1
    }
    else if (IDB_start1_bc > IDB_start1_end2 && IDB_start1_bc < IDB_start1_start2) {
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
    # disp
    disp("Error: bndcorners do not belong to either side.")
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
  # disp
  disp("Error: bndcorners count is not right.")
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
  if (count_pnt == 0) { return(NA) }
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
      # disp
      disp("IDB_step is NA")
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
    # disp
    disp("IDB_step is NA")
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
  if (is.na(bndcorners) || is.null(bndcorners)) {
    #
    # disp
    disp("bndcorners is NA or Null.")
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
    # get centerP
    centerP <- f_pnts_centerP_bbox(pnts)
    #
    # get xyz
    # cPx <- centerP[1]
    # cPy <- centerP[2]
    cPz <- centerP[3]
    #
    # get the point with maximum Cz
    lPntX <- pnts[which.max(pnts[, "Cz"]), "Cx"]
    lPntY <- pnts[which.max(pnts[, "Cz"]), "Cy"]
    #
    # get anchors, for start and end
    list_anchors <- f_pnts_anchor_end(pnts, c(lPntX, lPntY))
    anchor_end_se <- list_anchors[[1]]
    #
    # initial
    anchor_start <- c()
    anchor_end <- c()
    #
    # return two (anchor_end) or not
    if (nrow(anchor_end_se) == 2) {
      #
      # split start and end
      if (anchor_end_se[1, 3] >= cPz) {
        #
        # use the first one
        anchor_start <- anchor_end_se[1, ]
        anchor_end <- anchor_end_se[2, ]
      }
      else {
        #
        # use the first one
        anchor_start <- anchor_end_se[2, ]
        anchor_end <- anchor_end_se[1, ]
      }
    }
    else if (nrow(anchor_end_se) == 1) {
      #
      # it is possible that only return one anchor_end
      # disp("Warning: The initial group of landslide elevation points does not have start and end anchors.")
      #
      # split start and end
      if (anchor_end_se[1, 3] >= cPz) {
        #
        # use the first one
        anchor_start <- anchor_end_se[1, ]
        anchor_end <- c()
      }
      else {
        #
        # use the first one
        anchor_start <- c()
        anchor_end <- anchor_end_se[1, ]
      }
    }
    else {
      #
      # according to Cz
      anchor_start <- anchor_end_se[which.max(anchor_end_se[, "Cz"]), ]
      anchor_end <- anchor_end_se[which.min(anchor_end_se[, "Cz"]), ]
    }
    #
    # get column_anchors
    column_anchors <- rbind(anchor_start, centerP, anchor_end)
    column_anchors <- data.frame(column_anchors)
    #
    # colnames
    # colnames(column_anchors) <- c("Cx", "Cy", "zip")
    colnames(column_anchors) <- c("Cx", "Cy", "Cz")
    #
    # return
    return(column_anchors)
  }
  #
  # initial
  column_centerL <- c()
  #
  # for start pnts
  pnts_start <- list_pnts[[1]]
  #
  # get center of split "line"
  centerL_start <- f_pnts_centerL_end(pnts_start, IDB_max)
  column_centerL <- rbind(centerL_start, column_centerL)
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
        disp("Error: pnts do not have 4 bnd corners")
        #
        # return
        return(NA)
      }
    }
  }
  #
  # for end pnts
  pnts_end <- list_pnts[[count_pnts]]
  #
  # get center of split "line"
  centerL_end <- f_pnts_centerL_end(pnts_end, IDB_max)
  column_centerL <- rbind(column_centerL, centerL_end)
  # #
  # # output column_centerL for show
  # file_out <- paste("column_centerL_grp", as.character(pnts_end[nrow(pnts_end), "grp"]), ".xlsx", sep = "", collapse = NULL)
  # f_pnts_save_xls(column_centerL, file_out)
  #
  # initial
  column_centerG <- c()
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
  }
  #
  # initial
  column_anchors <- column_centerG
  #
  # colnames
  # colnames(column_anchors) <- c("Cx", "Cy", "zip")
  colnames(column_anchors) <- c("Cx", "Cy", "Cz")
  #
  # get anchors, for start
  list_anchors <- f_pnts_anchor_end(pnts_start, c(column_anchors[1, 1], column_anchors[1, 2]))
  anchor_end_start <- list_anchors[[1]]
  #
  # is possible that f_pnts_anchor_end return NA?
  if (is.na(anchor_end_start)) {
    #
    # return
    return(NA)
  }
  #
  # possible return more than one?
  if (nrow(anchor_end_start) >= 2) {
    # #
    # # return
    # return(NA)
    #
    # use the one with the maximum distance from anchor
    anchor_end_start <- list_anchors[[2]]
  }
  #
  # update anchors
  column_anchors <- rbind(anchor_end_start, list_anchors[[3]], column_anchors)
  #
  # get anchors, for end
  list_anchors <- f_pnts_anchor_end(pnts_end, c(column_anchors[nrow(column_anchors), 1], column_anchors[nrow(column_anchors), 2]))
  anchor_end_end <- list_anchors[[1]]
  #
  # is possible that f_pnts_anchor_end return NA?
  if (is.na(anchor_end_end)) {
    #
    # return
    return(NA)
  }
  #
  # possible return more than one?
  if (nrow(anchor_end_end) >= 2) {
    # #
    # # return
    # return(NA)
    #
    # use the one with the maximum distance from anchor
    anchor_end_start <- list_anchors[[2]]
  }
  #
  # update anchors
  column_anchors <- rbind(column_anchors, list_anchors[[3]], anchor_end_end)
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
    disp("Error: (end) pnts do not have 2 bnd corners")
    #
    # return
    return(NA)
  }
}
#
#
#
f_pnts_centerP_bbox <- function(pnts) {
  #
  # get xy
  xy <- cbind(pnts[, "Cx"], pnts[, "Cy"])
  #
  # getMinBBox
  # Minimum-area bounding box for a set of 2D-points
  # Calculates the vertices of the minimum-area, 
  # possibly oriented bounding box given a set of 2D-coordinates.
  library(shotGroups)
  bb <- getMinBBox(xy)
  #
  bb_pts <- bb$pts
  #
  # get center xy
  center_x <- mean(bb_pts[, 1])
  center_y <- mean(bb_pts[, 2])
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
#
#
#
f_pnts_anchor_end <- function(pnts, lPnt) {
  #
  # for the two end pnts (start and end)
  # and also
  # for the initial landslide elevation points
  #
  # get pnts_bnd
  pnts_bnd <- pnts[which(pnts[, "IDB"] != 0), ]
  #
  # get sorted
  pnts_bnd_ascending <- pnts_bnd[order(pnts_bnd[, "IDB"]), ]
  #
  # get count
  count_bnd <- nrow(pnts_bnd_ascending)
  #
  # initial
  ID_break <- 1
  #
  # get
  for (i in 2:count_bnd) {
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
  # sort
  pnts_bnd_ascending <- rbind(pnts_bnd_ascending[ID_break:count_bnd, ], pnts_bnd_ascending[1:ID_break-1, ])
  #
  # get pnts_curve
  pnts_curve_bnd <- pnts_bnd_ascending[, c("Cx", "Cy")]
  #
  # get cBBOX
  cBBOX <- f_pnts_centerP_bbox(pnts)
  #
  # get xyz
  cBBOXX <- cBBOX[1]
  cBBOXY <- cBBOX[2]
  #
  # get x and y
  lPntX <- lPnt[1]
  lPntY <- lPnt[2]
  #
  #
  #
  # DO NOT consider the situation
  # that cBBOX and lPnt are THE SAME point
  #
  #
  #
  # get line
  lX <- cBBOXX - lPntX
  lY <- cBBOXY - lPntY
  #
  # get min and max for x
  pnts_bnd_Xmin <- min(pnts_bnd[, "Cx"])
  pnts_bnd_Xmax <- max(pnts_bnd[, "Cx"])
  #
  # get min and max for y
  pnts_bnd_Xmin_Y <- lPntY + (pnts_bnd_Xmin - lPntX) * lY / lX
  pnts_bnd_Xmax_Y <- lPntY + (pnts_bnd_Xmax - lPntX) * lY / lX
  #
  # get pnts_curve
  pnts_curve_segment <- t(c(pnts_bnd_Xmin, pnts_bnd_Xmin_Y))
  pnts_curve_segment <- rbind(pnts_curve_segment, t(c(pnts_bnd_Xmax, pnts_bnd_Xmax_Y)))
  #
  # if vertical
  if (lX == 0) {
    #
    # get min and max for x
    pnts_bnd_Ymin <- min(pnts_bnd[, "Cy"])
    pnts_bnd_Ymax <- max(pnts_bnd[, "Cy"])
    #
    # get pnts_curve
    pnts_curve_segment <- t(c(lPntX, pnts_bnd_Ymin))
    pnts_curve_segment <- rbind(pnts_curve_segment, t(c(lPntX, pnts_bnd_Ymax)))
  }
  #
  # colnames
  colnames(pnts_curve_segment) <- c("Cx", "Cy")
  #
  # get data.frame
  pnts_curve_segment <- data.frame(pnts_curve_segment)
  #
  # get intersects
  ints <- f_intersection_curves(pnts_curve_bnd, pnts_curve_segment)
  #
  # if null, use bbox
  if (is.null(ints)) {
    #
    # disp
    disp("Error: no intscts gotten.")
    #
    # return
    return(NA)
  }
  #
  # get intersects x and y
  intsctsX <- data.frame(ints[, 1])
  intsctsY <- data.frame(ints[, 2])
  #
  # get count
  count_intsct <- nrow(intsctsX)
  #
  # debugging
  if (is.na(count_intsct) || is.nan(count_intsct) || is.null(count_intsct) || isempty(count_intsct)) {
    #
    # disp
    disp("Warning: empty intscts gotten.")
    #
    # return
    return(NA)
  }
  #
  # initial
  intsctsZ <- c()
  #
  # initial
  ID_MaxDist_intsct_lPnt <- 1
  MaxDist_intsct_lPnt <- -Inf
  #
  # get z (z in plane)
  for (i in 1:count_intsct) {
    #
    # get z (z in plane)
    intsctsZ <- rbind(intsctsZ, f_pnts_zip(pnts, intsctsX[i, 1], intsctsY[i, 1]))
    #
    # get dist
    Dist_intsct_lPnt <- sqrt((intsctsX[i, 1]-lPntX)^2 + (intsctsY[i, 1]-lPntY)^2)
    #
    # update
    if (Dist_intsct_lPnt > MaxDist_intsct_lPnt) {
      #
      # update
      MaxDist_intsct_lPnt <- Dist_intsct_lPnt
      ID_MaxDist_intsct_lPnt <- i
    }
  }
  #
  # get intersects
  intscts <- data.frame(cbind(intsctsX, intsctsY, intsctsZ))
  colnames(intscts) <- c("Cx", "Cy", "Cz")
  #
  # get x and y, anchor middle
  anchor_mid_x <- cBBOXX
  anchor_mid_y <- cBBOXY
  #
  # get dot
  cBBOX_lPnt_x <- lPntX - cBBOXX
  cBBOX_lPnt_y <- lPntY - cBBOXY
  cBBOX_ints_x <- intscts[ID_MaxDist_intsct_lPnt, "Cx"] - cBBOXX
  cBBOX_ints_y <- intscts[ID_MaxDist_intsct_lPnt, "Cy"] - cBBOXY
  # get dot
  dot_cBBOX <-  cBBOX_lPnt_x * cBBOX_ints_x + cBBOX_lPnt_y * cBBOX_ints_y
  #
  # cBBOX is outside the ints???
  # if (dot_cBBOX > 0) {
  # adjust anyway, or sometimes intsct and cBBOX will be too close
  # which in turn let the minimum group distance criterion not fair
  if (TRUE) {
    #
    # adjust
    anchor_mid_x <- (intscts[ID_MaxDist_intsct_lPnt, "Cx"] + lPntX) / 2
    anchor_mid_y <- (intscts[ID_MaxDist_intsct_lPnt, "Cy"] + lPntY) / 2
  }
  #
  # get z
  anchor_mid_z <- f_pnts_zip(pnts, anchor_mid_x, anchor_mid_y)
  #
  # get anchor middle
  anchor_mid <- data.frame(cbind(anchor_mid_x, anchor_mid_y, anchor_mid_z))
  colnames(anchor_mid) <- c("Cx", "Cy", "Cz")
  #
  # return
  return(list(intscts, intscts[ID_MaxDist_intsct_lPnt, ], anchor_mid))
}
#
#
#
f_intersection_curves <- function(pnts1, pnts2) {
  #
  # get x and y
  x1 <- pnts1[, 1]
  y1 <- pnts1[, 2]
  #
  # get x and y
  x2 <- pnts2[, 1]
  y2 <- pnts2[, 2]
  #
  # convert to a sp object (spatial lines)
  l1 <- Line(matrix(c(x1, y1), nc = 2, byrow = F))
  l2 <- Line(matrix(c(x2, y2), nc = 2, byrow = F))
  ll1 <- Lines(list(l1), ID = "1")
  ll2 <- Lines(list(l2), ID = "1")
  sl1 <- SpatialLines(list(ll1), proj4string = CRS("+init=epsg:4269"))
  sl2 <- SpatialLines(list(ll2), proj4string = CRS("+init=epsg:4269"))
  #
  # Calculate locations where spatial lines intersect
  int.pts <- gIntersection(sl1, sl2, byid = TRUE)
  #
  # debugging
  # if (is.na(int.pts) || is.nan(int.pts) || is.null(int.pts) || isempty(int.pts)) {
  if (is.null(int.pts) || isempty(int.pts)) {
    #
    # debugging
    f_pnts_save_shp(pnts1, CRS("+init=epsg:32643"), "pnts1")
    f_pnts_save_shp(pnts2, CRS("+init=epsg:32643"), "pnts2")
    #
    # disp
    disp("Warning: invalid intersection.")
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
    # get coords
    int.coords <- int.pts@pointobj@coords
  }
  #
  # other type?
  else {
    #
    #
    return(NA)
  }
  # #
  # # Plot line data and points of intersection
  # plot(x1, y1, type = "l")
  # lines(x2, y2, type = "l", col = "red")
  # points(int.coords[, 1], int.coords[, 2], pch = 20, col = "blue")
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
  if (is.na(CoordZ) || is.nan(CoordZ) || is.null(CoordZ) || isempty(CoordZ)) {
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
  library(xlsx)
  write.xlsx(pnts, file = file_out, sheetName = "pnts", row.names = FALSE, col.names = TRUE)
}
#
#
#
f_pnts_save_shp <- function(pnts, CRS_out, file_out) {
  #
  # library
  library(maptools)
  library(rgdal)
  library(sp)
  #
  # get coords
  PntsCoords <- pnts
  # coordinates(PntsCoords) = ~ Cx + Cy + Cz
  coordinates(PntsCoords) = ~ Cx + Cy
  #
  # get prj
  proj4string(PntsCoords) = CRS_out
  #
  # get SpatialPointsDataFrame
  PntsCoords_spdf <- SpatialPointsDataFrame(PntsCoords, pnts)
  #
  # write
  writeOGR(PntsCoords_spdf, dsn = ".", layer = file_out, overwrite = T, driver="ESRI Shapefile")
  # #
  # # write 3D, will produce extra Cx, Cy and Cz columns in shapefile???
  # # even, "Warning message: writePointsShape is deprecated; use rgdal::writeOGR or sf::st_write"
  # writePointsShape(PntsCoords_spdf, file_out, factor2char = TRUE, max_nchar = 254)
  # #
  # # indeed, writePointsShape() does not write the projection file
  # # but, using function showWKT from rgdal, you can also create one like that
  # cat(showWKT(proj4string(PntsCoords_spdf)), file = paste(file_out, ".prj", sep = "", collapse = NULL))
  #
  # return
  return(PntsCoords_spdf)
}
#
#
#
f_pnts_save_shp_pln <- function(pnts, CRS_out, file_out, LineID = 1) {
  #
  # library
  library(maptools)
  library(rgdal)
  library(sp)
  #
  # get spldf_data
  # no R library supports 3d line now?
  spldf_data <- cbind(pnts[, "Cx"], pnts[, "Cy"], pnts[, "Cz"])
  spldf_data <- cbind(pnts[, "Cx"], pnts[, "Cy"])
  rownames(spldf_data) <- c(1:nrow(spldf_data))
  #
  # get split line
  spldf <- f_pnts2spldf(spldf_data, LineID)
  #
  # get prj
  proj4string(spldf) = CRS_out
  #
  # write
  writeOGR(spldf, dsn = ".", layer = file_out, overwrite = T, driver="ESRI Shapefile")
  #
  # return
  return(spldf)
}
#
#
#
f_pnts2spldf <- function(pnts, IDin) {
  #
  # get 
  pnts_Line <- Line(pnts)
  pnts_Lines <- Lines(list(pnts_Line), ID = "a")
  pnts_SpatialLines <- SpatialLines(list(pnts_Lines))
  #
  # get
  pnts_df <- data.frame(ID = IDin)
  rownames(pnts_df) <- c("a")
  #
  # get
  spldf <- SpatialLinesDataFrame(pnts_SpatialLines, data = pnts_df)
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
  data_ID <- data.frame()
  #
  # get stack
  for (i in 1:count_spldf) {
    #
    # get spldf
    spldf_i <- list_spldf[[i]]
    #
    # get Lines
    pLine <- Line(spldf_i@lines[[1]]@Lines[[1]]@coords)
    pLines <- Lines(list(pLine), ID = i)
    #
    # update list
    list_SpatialLines[[i]] <- pLines
    #
    # get ID
    data_ID <- rbind(data_ID, spldf_i@data$ID)
  }
  #
  # get SpatialLines
  pSpatialLines = SpatialLines(list_SpatialLines)
  #
  # set names
  colnames(data_ID) <- "ID"
  rownames(data_ID) <- c(1:nrow(data_ID))
  #
  # get spldf
  spldf <- SpatialLinesDataFrame(pSpatialLines, data = data_ID)
  # plot(spldf, col = c("red"))
  #
  # get prj
  proj4string(spldf) = CRS_out
  #
  # write
  writeOGR(spldf, dsn = ".", layer = file_out, overwrite = T, driver="ESRI Shapefile")
}