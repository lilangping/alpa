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
#                          February 3th, 2020                            #
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
# clear
rm(list = ls())
#
# install all required packages, if not available.
if (!require("dendextend"))   install.packages("dendextend", dependencies = TRUE)
if (!require("rgdal"))   install.packages("rgdal", dependencies = TRUE)
if (!require("raster"))   install.packages("raster", dependencies = TRUE)
if (!require("pracma"))   install.packages("pracma", dependencies = TRUE)
if (!require("dbscan"))   install.packages("dbscan", dependencies = TRUE)
if (!require("shotGroups"))   install.packages("shotGroups", dependencies = TRUE)
if (!require("xlsx"))   install.packages("xlsx", dependencies = TRUE)
if (!require("maptools"))   install.packages("maptools", dependencies = TRUE)
if (!require("sp"))   install.packages("sp", dependencies = TRUE)
if (!require("rgeos"))   install.packages("rgeos", dependencies = TRUE)
# 
# set working directory
# please put all the R scripts and data files in the working directory
# the following directory is just an example
# please change it to the working directory on your computer
setwd("D:\\ALPA")
#
# source the main R script
source('ALPA_main.R')
#
#
#
# INSTRUCTIONS
#
# 0. Please first execute the above command lines, before applying the main function
#
# 1. Please make sure all data have the same coordinate system
#
# 2. The main function ALPA has four input parameters
#    -> FileDEM: input DEM file, tif raster file is recommended
#    -> FileLandslides: input landslide polygon, shapefile vector file
#    -> MinGrpPntCnt: minimum group point count
#    -> MinGrpAncDist: minimum group achor distance
#    -> FilePrefix: prefix for output files
#    -> OutputTemp: whether output temporary results/files, T or F (default)
#
# 3. The output points and profiles are two dimensional (2D)
#    3D (z-aware) points and profiles can be produced using software like ArcGIS
#    based on the output excel files
#
#
#
# EXAMPLES OF APPLICATION
# 
# If you do not want to output temporal results/files (recommended)
# ALPA("lasd_dem.tif", "lasd_polygon.shp", 3, 0, "lasd")
# 
# If you want to output temporal results/files
# ALPA("lasd_dem.tif", "lasd_polygon.shp", 3, 0, "lasd", T)
#
#
#
