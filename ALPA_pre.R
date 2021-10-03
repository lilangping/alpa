#
# ###################################################################### #
# ###################################################################### #
#                                                                        #
#                                ALPA                                    #
#                                                                        #
#      An R-Script for for automatic analysis of landslide profile       #
#                                                                        #
#                             Version 2.6                                #
#                                                                        #
#                         September 26th, 2021                           #
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
# === INSTALLATION ===
#
# 1. Please make sure that,
#    you have installed the latest and compatible versions of R and RStudio.
# 2. Please download and install Rtools, 
#    and also put the location of the Rtools make utilities on the PATH,
#    following official instructions.
# 3. Please execute the script below,
#    to install all required packages.
#    And, be aware that the installation will take quite a while.
# 4. Please try a manual installation,
#    if you have difficulties in the installation of any package.
# 5. If error messages tell that some packages are missing or can't find,
#    try to manually install those packages.
# 6. Please add the paths of proj of "rgdal" and "sf" libraries,
#    to the system environment variable "PROJ_LIB",
#    for avoiding the possible warning message "Cannot find proj.db".
#    The paths depend on where you install the "rgdal" and "sf" packages.
# 7. Restart RStudio.
#
# The script for installing all required packages
if (!require("dendextend"))   install.packages("dendextend", dependencies = TRUE)
if (!require("pracma"))   install.packages("pracma", dependencies = TRUE)
if (!require("xlsx"))   install.packages("xlsx", dependencies = TRUE)
if (!require("shotGroups"))   install.packages("shotGroups", dependencies = TRUE)
if (!require("circular"))   install.packages("circular", dependencies = TRUE)
if (!require("lwgeom"))   install.packages("lwgeom", dependencies = TRUE)
if (!require("decido"))   install.packages("decido", dependencies = TRUE)
if (!require("raster"))   install.packages("raster", dependencies = TRUE)
if (!require("rgdal"))   install.packages("rgdal", dependencies = TRUE)
if (!require("rgeos"))   install.packages("rgeos", dependencies = TRUE)
if (!require("sp"))   install.packages("sp", dependencies = TRUE)
if (!require("sf"))   install.packages("sf", dependencies = TRUE)
#
#
#
# === INSTRUCTIONS ===
#
# 1. Before applying the main function:
#    Please execute the script below,
#    to set working directory and source the main R script.
#
# 2. Please make sure all data have the same coordinate system.
#
# 3. The main function ALPA has the following inputs:
#    -> FileDEM:              input DEM file, tif raster file is recommended.
#    -> FileLandslides:       input polygon of landslide, shapefile vector file.
#    -> MinGrpPntCnt:         minimum group point count.
#    -> MinGrpAcrDst:         minimum group anchor distance.
#    -> MinEndAptRto:         minimum aspect ratio for end group (before split).
#    -> MaxEndDvtAgl:         maximum deviation angle for end group (after split).
#    -> EndHdlAgrCfg:         algorithm configuration for handling end group (anchor and direction).
#    -> MinStpHrzLen:         minimum strip horizontal length.
#    -> FilePrefix:           prefix for output files.
#    -> OutputTemp:           whether output temporary results/files, T or F (default).
#
# 4. The output points and profiles are two dimensional (2D).
#    3D (z-aware) points and profiles can be produced based on the output excel files,
#    using software like ArcGIS.
# 
# The script for setting working directory
# The following directory is just an example.
# Please change it to the working directory on your computer.
# Please put all the R scripts (ALPA_*.R) and data files in the working directory.
setwd("D:\\ALPA")
#
# The script for sourcing the main R script
rm(list = ls())
source('ALPA_main.R')
#
#
#
# === EXAMPLES ===
# 
# Assume that the following inputs will be adopted:
# -> FileDEM:               "dem.tif"
# -> FileLandslides:        "lsd.shp"
# -> MinGrpPntCnt:          3
# -> MinGrpAcrDst:          0
# -> MinEndAptRto:          1.5
# -> MaxEndDvtAgl:          15
# -> EndHdlAgrCfg:          1
# -> MinStpHrzLen:          30
# -> FilePrefix:            "case_c003_d000_r1.5_a015_e001_s030"
# 
# If you process an inventory consisting of many landslides,
# recommend not to output temporal results/files, and execute:
ALPA("dem.tif", "lsd.shp", 3, 0, 1.5, 15, 1, 30, "case_c003_d000_r1.5_a015_e001_s030")
# 
# If you process an individual landslide,
# you can try to output temporal results/files, and execute:
ALPA("dem.tif", "lsd.shp", 3, 0, 1.5, 15, 1, 30, "case_c003_d000_r1.5_a015_e001_s030", T)
#
#
#
# === ExPERIENCES (ON PARAMETER SELECTION) ===
# 
# 1. Please be aware that there are no a set of universally appropriate parameters,
#    for all landslide inventories, or for all landslides in an inventory.
#    Parameters should be selected according to every individual landslide cases.
# 2. According to personal experiences obtained from experiments on limited landslide inventories,
#    the following recommendation is a good start.
# 3. Set MinGrpPntCnt and MinGrpAcrDst to be 3 and 0, respectively.
#    That is to say, no constraints from group point count and group anchor distance.
# 4. Select MinEndAptRto and MaxEndDvtAgl applicable to your case study.
#    0 and 180 mean no constraints from deviation angle for end group.
#    Personal experiences suggest that 1.5 and 15 are good for most landslide cases.
# 5. For some landslide cases, adjusting MinEndAptRto and MaxEndDvtAgl is not adequate.
#    Try to adjust the MinGrpAcrDst (minimum group anchor distance).
# 6. EndHdlAgrCfg can be 1 ~ 12.
#    The recommendation for EndHdlAgrCfg is 1.
#    If the initial and distal anchors of profiles of some landslides are not satisfied,
#    try to re-process those landslides individually using an EndHdlAgrCfg value of 7 and 11.
#    If the directions of strips of end groups for some landslides are not satisfied, 
#    try to re-process those landslides individually,
#    using an EndHdlAgrCfg value of 2 (if profile is generated using an EndHdlAgrCfg value of 1),
#    or using 10 and 8 (if profile is generated using an EndHdlAgrCfg value of 7 and 11).
#    Other choices of EndHdlAgrCfg are not commonly used, and might be very time consuming.
# 7. The spatial resolution of the used DEM is a reasonable value for MinStpHrzLen.
#    MinStpHrzLen controls number of strips, and will not influence the generation of profile.
#    Smaller MSHL values will in general produce a larger number of strips.
#    Be aware that the smallest count of strip is 2, no matter how large MinStpHrzLen is.
# 8. The input DEM is expected to represent the sliding surface.
#    If the input DEM represents the post-sliding or pre-sliding surface,
#    which will be most likely in applications, a profile will be still generated;
#    but, be aware that, elevation of profile will not represent elevation of sliding surface.
# 9. Some problems of running ALPA 2.0 might be owing to the environment of user computer.
#    Try to firstly solve them by searching solutions regarding R and RStudio.
#
#
#
# === (major) UPDATES ===
# 
# Version 1.1 (January 12th, 2021)
# 1. Improve the efficiency of raster cropping (clipping) significantly.
# 2. Fix a bug in handling anchors of end groups.
#
# Version 2.0 (February 3rd, 2021)
# 1. Add the module for measuring landslide width and aspect ratio.
# 2. Add a criterion for stopping the segmentation of initial and distal groups.
# 3. Improve the approach for obtaining anchors for initial and distal groups.
# 4. Improve the efficiency of checking spatial contiguity significantly.
# 5. Format names of output files further.
# 6. Log time for processing.
#
# Version 2.1 (May 5th, 2021)
# 1. Improve the module for measuring landslide width and aspect ratio.
# 2. Improve the approach for checking spatial impartiality.
# 3. Improve the algorithm for generating anchors.
#
# Version 2.2 (May 12th, 2021)
# 1. Improve the algorithm for measuring overall length.
# 2. Prepare for categorizing landslide longitudinal shape.
#
# Version 2.3 (August 14th, 2021)
# 1. Improve the algorithm for generating strips of end group.
# 2. Improve the algorithm for generating initial and distal anchors.
#
# Version 2.4 (August 19th, 2021)
# 1. Formalize the algorithm configuration for handling end group.
#
# Version 2.5 (September 15th, 2021)
# 1. Consider the situation that DEM for the sliding surface is available.
#
# Version 2.6 (September 26th, 2021)
# 1. Improve the algorithm for generating profile in three-dimensional space.
#
#
#
