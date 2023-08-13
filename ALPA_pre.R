#
# ###################################################################### #
# ###################################################################### #
#                                                                        #
#                                ALPA                                    #
#                                                                        #
#      An R-Script for for automatic analysis of landslide profile       #
#                                                                        #
#                             Version 3.1                                #
#                                                                        #
#                           July 20th, 2023                              #
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
#    following official instructions of Rtools.
# 3. Important: compatibility of version.
#    This version of ALPA works well on the Windows PC of the authors,
#    with R-4.2.2, Rstudio 2022.07.2+576, Rtools42.
#    And, all packages were updated to the latest version on NOV. 21th, 2022.
# 4. Please execute the script below,
#    to install all required packages.
#    And, be aware that the installation will take quite a while.
# 5. Please try a manual installation,
#    if you have difficulties in the installation of any package.
# 6. If error messages tell that some packages are missing or can't find,
#    try to manually install those packages.
# 7. Please add the paths of proj of "rgdal" and "sf" libraries,
#    to the system environment variable "PROJ_LIB",
#    for avoiding the possible warning message "Cannot find proj.db".
#    The paths depend on where you install the "rgdal" and "sf" packages.
# 8. Restart RStudio.
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
# 3. Please make sure the resolution of the input DEM is integer.
#
# 4. The main function ALPA has the following inputs:
#    -> FileDEM:              input DEM file, tif raster file is recommended.
#    -> FileLandslides:       input polygon of landslide, shapefile vector file.
#    -> MinGrpPntCnt:         minimum group point count.
#    -> MinGrpAcrDst:         minimum group anchor distance.
#    -> MinEndRtoInt:         minimum ratio for end group, initial.
#    -> MinEndRtoDit:         minimum ratio for end group, distal.
#    -> EndAgrCfgInt:         algorithm configuration for handling end group (anchor and strip), initial.
#    -> EndAgrCfgDit:         algorithm configuration for handling end group (anchor and strip), distal.
#    -> MinStpHrzLen:         minimum strip horizontal length.
#    -> FilePrefix:           prefix for output files.
#    -> OutputHrcl:           whether output hierarchical results/files, T or F (default).
#
# 5. The output points and profiles are two dimensional (2D).
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
# -> MinEndRtoInt:          0.00
# -> MinEndRtoDit:          0.00
# -> EndAgrCfgInt:          1
# -> EndAgrCfgDit:          1
# -> MinStpHrzLen:          30
# -> FilePrefix:            "case_c003_d000_ri0.00_rd0.00_ei01_ed01_s030"
# 
# Recommend not to output hierarchical results/files, 
# especially when processing an inventory consisting of many landslides,
# and to execute:
ALPA("dem.tif", "lsd.shp", 3, 0, 0.00, 0.00, 1, 1, 30, "case_c003_d000_ri0.00_rd0.00_ei01_ed01_s030")
# 
# If you process an individual landslide,
# and want to inspect the hierarchical sub- grouping steps,
# you can try to output hierarchical results/files, and to execute:
ALPA("dem.tif", "lsd.shp", 3, 0, 0.00, 0.00, 1, 1, 30, "case_c003_d000_ri0.00_rd0.00_ei01_ed01_s030", T)
# Please NOTE: only hierarchical sub- groups and anchors will be output,
# and, will be very time consuming, if set EndAgrCfg to be 1 or 4.
#
#
#
# === ExPERIENCES (ON PARAMETER SELECTION) ===
# 
# 01. Please be aware that there are no a set of universally appropriate parameters,
#     for all landslide inventories, or for all landslides in an inventory.
#     Parameters should be selected according to every individual landslide cases.
# 02. According to personal experiences obtained from experiments on limited landslide inventories,
#     the following recommendation is a good start.
# 03. Set MinGrpPntCnt and MinGrpAcrDst to be 3 and 0, respectively.
#     That is to say, no constraints from group point count and group anchor distance.
# 04. Select MinEndRtoInt and MinEndRtoDit to be 0 and 0, respectively.
#     That is to say, no constraints from point count ratio of end groups.
# 05. EndAgrCfgInt and EndAgrCfgDit, define algorithm for generating end anchors and strips, and can be 1 ~ 6.
#     1: quad + quad for anchor and strip respectively, slow.
#     2:  mbb + quad for anchor and strip respectively, slow.
#     3: even + quad for anchor and strip respectively, slow.
#     4: quad + mbb  for anchor and strip respectively, slow.
#     5:  mbb + mbb  for anchor and strip respectively, 1st fast.
#     6: even + mbb  for anchor and strip respectively, 2nd fast.
#     "quad", "mbb", and "even" stand for quadrilateral, minimum bounding box, and evenly splitting, respectively.
#     The recommendation for EndAgrCfg is 1.
# 06. TIPS: if the profile (sub- grouping) of some landslides are not satisfied,
#     try to re-process those landslides individually.
#     First, by checking initial and distal end groups respectively, determine MinEndRtoInt and MinEndRtoDit,
#     which will be a little bit larger than the ratio of the to-be-merged end sub- groups to all points.
#     Then, determine EndAgrCfgInt and EndAgrCfgDit respectively by optimizing end anchors and strips.
#     Try 1, 2, and 3 to see which algorithm (quad, mbb, or even) is the best for generating end anchors.
#     If end strips are not satisfied, try correspondingly 4, 5, and 6 to see if mmb is better for generating end strips.
#     That is to say, if quad yields the best end anchor, try 4; if mbb, try 5; if even, try 6.
#     EndAgrCfg having quad might be very time consuming.
#     Therefore, if the direction of end strip is not the concern,
#     the recommended order of EndAgrCfg is 5, 6, and then 4.
# 07. For some landslide cases, adjusting MinEndRtoInt, MinEndRtoDit, EndAgrCfgInt and EndAgrCfgDit is not adequate.
#     Try to adjust the MinGrpAcrDst (minimum group anchor distance).
#     Use MinGrpAcrDst as a last resort. Whenever possible, don't use MinGrpAcrDst.
#     The use of MinGrpAcrDst will reduce the resolution of profile (count of anchors).
# 08. The spatial resolution of the used DEM is a reasonable value for MinStpHrzLen.
#     MinStpHrzLen controls number of strips, and will not influence the generation of profile.
#     Smaller MSHL values will in general produce a larger number of strips.
#     Be aware that the smallest count of strip is 2, no matter how large MinStpHrzLen is.
# 09. The input DEM is expected to represent the sliding surface.
#     If the input DEM represents the post-sliding or pre-sliding surface,
#     which will be most likely in applications, a profile will be still generated;
#     but, be aware that, elevation of profile will not represent elevation of sliding surface.
# 10. Some problems of running ALPA might be owing to the environment of user computer.
#     Try to firstly solve them by searching solutions regarding R and RStudio.
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
# Version 2.7 (November 21th, 2022)
# 1. Fix generating profile in three-dimensional space.
# 2. Fix bugs owing to compatibility of version.
# 3. Temporarily disable output of hierarchical results/files.
#
# Version 3.0 (June 09th, 2023)
# 1. Remove the stop criterion based on deviation angle in sub- grouping end groups.
# 2. Add a stop criterion based on point count ratio in sub- grouping end groups.
# 3. Fix a bug in checking spatial contiguity
# 4. Fix a bug in fitting plane for sub- group pnts.
# 5. Fix bugs in differentiating right and left boundaries.
# 6. Consider the situation that anchors or profile are outside landslide polygon.
# 7. Update the algorithm configuration for handling the initial and distal end group.
# 8. Display warning and error messages.
#
# Version 3.1 (July 20th, 2023)
# 1. Update the algorithm for minimizing point count ratio of end groups.
# 2. Fix a bug in handling the end group with the quad algorithm.
#
#
#