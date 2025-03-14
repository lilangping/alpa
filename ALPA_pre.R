#
# ###################################################################### #
# ###################################################################### #
#                                                                        #
#                                ALPA                                    #
#                                                                        #
#       An R-Script for for automatic analysis of landslide path         #
#                                                                        #
#                             Version 4.4                                #
#                                                                        #
#                          February 24th, 2024                           #
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
#
# 2. Please download and install Rtools,
#    following official instructions of Rtools.
#
# 3. Important: compatibility of version.
#    This version of ALPA works well on the Windows PC of the authors,
#    with R-4.2.2, Rstudio 2022.07.2+576, Rtools42.
#    And, all packages were updated to the latest version on NOV. 21th, 2022.
#
# 4. Please execute the script below,
#    to install all required packages.
#    And, be aware that the installation will take quite a while.
#
# 5. Please try a manual installation,
#    if you have difficulties in the installation of any package.
#
# 6. If error messages tell that some packages are missing or can't find,
#    try to manually install those packages.
#
# 7. Please add the paths of proj of "rgdal" and "sf" libraries,
#    to the system environment variable "PROJ_LIB",
#    for avoiding the possible warning message "Cannot find proj.db".
#    The paths depend on where you install the "rgdal" and "sf" packages.
#
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
# 1. The main function ALPA has the following inputs:
#    -> FileDEM:              input DEM, tif raster file is recommended.
#    -> FileLandslides:       input landslide polygons, shapefile vector file.
#    -> FileParameters:       input file for several parameters, xlsx file.
#    -> CfgEndAnchorsInt:     input configuration for end anchors (coordinates), initial.
#    -> CfgEndAnchorsDit:     input configuration for end anchors (coordinates), distal.
#    -> CfgEndStripsInt:      input configuration for end strips (directions), initial.
#    -> CfgEndStripsDit:      input configuration for end strips (directions), distal.
#    -> OutputHrcl:           whether output hierarchical results/files, T or F (default).
#
# 2. The parameters file (FileParameters) must indicate the following 6 inputs:
#    -> FilePrefix:           prefix for output files, usually the name of landslide.
#    -> MinGrpPntCnt:         minimum group point count.
#    -> MinGrpAcrDst:         minimum group anchor distance.
#    -> MinStpHrzLen:         minimum strip horizontal length.
#    -> MinEndRtoInt:         minimum ratio for end group, initial.
#    -> MinEndRtoDit:         minimum ratio for end group, distal.
#    ** Assume there are n landslides to be treated, 
#       then, the parameters file should be a xlsx table with (n+1) row and 6 column,
#       and, the first row is for the headings of the 6 parameter inputs,
#       and, the parameters in the ith row will applied to the ith landslide.
#       Users can just set FileParameters as the filename, e.g., "FileParameters.xlsx".
#    ** Please note: The parameters file is not obligatory.
#       If no parameters file is specified, set FileParameters as 0, 
#       and all landslides will apply the same default parameters as follows.
#    -> FilePrefix:           "lsdxxxxxx"
#    -> MinGrpPntCnt:         3
#    -> MinGrpAcrDst:         0
#    -> MinStpHrzLen:         30
#    -> MinEndRtoInt:         0.0000
#    -> MinEndRtoDit:         0.0000
#    ** Please note: "lsdxxxxxx" indicates the order of a landslide in FileLandslides.
#
# 3. The configuration for end anchors, indicates 
#    either the coordinates of end anchors, 
#    or the algorithm used to obtain end anchors.
#    ** The recommendation is using point shapefiles to indicate the coordinates of end anchors.
#       A point shapefile for initial end groups, and another one for distal end groups.
#       The point shapefile must contains n point features,
#       each of them indicate the location (coordinates) of the end anchor of a landslide.
#       Users can just set CfgEndAnchorsInt and CfgEndAnchorsDit as the filenames, 
#       e.g., "FileEndAnchorsInt.shp" and "FileEndAnchorsDit.shp".
#    ** If users do not want to define point shapefiles, 
#       predefined algorithms can be used to obtain end anchors.
#       The code 1, 2 and 3 indicate "mbb", "quad" and "even" algorithms, respectively,
#       which stand for using "minimum bounding box", "quadrilateral" and "even division" 
#       to split end groups, respectively.
#       The "mbb" algorithm is the most fast one, but is also expected to be the most rough one.
#       The "quad" algorithm is expected to be the most sophisticated one, but is also very time-consuming.
#       The speed of the "even" algorithm is between those of "mbb" and "quad".
#    ** Users can just set CfgEndAnchorsInt and CfgEndAnchorsDit as 1, 2 or 3, 
#       to indicate the algorithms used for initial and distal end groups.
#       In this way, all landslides will apply the same algorithm.
#       The recommendation is to use "mbb" (1) as a start, which is also the default choice.
#
# 4. The configuration for end strips, indicates 
#    either the directions of end strips, 
#    or the algorithm used to obtain end strips.
#    ** The recommendation is using polyline shapefiles to indicate the directions of end strips.
#       A polyline shapefile for initial end groups, and another one for distal end groups.
#       The polyline shapefile must contains n polyline features,
#       each of them indicate the orientation (direction) of the end strip of a landslide.
#       Users can just set CfgEndStripsInt and CfgEndStripsDit as the filenames, 
#       e.g., "FileEndStripsInt.shp" and "FileEndStripsDit.shp".
#    ** If users do not want to define polyline shapefiles, 
#       predefined algorithms can be used to obtain end strips.
#       The code 1, 2 and 3 indicate "mbb", "quad" and "prll" algorithms, respectively,
#       which stand for using "minimum bounding box", "quadrilateral" and "parallel edge" 
#       to determine the directions of end strips, respectively.
#       The "prll" algorithm means the end strip direction is "parallel" to the direction of the end group edge,
#       which is the most fast one, but is also expected to be the most rough one.
#       The "quad" algorithm is expected to be the most sophisticated one, but is also very time-consuming.
#       The speed of the "mbb" algorithm is between those of "prll" and "quad".
#    ** Users can just set CfgEndStripsInt and CfgEndStripsDit as 1, 2 or 3, 
#       to indicate the algorithms used for initial and distal end groups.
#       In this way, all landslides will apply the same algorithm.
#       The recommendation is to use "mbb" (1) as a start, which is also the default choice.
#
# 5. Please make sure, in the input landslide, parameter, and configuration files, 
#    the orders (of landslides/features) in the xlsx and shapefile dbf tables are consistent.
#
# 6. Please make sure, all input GIS data (DEM file and shapefiles), 
#    have the same coordinate system.
#
# 7. Please make sure, the resolution of the input DEM, 
#    is integer.
#
# 8. It is recommended not to output hierarchical results/files, 
#    i.e. to set OutputHrcl as "F" (also the default choice),
#    especially when processing an inventory consisting of many landslides.
#    If OutputHrcl is "T", only hierarchical sub- groups and anchors will be output.
#
# 9. The output points and profiles are two dimensional (2D).
#    3D (z-aware) points and profiles can be produced based on the output excel files,
#    using software like ArcGIS.
#
# 10. Before applying the main function,
#    please execute the script below,
#    to set working directory and source the main R script.
# 
# The script for setting working directory.
# The following directory is just an example.
# Please change it to the working directory on your computer.
# Please put all the R scripts (ALPA_*.R) and data files in the working directory.
setwd("D:\\ALPA")
#
# The script for sourcing the main R script.
rm(list = ls())
source('ALPA_main.R')
#
#
#
# === EXAMPLE ===
#
# DATA for the example are AVAILABLE in this repository.
# Four landslides in the Central Asian area are used for illustrations.
# 
# Assume that we have the following inputs of DEM and landslides:
# -> FileDEM:              "In_DEM.tif"
# -> FileLandslides:       "In00_Landslides.shp"
# We can execute the following command, by applying pure default settings:
ALPA("In_DEM.tif", "In00_Landslides.shp")
# The above command is equivalent to the following ones:
ALPA("In_DEM.tif", "In00_Landslides.shp", 0)
ALPA("In_DEM.tif", "In00_Landslides.shp", 0, 1, 1, 1, 1)
ALPA("In_DEM.tif", "In00_Landslides.shp", 0, 1, 1, 1, 1, F)
# Please note: the above command is mostly an trail for new users.
# Usually, it is not recommended to apply pure default settings.
# In realistic applications, in the first trail, users are encouraged 
# to directly apply files to define parameters and configurations for individual landslides. 
#
# If we are ready to apply different parameters and configurations for different landslides,
# and, assume we have the following input files:
# -> FileLandslides:       "In01_Landslides.shp"
# -> FileParameters:       "In01_Parameters.xlsx".
# -> CfgEndAnchorsInt:     "In01_EndAnchorsInt.shp".
# -> CfgEndAnchorsDit:     "In01_EndAnchorsDit.shp".
# -> CfgEndStripsInt:      "In01_EndStripsInt.shp".
# -> CfgEndStripsDit:      "In01_EndStripsDit.shp".
# We can execute the following command:
ALPA("In_DEM.tif", "In01_Landslides.shp", "In01_Parameters.xlsx", "In01_EndAnchorsInt.shp", "In01_EndAnchorsDit.shp", "In01_EndStripsInt.shp", "In01_EndStripsDit.shp")
# Please note: "In01_Landslides.shp" and "In00_Landslides.shp" are identical, 
# and, in "In01_Parameters.xlsx", default parameters are actually used for all landslides,
# except the FilePrefix, which is defined by the names of landslides.
# The trick is to set appropriate EndAnchors and EndStrips configurations for individual landslides in this step, 
# and, to adjust MinEndRto parameters for individual landslides (with unsatisfactory results), in the next step.
#
# We will check the results obtained by the above command, 
# detect landslides with unsatisfactory results (anchors, path and strips), and make adjustments.
# Usually, appropriate EndAnchors and EndStrips configurations have been already defined in the above step,
# and, in this step, MinEndRto parameters will be adjusted for landslides with unsatisfactory results.
# We found that, landslides "Kyrgyz_247R" and "Kyrgyz_248R" had obtained satisfying results in the above step,
# so, in this step, we will adjust MinEndRto parameters for "Kyrgyz_314aR" and "Kyrgyz_314R".
# After adjustments, assume we have the following input files:
# -> FileLandslides:       "In02_Landslides.shp"
# -> FileParameters:       "In02_Parameters.xlsx".
# -> CfgEndAnchorsInt:     "In02_EndAnchorsInt.shp".
# -> CfgEndAnchorsDit:     "In02_EndAnchorsDit.shp".
# -> CfgEndStripsInt:      "In02_EndStripsInt.shp".
# -> CfgEndStripsDit:      "In02_EndStripsDit.shp".
# We can execute the following command:
ALPA("In_DEM.tif", "In02_Landslides.shp", "In02_Parameters.xlsx", "In02_EndAnchorsInt.shp", "In02_EndAnchorsDit.shp", "In02_EndStripsInt.shp", "In02_EndStripsDit.shp")
# Please note: in this example, the "In02_xxx" input files are only for "Kyrgyz_314aR" and "Kyrgyz_314R".
# In this way, "Kyrgyz_247R" and "Kyrgyz_248R" with satisfying results will not be treated again, 
# which can save time especially when the processing of some landslides costs much time.
# If time is not a problem, users can also adjust MinEndRto parameters 
# for landslides with unsatisfactory results, directly in "In01_Parameters.xlsx".
#
#
#
# === EXPERIENCES (ON SETTINGS) ===
# 
# 01. Please be aware that there are no a set of universally appropriate settings,
#     for all landslide inventories, or for all landslides in an inventory.
#     Settings should be defined according to every individual landslide cases.
#
# 02. According to personal experiences obtained from experiments on limited landslide inventories,
#     the following recommendation is a good start.
#
# 03. Set MinGrpPntCnt and MinGrpAcrDst to be 3 and 0, respectively.
#     That is to say, no constraints from group point count and group anchor distance.
#
# 04. The spatial resolution of the used DEM is a reasonable value for MinStpHrzLen.
#     MinStpHrzLen controls number of strips, and will not influence the generation of path
#     Smaller MSHL values will in general produce a larger number of strips.
#     Be aware that the smallest count of strip is 2, no matter how large MinStpHrzLen is.
#
# 05. Set MinEndRtoInt and MinEndRtoDit to be 0 and 0, respectively.
#     That is to say, no constraints from point count ratio of end groups.
#
# 06. Define a parameters file, using the above suggested parameters, 
#     and, using the names of landslides as FilePrefix.
#
# 07. Define EndAnchors and EndStrips configurations shapefiles, 
#     according to your experts' knowledge.
#
# 08. TIPS: if the path (sub- grouping) of some landslides are not satisfying,
#     try to re-process those landslides using respective appropriate settings.
#     Usually, the adjustment in this step is for MinEndRto parameters.
#     By checking end groups, adjust MinEndRtoInt and MinEndRtoDit,
#     which will be a little bit larger than, 
#     the ratio of the points in the to-be-merged end sub- groups to all points.
#
# 09. It is not recommended to use predefined algorithms to obtain end anchors and strips.
#     If predefined algorithms ("mmb", "quad", "even", and "prll") are used, 
#     users can try different combinations of predefined algorithms to get satisfying results,
#     but, the results will be less controllable.
#
# 10. Normally, we do not need to use MinGrpAcrDst, and it is a last resort.
#     Do not use MinGrpAcrDst (minimum group anchor distance), Whenever possible.
#     The use of MinGrpAcrDst will reduce the resolution of profile (count of anchors).
#     
# 11. The input DEM is expected to represent the sliding surface.
#     If the input DEM represents the post-sliding or pre-sliding surface,
#     which will be most likely in applications, a profile will be still generated;
#     but, be aware that, elevation of profile will not represent elevation of sliding surface.
#
# 12. Some problems of running ALPA might be owing to the environment of user computer.
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
# Version 4.0 (November 4th, 2023)
# 1. Update the algorithm for obtaining end anchors and end strips.
# 2. Fix a bug in splitting end groups with small point counts.
#
# Version 4.1 (December 18th, 2023)
# 1. Fix a bug in handling manually-defined end anchors out end groups.
# 2. Disable outputting console to file using sink() because it does not work appropriately.
#
# Version 4.2 (January 15th, 2024)
# 1. Fix a bug in the situation that no sub-grouping is valid.
# 2. Update the path generation algorithm for the situation that no sub-grouping is valid.
#
# Version 4.3 (January 23th, 2024)
# 1. Fix a bug in handling the entire landslide points (no sub-groups).
#
# Version 4.4 (February 24th, 2024)
# 1. Fix a bug in merging (replacing) initial and distal groups.
#
#
#