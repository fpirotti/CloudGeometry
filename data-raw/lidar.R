## code to prepare `lidar` dataset goes here
# lidar <- as.matrix( readr::read_delim("data-raw/fusa.txt",delim = "\t")
# download.file("https://github.com/LAStools/LAStools/raw/master/data/fusa.laz", "data-raw/fusa.laz")
# lidar <- ( lidR::readLAS("data-raw/fusa.laz", select = "xyz"))
# lidar  <- as.matrix( lidar @data )
# # inMatrix <- lidar2
#
# usethis::use_data(lidar, overwrite = TRUE)
