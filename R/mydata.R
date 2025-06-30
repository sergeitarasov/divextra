#' Example YAML Data: My Data NEW
#'
#' A description of the data stored in the YAML file.
#'
#' @details This YAML file contains ...
#' @format yaml
#' @examples
#' file_path <- system.file("extdata", "geosse3_td_test.yml", package = "divextra")
#' data <- yaml::read_yaml(file_path)
#' @name yaml_files
NULL


#"2-regions-td-test"


#' Geosse tree
#'
#' A description of the data stored in the YAML file.
#'
#' @details This YAML file contains ...
#' @format yaml
#' @examples
#' file_path <- system.file("extdata", "geosse_tb_tree.rds", package = "divextra")
#' phy <- readRDS(file_path)
#' @name geosse_tree
NULL


#
# file_path <- system.file("extdata", "geosse3_td_test.yml", package = "divextra")
# par.categories.td <- read_yaml_pars_td(file_path)
# formula.td <- make_constraints_sse_td(par.categories.td)


#' Toy Tree Data
#'
#' A phylogenetic tree with five tips
#'
#' @format phylo object
#' @source Generated using diversitree simulation functions
#'
#' @examples
#' data(phy5)
#' print(phy5)
#'
"phy5"
