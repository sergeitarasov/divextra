# Functions to set up R package

# usethis::use_roxygen_md()
# devtools::document()

# usethis::use_package("diversitree")
# usethis::use_package("yaml")
# usethis::use_package("ape")
# usethis::use_package("stringr")
# usethis::use_pipe()

# utils::getFromNamespace("make.initial.conditions.classe", "diversitree")

# devtools::load_all()

# devtools::check()
# devtools::build()

# usethis::use_gpl_license(version = 3, include_future = TRUE)

#---- Data
# usethis::use_data(mydataset, demographics, overwrite = T)
#
# ## Not run:
# x <- 1:10
# y <- 1:100
#
# use_data(x, y) # For external use
# use_data(x, y, internal = TRUE) # For internal use


#-------- Vignette

# usethis::use_vignette("Episodic-GeoSSE")

# devtools::build_vignettes()

#------ Testing
# report_p(-1)

# usethis::use_test("report_p")
# devtools::test()


#---- from ontoFAST


# #utils::globalVariables("shiny_in", package="ontoFAST")
# utils::globalVariables("ontofast", package="ontoFAST")
#
# usethis::use_package("shiny", type = "Depends")
# usethis::use_package("ontologyIndex", type = "Depends")
#
#
# # devtools::clean_dll()

# # devtools::build_vignettes()

# # devtools::document()
# # devtools::build()
# # devtools::install(build_vignettes = TRUE)
# # devtools::check()

