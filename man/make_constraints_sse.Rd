% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-yaml.R
\name{make_constraints_sse}
\alias{make_constraints_sse}
\title{Create parameter constraints for time-homogeneous ClaSSE}
\usage{
make_constraints_sse(par.categories)
}
\arguments{
\item{par.categoriesa}{a named list of parameters, for example, from read_yaml_pars_td()}
}
\value{
list of formulas
}
\description{
Create parameter constraints for time-homogeneous ClaSSE
}
\examples{
file_path <- system.file("extdata", "geosse3_ti_maxpar.yml", package = "divextra")
par.categories <- read_yaml_pars(file_path)
print(par.categories)
pars <- make_constraints_sse(par.categories)
print(pars)
}
