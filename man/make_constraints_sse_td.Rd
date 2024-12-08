% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-yaml.R
\name{make_constraints_sse_td}
\alias{make_constraints_sse_td}
\title{Create parameter constraints for episodic ClaSSE}
\usage{
make_constraints_sse_td(par.categories.td)
}
\arguments{
\item{par.categories.td}{a named list of parameters, for example, from read_yaml_pars_td()}
}
\value{
list of formulas
}
\description{
Create parameter constraints for episodic ClaSSE
}
\examples{
file_path <- system.file("extdata", "geosse3_td_maxpar.yml", package = "divextra")
par.categories.td <- read_yaml_pars_td(file_path)
print(par.categories.td)
formula.td <- make_constraints_sse_td(par.categories.td)
print(formula.td)
}
