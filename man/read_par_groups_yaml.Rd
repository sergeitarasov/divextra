% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-yaml.R
\name{read_par_groups_yaml}
\alias{read_par_groups_yaml}
\title{Read Parameter Groups from YAML File}
\usage{
read_par_groups_yaml(yaml_file)
}
\arguments{
\item{yaml_file}{Path to YAML file containing parameter groups}
}
\value{
List of parameter group lists, each named according to YAML structure
}
\description{
Reads multiple named lists of parameter groups from a YAML file.
Each list in the YAML file becomes a separate named element in the returned list.
}
\examples{
\dontrun{
# Read parameter groups from YAML file
groups <- read_par_groups_yaml("parameter_groups.yml")

# Access specific group lists
groups$parameter_groups    # First group list
groups$parameter_groups2   # Second group list
}
}
