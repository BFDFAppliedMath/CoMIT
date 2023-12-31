% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-DB.R
\name{buildDBFromSchema}
\alias{buildDBFromSchema}
\title{Build a database from a given SQL schema file}
\usage{
buildDBFromSchema(schemaFile = NULL, dbName, force = FALSE, silently = FALSE)
}
\arguments{
\item{schemaFile}{a .sql schema file (if no value given, defaults to the base Covid-19 schema used in build_COVID_DB())}

\item{dbName}{name of database (with file path)}

\item{force}{if true, replace DB if already exists; if false, will not overwrite file}

\item{silently}{whether or not to print tables made}
}
\description{
This is a generic database build function that can create a database file from any .sql schema. If you use a schema other than the default one, later functions in the pipeline may not work.
}
