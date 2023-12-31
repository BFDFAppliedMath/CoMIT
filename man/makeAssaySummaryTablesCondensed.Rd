% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/table-assay_summary_condensed.R
\name{makeAssaySummaryTablesCondensed}
\alias{makeAssaySummaryTablesCondensed}
\title{Make Assay Summary Tables (Condensed Version)}
\usage{
makeAssaySummaryTablesCondensed(
  DB,
  variantFile,
  extraNotes = NULL,
  saveTable = FALSE,
  saveFolder = NULL,
  height = NULL,
  width = NULL,
  colWidth1 = NULL,
  colWidth2 = NULL,
  colWidth3 = NULL
)
}
\arguments{
\item{DB}{File path to subset database file. In the case of making summary tables, we use the
3M version from CoMIT, returned from \code{pull_subset_DB()}}

\item{variantFile}{csv file containing "Displayed_Lineage_Order" column for ordering row names/row groups in the GT table}

\item{extraNotes}{vector of strings that are extra source notes to be displayed at the bottom of the table}

\item{saveTable}{logical value indicating whether to save the tables to the specified saveFolder}

\item{saveFolder}{file path containing folder the finished tables will be saved in}

\item{height}{desired height of table .png images}

\item{width}{desired width of table .png images; if the width is not large enough, text will wrap or table will be cut off}

\item{colWidth1}{width of "5" column; if NULL, widths are automatically assigned to try to fit all text}

\item{colWidth2}{width of ">= 6" column; if NULL, widths are automatically assigned to try to fit all text}

\item{colWidth3}{width of "Sequences by Lineage" column; if NULL, widths are automatically assigned to try to fit all text}
}
\description{
Make 2 tables, a version with all mutations from the subset database, and a version with high-risk mutations; displays mutated sequence counts categorized by assay and number of assays affected on top and by variant on the left
}
