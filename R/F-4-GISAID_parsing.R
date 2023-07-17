
#' Parse metadata from a Dataframe
#'
#' @param accClade Dataframe containing only rows with relevant data. The first row will be used.
#'
#' @return List with keys: ("clade", "lineage", "subDate")
#'          if one of the keys is not present it will be NA
#'
getMetaData <- function(accClade) {
  metaD <- list()

  if (s$GISAID_CLADE %in% colnames(accClade)) {
    metaD[["clade"]] <- accClade[[s$GISAID_CLADE]][1]
  } else {
    metaD[["clade"]] <- NA
  }

  if (s$PANGO_LINEAGE %in% colnames(accClade)) {
    metaD[["lineage"]] <- accClade[[s$PANGO_LINEAGE]][1]
  } else {
    metaD[["lineage"]] <- NA
  }

  return(metaD)
}


#' Parse GISAID FASTA line for metadata and check for sequences with no accession number,
#' the wrong host, and/or invalid collection dates
#'
#'
#' @param line A string containing the initial line of a FASTA item (line after >).
#'             This must be pulled from the GISAID database
#'
#' @return List with keys: ("accNum", "location", "colDate", "continue") containing
#'         the Accession ID, location, collection date, and error status receptively.
#'         If a host other than human is found or an accession number is not identified
#'         the "continue" field will be marked false and a "Fail_Type" key will be added
#'         containing failure mode details
#'
parseFastaLine <- function(line) {
  metaInfo <- list()
  tryCatch(
    expr = {
      metaInfo[["continue"]] <- TRUE

      seqData <- strsplit(line, split = s$STRING_SPLIT_EXPRESSION)
      metaInfo[["location"]] <- seqData[[1]][s$LOCATION_INDEX]
      metaInfo[["accNum"]] <- stringr::str_match(line, s$ACCESSION_REG_EXP)[, 2]
      metaInfo[["colDate"]] <- stringr::str_extract(line, s$COL_DATE_REG_EXP)
      metaInfo[["colDay"]] <- NA

      if (is.na(metaInfo[["accNum"]])) {
        print(paste("Found NA Accession:", line))
        metaInfo[["continue"]] <- FALSE
        metaInfo[["Fail_Type"]] <- "NA Accession"
        return(metaInfo)
      }

      # Check for human
      if (metaInfo[["location"]] %in% s$NON_HUMAN_CHECK) {
        print(paste("Host Error:", metaInfo[["location"]]))
        metaInfo[["continue"]] <- FALSE
        metaInfo[["Fail_Type"]] <- "Host Error"
        return(metaInfo)
      }

      # Deal with date
      if (nchar(metaInfo[["colDate"]]) > 7) {
        cDate <- lubridate::parse_date_time(c(metaInfo[["colDate"]]), orders = "ymd")
        metaInfo[["colDate"]] <- toString(format(cDate, "%Y-%m"))
        metaInfo[["colDay"]] <- toString(format(cDate, "%d"))
      } else if (nchar(metaInfo[["colDate"]]) > 5) {
        cDate <- lubridate::parse_date_time(c(metaInfo[["colDate"]]), orders = "ym")
        metaInfo[["colDate"]] <- toString(format(cDate, "%Y-%m"))
        metaInfo[["colDay"]] <- NA
      } else {
        metaInfo[["continue"]] <- FALSE
        metaInfo[["Fail_Type"]] <- "Date is too short"
        return(metaInfo)
      }
    },
    error = function(e) {
      metaInfo[["continue"]] <- FALSE
      metaInfo[["Fail_Type"]] <- paste("Parsing Error:", line)
      print(metaInfo[["Fail_Type"]])
    },
    warning = function(w) {
      metaInfo[["continue"]] <- FALSE
      metaInfo[["Fail_Type"]] <- paste("Parsing Warning:", line)
      print(metaInfo[["Fail_Type"]])
    }
  )

  return(metaInfo)
}
