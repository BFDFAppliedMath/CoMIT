
clean_lineage_col <- function() {
  query = c(glue::glue("ALTER TABLE {s$SEQ_INFO_TABLE} ADD Cleaned_Lineage TEXT;"),
            glue::glue("UPDATE {s$SEQ_INFO_TABLE} SET Cleaned_Lineage = SUBSTR({s$LINEAGE}, 1, INSTR({s$LINEAGE}, ' ') -1);"),
            glue::glue("UPDATE {s$SEQ_INFO_TABLE} SET Cleaned_Lineage = NULLIF(Cleaned_Lineage, '');"),
            glue::glue("UPDATE {s$SEQ_INFO_TABLE} SET {s$LINEAGE} = Cleaned_Lineage WHERE Cleaned_Lineage IS NOT NULL;"),
            glue::glue("ALTER TABLE {s$SEQ_INFO_TABLE} DROP COLUMN Cleaned_Lineage;"))

  return(paste(query, collapse = " "))

}


build_who_des_sql <- function(linDesDF) {
  query <- c(glue::glue("UPDATE {s$SEQ_INFO_TABLE} SET {s$GROUP} = CASE"))
  if (nrow(linDesDF) > 0) {
    searchAll <- c(
      glue::glue("WHEN ({s$LINEAGE} IN (SELECT DISTINCT {s$LINEAGE}
                 FROM {s$SEQ_INFO_TABLE} WHERE {s$LINEAGE} like '"),
      "%')) THEN '", "'"
    )
    searchOne <- c(glue::glue("WHEN ({s$LINEAGE} = '"), "') THEN '", "'")
    for (row in 1:nrow(linDesDF)) {
      if (linDesDF$`Check_Descendants?`[row]) {
        searchType <- searchAll
      } else {
        searchType <- searchOne
      }
      lin <- linDesDF$Lineage[row]
      des <- linDesDF$WHO_Designation[row]
      query <- c(query, paste0(searchType[1], lin, searchType[2], des, searchType[3]))
    }
    query <- c(query, "WHEN (Lineage IS NULL) THEN 'Unassigned'", "ELSE 'Non-Omicron' END;")
  } else {
    print("No rows")
  }
  return(paste(query, collapse = " "))
}


build_lin_subgroup_sql <- function(linDesDF) {
  query <- c(glue::glue("UPDATE {s$SEQ_INFO_TABLE} SET {s$SUBGROUP} = CASE"))
  if (nrow(linDesDF) > 0) {
    searchAll <- c(
      glue::glue("WHEN ({s$LINEAGE} IN (SELECT DISTINCT {s$LINEAGE}
                 FROM {s$SEQ_INFO_TABLE} WHERE {s$LINEAGE} like '"),
      "%')) THEN '", "'"
    )
    searchOne <- c(glue::glue("WHEN ({s$LINEAGE} = '"), "') THEN '", "'")
    for (row in 1:nrow(linDesDF)) {
      if (linDesDF$`Check_Descendants?`[row]) {
        searchType <- searchAll
      } else {
        searchType <- searchOne
      }
      lin <- linDesDF$Lineage[row]
      subgroup <- linDesDF$Subgroup[row]
      query <- c(query, paste0(searchType[1], lin, searchType[2], subgroup, searchType[3]))
    }
    query <- c(query, "WHEN (Lineage IS NULL) THEN 'Unassigned'", "ELSE 'Non-Omicron' END;")
  } else {
    print("No rows")
  }
  return(paste(query, collapse = " "))
}


build_display_lin_sql <- function(linDesDF) {
  query <- c(glue::glue("UPDATE {s$SEQ_INFO_TABLE} SET {s$DISPLAY_LIN} = CASE"))

  if (nrow(linDesDF) > 0) {
    searchAll <- c(
      glue::glue("WHEN ({s$LINEAGE} IN (SELECT DISTINCT {s$LINEAGE}
                 FROM {s$SEQ_INFO_TABLE} WHERE {s$LINEAGE} like '"),
      "%')) THEN '",
      "'"
    )
    searchOne <- c(
      glue::glue("WHEN ({s$LINEAGE} = '"),
      "') THEN '",
      "'"
    )
    for (row in 1:nrow(linDesDF)) { # for every row in lineage designation file
      if (linDesDF$`Check_Descendants?`[row]) {
        searchType <- searchAll
      } else {
        searchType <- searchOne
      }
      lin <- linDesDF$Lineage[row] # lineage value at this row
      display_lin <- linDesDF$Displayed_Lineage[row]

      query <- c(query, paste0(searchType[1], lin, searchType[2], display_lin, searchType[3]))
    }
    query <- c(query, "WHEN (Lineage IS NULL) THEN 'Unassigned'", "ELSE 'Non-Omicron' END;")
  } else {
    print("No rows")
  }
  return(paste(query, collapse = " "))
}


build_display_grp_sql <- function(linDesDF) {
  query <- c(glue::glue("UPDATE {s$SEQ_INFO_TABLE} SET {s$DISPLAY_GROUP} = CASE"))

  if (nrow(linDesDF) > 0) {
    searchAll <- c(
      glue::glue("WHEN ({s$LINEAGE} IN (SELECT DISTINCT {s$LINEAGE}
                 FROM {s$SEQ_INFO_TABLE} WHERE {s$LINEAGE} like '"),
      "%')) THEN '",
      "'"
    )
    searchOne <- c(
      glue::glue("WHEN ({s$LINEAGE} = '"),
      "') THEN '",
      "'"
    )
    for (row in 1:nrow(linDesDF)) { # for every row in lineage designation file
      if (linDesDF$`Check_Descendants?`[row]) {
        searchType <- searchAll
      } else {
        searchType <- searchOne
      }

      lin <- linDesDF$Lineage[row]
      display_grp <- linDesDF$Displayed_Group[row]

      query <- c(query, paste0(searchType[1], lin, searchType[2], display_grp, searchType[3]))
    }
    query <- c(query, "ELSE 'Variants Not Currently Under Monitoring by the WHO' END;")
  } else {
    print("No rows")
  }
  return(paste(query, collapse = " "))
}


#' Pull subset DB
#'
#' @description Pull a subset database from a given "archive" database based on a given start date and end date and a variant file
#' @param archiveFile string name of archive db file or db file to subset from, e.g. "databases/All_Assay_Archive.db"
#' @param startDate string start date of data to be included in subset db, e.g. "2022-05"
#' @param endDate string end date of data to be included in subset db, e.g. "2022-07"
#' @param NewDBFile string name of new subset db file, e.g. "DB_files/DB_3M/all_May-July_SD08-15_subgroups.db"
#' @param variantFile string name of .csv file containing breakdown of lineage/subgroup names of sequences, e.g. "who_lineages_subgroups.csv"
#' @param force logical value indicating whether or not to force override already existing subset db with same name
#' @param silently logical value indicating whether to print tables created in subset db
#'
#' @export
#'
pull_subset_db <- function(archiveFile, startDate, endDate, NewDBFile,
                           variantFile,
                           force = FALSE, silently = TRUE) {
  setStrings2AssayDB()
  setSubsetStrings()

  # Write schema for new DB
  base <- system.file("sql", package = "CoMIT")
  schemaFile <- dir(base, "subset_comit_schema.sql", full.names = TRUE)

  status <- buildDBFromSchema(schemaFile, NewDBFile,
    force = force, silently = silently
  )

  subsetCon <- DBI::dbConnect(RSQLite::SQLite(), NewDBFile)
  DBI::dbExecute(subsetCon, glue::glue("attach database '{archiveFile}' as archive"))
  seqs <- paste0(subsetStrings$build_seq_info, " WHERE Collection_Date BETWEEN '", startDate, "' AND '", endDate, "';")
  clean_lineages <- clean_lineage_col()

  updates <- c(paste(subsetStrings$build_primers, seqs), clean_lineages)

  executeStatus <- runSQLstring(updates, subsetCon, silently = FALSE)


  # Build Variable statements

  WHO_lin_des <- readr::read_csv(variantFile, show_col_types = FALSE)
  addWhoDes <- build_who_des_sql(WHO_lin_des)
  addSubgroup <- build_lin_subgroup_sql(WHO_lin_des)
  addDisplayLin <- build_display_lin_sql(WHO_lin_des)
  addDisplayGrp <- build_display_grp_sql(WHO_lin_des)

  # Fill subset DB
  updates <- c(
    subsetStrings$addCol, addWhoDes,
    subsetStrings$addCol2, addSubgroup, subsetStrings$addIDmap, subsetStrings$insertOtherVarInfo
  )

  executeStatus <- runSQLstring(updates, subsetCon, silently = FALSE)

  DBI::dbDisconnect(subsetCon)

  # add display lineage and display group (dependent on subgroup column that was just added)
  subsetCon <- DBI::dbConnect(RSQLite::SQLite(), NewDBFile)
  updates <- c(subsetStrings$addCol3, addDisplayLin, subsetStrings$addCol4, subsetStrings$addCol5, addDisplayGrp, subsetStrings$addCol6, subsetStrings$addDisplayIDmap)
  executeStatus <- runSQLstring(updates, subsetCon, silently = FALSE)
  DBI::dbDisconnect(subsetCon)
}
