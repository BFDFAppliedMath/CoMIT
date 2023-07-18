#' Insert location and sequencing information
#'
#' @param DB filepath of CoMIT database
#' @param metadata filepath of metadata
#' @noRd
#'
add_extra_metadata <- function(DB, metadata) {
  `1` <- `2` <- `3` <- `Accession ID` <- `Sequencing technology` <- `Assembly method` <-
    `Country` <- `Location` <- `Continent` <- `Location_Details` <- `Pull_Date` <- NULL

  meta <- readr::read_tsv(metadata)

  # Filter for desired columns------------------------------------------
  if ("Sequencing technology" %in% names(meta)) {
    meta <- meta %>% dplyr::select("Accession ID", "Location", "Sequencing technology", "Assembly method")
  } else {
    meta <- meta %>%
      dplyr::select("Accession ID", "Location") %>%
      dplyr::mutate("Sequencing technology" = NA, "Assembly method" = NA)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)

  # Retrieve accessions to subset by------------------------------------
  acc_seq_info <- dplyr::tbl(con, "Seq_Info") %>%
    dplyr::filter(Pull_Date >= "2022-01-10") %>%
    dplyr::select(s$ACC_ID) %>%
    dplyr::collect()
  acc_high_n <- dplyr::tbl(con, "High_N_Seqs") %>%
    dplyr::select(s$ACC_ID) %>%
    dplyr::collect()

  # DBI::dbListTables(con)

  # Data cleaning--------------------------------------------------------
  split_location <- stringr::str_split_fixed(meta$Location, "/", 3)
  meta <- cbind(meta, split_location)
  meta <- meta %>% dplyr::rename("Continent" = `1`, "Country" = `2`, "Location_Details" = `3`, "Accession_ID" = `Accession ID`)

  high_n <- meta %>%
    dplyr::select(s$ACC_ID, "Sequencing_Method" = `Sequencing technology`, "Assembly_Method" = `Assembly method`, Country) %>%
    dplyr::inner_join(acc_high_n, by = "Accession_ID")

  seq_meta <- meta %>%
    dplyr::select(-Location) %>%
    dplyr::inner_join(acc_seq_info, by = "Accession_ID") %>%
    dplyr::select(s$ACC_ID, Continent, Country, Location_Details)

  # Add metadata to database---------------------------------------------
  for (row in 1:nrow(seq_meta)) {
    accessionNum <- seq_meta$Accession_ID[row]
    if (isAccInMetaDB(accessionNum, DB)) {
      # Next sequence
      next
    }
    insert2SeqMeta(DB, accessionNum, seq_meta$Continent[row], seq_meta$Country[row], seq_meta$Location_Details[row])
  }

  for (row in 1:nrow(high_n)) {
    accessionNum <- high_n$Accession_ID[row]

    rs <- DBI::dbSendStatement(
      conn = con, "UPDATE High_N_Seqs SET Sequencing_Method=?, Assembly_Method=?, Country=? WHERE Accession_ID=?;",
      params = list(
        high_n$Sequencing_Method[row],
        high_n$Assembly_Method[row],
        high_n$Country[row],
        accessionNum
      )
    )
    DBI::dbClearResult(rs)
  }

  DBI::dbDisconnect(con)
}


#' Update DB (Add sequence info to a database)
#'
#' @description CoMIT Database update. Requires an existing CoMIT database. Pulls from the CSV specified
# for fasta_list_file to find locations of fasta and metadata files to use.
#'
#' @param dbName string name of the database, with no filepath or .db extension (e.g. "All_Assay_Archive")
#' @param dbFile string file path of the database, with .db extension
#' @param saveCopy logical value; if true, will save a copy of the old DB called dbName to _DB_files/DB_Copies in the user's home directory
#' @param pullDate_1 string date for beginning of date selection used in the update (refers to pullDate column in fasta_list_file, only things >= this date will be included in comit_classify run), e.g. "2022-10-11"
#' @param pullDate_2 string date for end of date selection used in the update (refers to pullDate column in fasta_list_file, only things <= this date will be included in comit_classify run); if NA, use everything after pullDate_1
#' @param fasta_list_file string name of the csv file with fasta file information, including path to file; has columns folder (location of each fasta file), fastaFile (fasta file name), pullDate, metaDataFile (NA if same name as fasta file, name if different)
#' @param run_histories_filepath string path to .txt file used to record run histories for this database (file will be called run_histories_file + db name + .txt)
#'
#' @export
#'
update_db <- function(dbName, dbFile, saveCopy, pullDate_1, pullDate_2, fasta_list_file, run_histories_filepath) {
  # copy old db with the same name to file with date appended
  if (saveCopy) {
    if (file.exists(dbFile) && !file.exists(paste0("/home/", Sys.info()[["user"]], "/_DB_files/DB_Copies/", dbName, "_", Sys.Date(), ".db"))) {
      system(paste("cp", dbFile, paste0("/home/", Sys.info()[["user"]], "/_DB_files/DB_Copies")))
      Sys.setenv(TZ = "US/Mountain")
      file.rename(
        from = paste0("/home/", Sys.info()[["user"]], "/_DB_files/DB_Copies/", dbName, ".db"),
        to = paste0("/home/", Sys.info()[["user"]], "/_DB_files/DB_Copies/", dbName, "_", Sys.Date(), ".db")
      )
    }
  }

  # Get File Names
  fastaDF <- readr::read_csv(fasta_list_file)

  if (is.na(pullDate_2)) {
    fastaDF <- fastaDF %>% dplyr::filter(pullDate >= pullDate_1)
  } else {
    fastaDF <- fastaDF %>% dplyr::filter(pullDate >= pullDate_1 & pullDate <= pullDate_2)
  }

  # START RUN-------------------------------------------
  run_file_name <- paste0(run_histories_filepath, dbName, ".txt")
  write(paste("Update Run for", dbFile), file = run_file_name, append = TRUE)
  cache$TotalSeqsRun <- 0

  overall_start <- Sys.time()

  totalFiles <- length(fastaDF$fastaFile)

  for (i in 1:totalFiles) {
    FASTA <- paste0(fastaDF$folder[i], "/", fastaDF$fastaFile[i], ".fasta")
    pullDate <- fastaDF$pullDate[i]
    metaD <- paste0(fastaDF$folder[i], "/meta/", fastaDF$fastaFile[i], ".tsv")
    if (!is.na(fastaDF$metaDataFile[i])) {
      metaD <- fastaDF$metaDataFile[i]
    }

    myvars <- readr::cols_only(`Accession ID` = "c", Clade = "c", Lineage = "c")
    cladeDF <- readr::read_tsv(metaD, col_types = myvars)

    print(paste(i, "out of", totalFiles, "files"))

    outputvector <- comit_classify(FASTA, as.character(pullDate), dbFile, cladeDF)

    write(outputvector, file = run_file_name, append = TRUE)
    print("metadata file:")
    print(metaD)

    add_extra_metadata(dbFile, metaD)
  }

  overall_finish <- Sys.time()

  time_statment <- paste("Start:", overall_start, "End:", overall_finish)
  totalRunTime <- difftime(overall_finish, overall_start, units = "hours")[[1]]
  total_time <- paste("Total Run Time:", totalRunTime, "hours")

  classifyRate <- cache$TotalSeqsRun / difftime(overall_finish, overall_start, units = "mins")[[1]]
  rateStatement <- paste("Sequence Classification Rate:", classifyRate, "seqs/min")

  write(c(time_statment, total_time, rateStatement, " ", " ", " ", "_______________________________________________________________________"),
    file = run_file_name, append = TRUE
  )

}
