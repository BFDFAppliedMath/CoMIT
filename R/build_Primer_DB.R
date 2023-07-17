loadPrimersFromXL <- function(primerExcel, dbFile, numRows) { # primerXL, dbName, NumPrimerRowsInXL

  setStrings2AssayDB()

  outerRange <- paste0("B4:G", 4 + numRows)
  innerRange <- paste0("H4:K", 4 + numRows)

  gac_Outers <- readxl::read_excel(primerExcel, sheet = 1, range = outerRange, trim_ws = TRUE)
  gac_Inners <- readxl::read_excel(primerExcel, sheet = 1, range = innerRange, trim_ws = TRUE)

  con <- DBI::dbConnect(RSQLite::SQLite(), dbFile)

  lastOuterCospot <- NULL
  lastInnerCospot <- NULL

  assaySeen <- c()
  assayIns <- genInsQuery(s$ASSAY_TABLE, c(s$ASSAY_NAME, s$ASSAY_DEV_NAME))
  for (pRow in 1:nrow(gac_Outers)) {
    # Get Assay
    assay <- gac_Outers$`Assay Name on Pouch`[pRow]
    if (!(assay %in% assaySeen)) {
      rs <- DBI::dbSendStatement(
        conn = con, assayIns,
        params = list(assay, gac_Outers$`Assay Name During Development`[pRow])
      )
      DBI::dbClearResult(rs)

      assaySeen <- c(assaySeen, assay)
    }
  }

  # Set up
  assayMap <- DBI::dbGetQuery(conn = con, genSelect(s$ASSAY_TABLE))

  insCols <- c(
    s$ASSAY_NAME, s$ASSAY_DEV_NAME, s$ASSAY_ID, s$PRIMER_NAME,
    s$DIRECTION_COL, s$REACTION_COL, s$P_SEQ_COL,
    s$S_SEQ_COL, s$IS_COSPOT_STRING
  )
  primerUS <- genInsQuery(s$PRIMER_TABLE, insCols)

  primerName2IdQuery <- genSelect(s$PRIMER_TABLE, cols = s$PRIMER_ID, where = s$PRIMER_NAME)

  cospotIns <- genInsQuery(s$COSPOT_KEY, c(s$SPOT1, s$SPOT2))

  for (pRow in 1:nrow(gac_Outers)) {
    isCospot <- FALSE

    assay <- gac_Outers$`Assay Name on Pouch`[pRow]
    assayNumDF <- assayMap[which(assayMap$Assay_Name == assay), ]
    assayNum <- assayNumDF$Assay_ID[[1]]

    # Outer
    if (is.na(gac_Outers$Name[pRow])) {
      isCospot <- TRUE
    } else {
      searchable <- gac_Outers$Sequence[pRow]
      if (gac_Outers$Direction[pRow] == "Reverse") {
        searchable <- seqinr::c2s(rev(seqinr::comp(seqinr::s2c(searchable),
          ambiguous = TRUE,
          forceToLower = FALSE
        )))
      }

      primerAttributes <- list(
        assay, gac_Outers$`Assay Name During Development`[pRow],
        assayNum, gac_Outers$Name[pRow],
        gac_Outers$Direction[pRow], "Outer",
        gac_Outers$Sequence[pRow], searchable, gac_Outers$Cospot[pRow]
      )

      rs <- DBI::dbSendStatement(conn = con, primerUS, params = primerAttributes)
      DBI::dbClearResult(rs)

      if (gac_Outers$Cospot[pRow] == 1) {
        if (!is.null(lastOuterCospot)) {
          # Get curr ID
          cospot_ID <- DBI::dbGetQuery(
            conn = con, primerName2IdQuery,
            params = list(gac_Outers$Name[pRow])
          )
          currOuterCospot <- cospot_ID[[1]][1]

          # ADD to Cospot Key
          rs <- DBI::dbSendStatement(
            conn = con, genInsQuery(s$COSPOT_KEY, c(s$SPOT1, s$SPOT2)),
            params = list(currOuterCospot, lastOuterCospot)
          )
          DBI::dbClearResult(rs)

          rs <- DBI::dbSendStatement(
            conn = con, genInsQuery(s$COSPOT_KEY, c(s$SPOT1, s$SPOT2)),
            params = list(lastOuterCospot, currOuterCospot)
          )
          DBI::dbClearResult(rs)

          # Reset
          lastOuterCospot <- NULL
        } else {
          # Add for next hit
          cospot_ID <- DBI::dbGetQuery(
            conn = con, primerName2IdQuery,
            params = list(gac_Outers$Name[pRow])
          )
          lastOuterCospot <- cospot_ID[[1]][1]
        }
      }
    }


    # Inner
    searchable <- gac_Inners$Sequence[pRow]
    if (gac_Inners$Direction[pRow] == "Reverse") {
      searchable <- seqinr::c2s(rev(seqinr::comp(seqinr::s2c(searchable),
        ambiguous = TRUE,
        forceToLower = FALSE
      )))
    }
    primerAttributes <- list(
      assay, gac_Outers$`Assay Name During Development`[pRow],
      assayNum, gac_Inners$Name[pRow],
      gac_Inners$Direction[pRow], "Inner",
      gac_Inners$Sequence[pRow], searchable,
      gac_Inners$Cospot[pRow]
    )

    rs <- DBI::dbSendStatement(
      conn = con, primerUS,
      params = primerAttributes
    )
    DBI::dbClearResult(rs)

    if (gac_Inners$Cospot[pRow] == 1) {
      if (!is.null(lastOuterCospot)) {
        # Get curr ID
        cospot_ID <- DBI::dbGetQuery(
          conn = con, primerName2IdQuery,
          params = list(gac_Inners$Name[pRow])
        )
        currOuterCospot <- cospot_ID[[1]][1]

        # ADD to Cospot Key
        rs <- DBI::dbSendStatement(
          conn = con, cospotIns,
          params = list(currOuterCospot, lastOuterCospot)
        )
        DBI::dbClearResult(rs)

        rs <- DBI::dbSendStatement(
          conn = con, cospotIns,
          params = list(lastOuterCospot, currOuterCospot)
        )
        DBI::dbClearResult(rs)

        # Reset
        lastOuterCospot <- NULL
      } else {
        # Add for next hit
        cospot_ID <- DBI::dbGetQuery(
          conn = con, primerName2IdQuery,
          params = list(gac_Inners$Name[pRow])
        )
        lastOuterCospot <- cospot_ID[[1]][1]
      }
    }
  }

  DBI::dbDisconnect(con)
}


getPrimerLocations <- function(refSeqFile, DATABASE_FILE) {

  # GLOBAL VARS and Database Info
  cospotList <- list()

  # Functions---------
  checkCospotSeenList <- function(curr_ID) {
    # Check cospotList
    if (length(cospotList) == 0) {
      return(FALSE)
    }
    for (id in 1:length(cospotList)) {
      if (curr_ID == cospotList[[id]]) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  findCospot <- function(currPrimerID) {
    # Cospot
    cospotPair <- cospotMap[which(cospotMap$Spot == currPrimerID), ]
    coSpot_ID <- cospotPair$Cospot[1]

    cospot <- DBI::dbGetQuery(
      conn = con, genSelect(s$PRIMER_TABLE,
        cols = s$SEARCHABLE_SEQ,
        where = s$PRIMER_ID
      ),
      params = list(coSpot_ID)
    )


    primerSeq <- cospot$Searchable_Seq[1]

    # If primer found in sequence
    if (any(grep(primerSeq, currSeq, fixed = TRUE))) {
      cospotList <- c(cospotList, coSpot_ID)

      primerLoc <- stringr::str_locate_all(pattern = primerSeq, currSeq)

      # Add data to Seq_Info Table
      rs <- DBI::dbSendStatement(
        conn = con, genInsQuery(
          s$PRIMER_LOCS_TABLE,
          c(
            s$PRIMER_ID,
            s$LOC_START,
            s$LOC_END
          )
        ),
        params = list(
          currPrimerID,
          primerLoc[[1]][1],
          primerLoc[[1]][2]
        )
      )
      DBI::dbClearResult(rs)

      # print("IS a Cospot")
      return(TRUE)
    }
    return(FALSE)
  }

  # Database Extract--------
  refseq <- seqinr::read.fasta(refSeqFile, as.string = TRUE, forceDNAtolower = FALSE)

  con <- DBI::dbConnect(RSQLite::SQLite(), DATABASE_FILE)
  cols2Sel <- c(s$SEARCHABLE_SEQ, s$PRIMER_ID, s$IS_COSPOT_STRING)
  primerDF <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_TABLE, cols2Sel))

  cospotMap <- DBI::dbGetQuery(conn = con, genSelect(
    s$COSPOT_KEY,
    c(s$SPOT1, s$SPOT2)
  ))

  currSeq <- refseq[[1]][1]

  # Find Locations--------------
  for (primerRow in 1:nrow(primerDF)) {
    currPrimerID <- primerDF$Primer_ID[primerRow]

    if (checkCospotSeenList(currPrimerID)) {
      next
    }

    primerSeq <- primerDF$Searchable_Seq[primerRow]

    # If primer found in sequence
    if (any(grep(primerSeq, currSeq, fixed = TRUE))) {
      # Get Start Location and calc end
      primerLoc <- stringr::str_locate_all(pattern = primerSeq, currSeq)


      # Add data to Seq_Info Table
      rs <- DBI::dbSendStatement(
        conn = con, genInsQuery(
          s$PRIMER_LOCS_TABLE,
          c(
            s$PRIMER_ID,
            s$LOC_START,
            s$LOC_END
          )
        ),
        params = list(
          currPrimerID,
          primerLoc[[1]][1],
          primerLoc[[1]][2]
        )
      )
      DBI::dbClearResult(rs)
    }

    if (primerDF$is_Cospot[primerRow] == 1) {
      if (findCospot(currPrimerID)) {
        next
      }
    }
  }
  DBI::dbDisconnect(con)
}



# adds inner_amp_len, outer_amp_len, and inner_gc to assay info table based on ref seq
addAmpliconInfoToDB <- function(refSeqFile, DATABASE_FILE) {

  if (is.null(refSeqFile)) {
    base <- system.file("cov_data", package = "CoMIT")
    refSeqFile <- dir(base, "refseq.fasta", full.names = TRUE)
    refSeq <- seqinr::read.fasta(refSeqFile, as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)[[1]][1]
  } else {
    refSeq <- seqinr::read.fasta(refSeqFile, as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)[[1]][1]
  }


  # Get primer information from db
  con <- DBI::dbConnect(RSQLite::SQLite(), DATABASE_FILE)

  primerDF <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_TABLE)) %>% dplyr::collect()
  primerLocs <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_LOCS_TABLE)) %>% dplyr::collect()

  primerInfo <- primerDF %>% dplyr::left_join(primerLocs, by = "Primer_ID")

  DBI::dbDisconnect(con)


  # Loop through assays
  for (assay in unique(primerInfo[[s$ASSAY_ID]])) {
    # filter for assay
    assayDF <- primerInfo %>% dplyr::filter(get(s$ASSAY_ID) == assay)

    IF_start <- assayDF %>%
      dplyr::filter(get(s$REACTION_COL) == "Inner", get(s$DIRECTION_COL) == "Forward") %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(get(s$LOC_START))

    IR_end <- assayDF %>%
      dplyr::filter(get(s$REACTION_COL) == "Inner", get(s$DIRECTION_COL) == "Reverse") %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(get(s$LOC_END))

    inner_amp_len <- IR_end - IF_start + 1


    OF_start <- assayDF %>%
      dplyr::filter(get(s$REACTION_COL) == "Outer", get(s$DIRECTION_COL) == "Forward") %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(get(s$LOC_START))

    OR_end <- assayDF %>%
      dplyr::filter(get(s$REACTION_COL) == "Outer", get(s$DIRECTION_COL) == "Reverse") %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(get(s$LOC_END))

    outer_amp_len <- OR_end - OF_start + 1


    inner_amp_seq <- substr(refSeq, IF_start, IR_end)

    num_g <- stringr::str_count(inner_amp_seq, "G")
    num_c <- stringr::str_count(inner_amp_seq, "C")
    GC <- (num_g + num_c) / inner_amp_len


    # Add data to Seq_Info Table
    con <- DBI::dbConnect(RSQLite::SQLite(), DATABASE_FILE)
    insCols <- c(s$INNER_AMP_LEN, s$OUTER_AMP_LEN, s$INNER_GC)

    # this only works when database is set up with build covid db, not with other methods like build db from schema that are called by some tests like test-classify-syn-amplicons
    rs1 <- DBI::dbSendStatement(
      conn = con, genUpdate2Val(
        s$ASSAY_TABLE,
        s$INNER_AMP_LEN,
        s$ASSAY_ID
      ),
      params = list(inner_amp_len, assay)
    )
    DBI::dbGetRowsAffected(rs1)
    DBI::dbClearResult(rs1)


    rs2 <- DBI::dbSendStatement(
      conn = con, genUpdate2Val(
        s$ASSAY_TABLE,
        s$OUTER_AMP_LEN,
        s$ASSAY_ID
      ),
      params = list(outer_amp_len, assay)
    )
    DBI::dbGetRowsAffected(rs2)
    DBI::dbClearResult(rs2)


    rs3 <- DBI::dbSendStatement(
      conn = con, genUpdate2Val(
        s$ASSAY_TABLE,
        s$INNER_GC,
        s$ASSAY_ID
      ),
      params = list(GC, assay)
    )

    DBI::dbGetRowsAffected(rs3)
    DBI::dbClearResult(rs3)


    DBI::dbDisconnect(con)
  }
}



buildBaseDbFromXL <- function(dbName, primerXL, NumPrimerRowsInXL, refSeqFile,
                              force = FALSE, silently = TRUE) {
  # Build DB
  base <- system.file("sql", package = "CoMIT")
  schemaFile <- dir(base, "base_comit_schema.sql", full.names = TRUE)

  status <- buildDBFromSchema(schemaFile, dbName,
    force = force, silently = silently
  )
  if (!status) {
    warning("Database could not be built")
    return()
  }

  # Fill DB
  loadPrimersFromXL(primerXL, dbName, NumPrimerRowsInXL)
  getPrimerLocations(refSeqFile, dbName)
  addAmpliconInfoToDB(refSeqFile, dbName)
}



#' Build COVID-19 Analysis Database
#'
#' @description Create a database given desired name and location, with primer information from the BioFire Defense Covid-19 Test assays
#' @param DB_NAME Desired database name, without filepath or .db extension
#' @param locationPrefix Folder location for database, defaults to home
#' @param include_all_assays Include primer information from all 7 assays or only A,D, and E assays?
#' @param force Force overwrite if file exists?
#' @param silently Turn off printed updates?
#'
#' @export
#'
#' @return full database file location name
build_COVID_DB <- function(DB_NAME, locationPrefix = "~", include_all_assays = TRUE, force = FALSE, silently = TRUE) {

  base <- system.file("cov_data", package = "CoMIT")

  refSeqFile <- dir(base, "refseq.fasta", full.names = TRUE)


  # Permanent Reference Files for COVID---------------------------
  if (include_all_assays) {
    # All COVID Assays
    primerXL <- dir(base, "COVID_GAC_Form_2021.xlsx", full.names = TRUE)
    NumPrimerRowsInXL <- 16
  } else {
    # Only unmasked assays
    primerXL <- dir(base, "COVID_unmasked_GAC_Form.xlsx", full.names = TRUE)
    NumPrimerRowsInXL <- 7
  }

  DatabaseName <- paste0(locationPrefix, DB_NAME, ".db")

  buildBaseDbFromXL(
    DatabaseName, primerXL, NumPrimerRowsInXL,
    refSeqFile, force, silently
  )

  return(DatabaseName)
}
