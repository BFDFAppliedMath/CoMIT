
loadTestPrimers <- function(primercsv, dbFile) {
  gac <- readr::read_csv(primercsv)

  con <- DBI::dbConnect(RSQLite::SQLite(), dbFile)

  uAssays <- unique(gac$`Assay Name`)
  for (assayName in uAssays) {
    rs <- DBI::dbSendStatement(
      conn = con, genInsQuery(s$ASSAY_TABLE, c(s$ASSAY_NAME)),
      params = list(assayName)
    )
    DBI::dbClearResult(rs)
  }

  assayMap <- DBI::dbGetQuery(conn = con, genSelect(s$ASSAY_TABLE))


  insCols <- c(
    s$ASSAY_NAME, s$ASSAY_ID, s$PRIMER_NAME,
    s$DIRECTION_COL, s$REACTION_COL,
    s$S_SEQ_COL, s$IS_COSPOT_STRING
  )
  primerInfoIns <- genInsQuery(s$PRIMER_TABLE, insCols)

  for (pRow in 1:nrow(gac)) {
    assay <- gac$`Assay Name`[pRow]
    assayNumDF <- assayMap[which(assayMap$Assay_Name == assay), ]
    assayNum <- assayNumDF$Assay_ID[[1]]

    rs <- DBI::dbSendStatement(
      conn = con, primerInfoIns,
      params = list(
        assay, assayNum, gac$`Primer Name`[pRow],
        gac$Direction[pRow], gac$Type[pRow],
        gac$Sequence[pRow], gac$isCospot[pRow]
      )
    )
    DBI::dbClearResult(rs)
  }

  for (pRow in 1:nrow(gac)) {
    cosPrimer <- gac$cospotPrimer[pRow]
    if (!is.na(cosPrimer)) {
      primer_ID <- DBI::dbGetQuery(
        conn = con, genSelect(s$PRIMER_TABLE,
          cols = s$PRIMER_ID,
          where = s$PRIMER_NAME
        ),
        params = list(gac$`Primer Name`[pRow])
      )
      cospot_ID <- DBI::dbGetQuery(
        conn = con, genSelect(s$PRIMER_TABLE,
          cols = s$PRIMER_ID,
          where = s$PRIMER_NAME
        ),
        params = list(cosPrimer)
      )

      # ADD to Cospot Key
      rs <- DBI::dbSendStatement(
        conn = con, genInsQuery(s$COSPOT_KEY, c(s$SPOT1, s$SPOT2)),
        params = list(primer_ID[[1]], cospot_ID[[1]])
      )
      DBI::dbClearResult(rs)
    }
  }

  DBI::dbDisconnect(con)
}


getPrimerLocationsForTest <- function(refSeqFile, DATABASE_FILE) {
  # GLOBAL VARS and Database Info
  cospotList <- list()

  # Functions---------
  checkCospotSeenList.T <- function(curr_ID) {
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

  findCospot.T <- function(currPrimerID) {
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

      print("IS a Cospot")
      return(TRUE)
    }
    return(FALSE)
  }

  # Database Extract--------
  refSeq <- seqinr::read.fasta(refSeqFile, as.string = TRUE, forceDNAtolower = FALSE)

  con <- DBI::dbConnect(RSQLite::SQLite(), DATABASE_FILE)
  cols2Sel <- c(s$SEARCHABLE_SEQ, s$PRIMER_ID, s$IS_COSPOT_STRING)
  primerDF <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_TABLE, cols2Sel))

  cospotMap <- DBI::dbGetQuery(conn = con, genSelect(
    s$COSPOT_KEY,
    c(s$SPOT1, s$SPOT2)
  ))

  currSeq <- refSeq[[1]][1]

  # Find Locations--------------
  for (primerRow in 1:nrow(primerDF)) {
    currPrimerID <- primerDF[[s$PRIMER_ID]][primerRow]

    if (checkCospotSeenList.T(currPrimerID)) {
      next
    }

    primerSeq <- primerDF[[s$SEARCHABLE_SEQ]][primerRow]

    # If primer found in sequence
    if (any(grep(primerSeq, currSeq, fixed = TRUE))) {
      # Get Start Location and calc end
      primerLoc <- stringr::str_locate_all(pattern = primerSeq, currSeq)


      # Add data to Primer Loc Table
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

    if (primerDF[[s$IS_COSPOT_STRING]][primerRow] == 1) {
      if (findCospot.T(currPrimerID)) {
        next
      }
    }
  }
  DBI::dbDisconnect(con)
}

build_Synthetic_DB <- function(DB_NAME, locationPrefix, primerCSV, refSeq) {
  # Build DB
  base <- system.file("sql", package = "CoMIT")
  schemaFile <- dir(base, "base_comit_schema.sql", full.names = TRUE)

  DatabaseName <- paste0(locationPrefix, "/", DB_NAME, ".db")

  status <- buildDBFromSchema(schemaFile, DatabaseName,
    force = TRUE, silently = TRUE
  )

  loadTestPrimers(primerCSV, DatabaseName)
  getPrimerLocationsForTest(refSeq, DatabaseName)

  return(DatabaseName)
}
