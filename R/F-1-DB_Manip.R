#' Insert variant details into the Variant Type table
#'
#' @param DBfile The database to add the variant details to
#' @param varID The Variant ID
#' @param class The type of change (Mismatch, Insertion, Deletion)
#' @param changeString String describing the base change in the variant,
#'        (Reference->Variant i.e. "A->G", "Deletion from {start postion} to {end position]}", Insertion Chars i.e. "ACTT"
#' @param pos The position from the three prime end of the mismatch or the beginning position of an insertion or deletion
#' @param primerID The ID of the primer that the variant is associated with
#' @param typeString String describing the Mismatch in the variant, Primer/Template i.e. "A/C"
#'        Only used for mismatches, defaults to empty string
#'
#' @return NULL
#'
#' @examples
#' # Insert a mismatch
#' # insertToVariantTypeTableDB("Test.db", 2, "Mismatch", "T->G", 5, 3, "T/C")
#'
#' @noRd
insertToVariantTypeTableDB <- function(DBfile, varID, class, changeString, pos,
                                       primerID, typeString = "", indelLen = NA) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DBfile)

  varTypeCols <- c(
    s$P_VAR_ID, s$VAR_TYPE, s$SEQ_CHANGE,
    s$MISMATCH_TYPE, s$POS_FROM_SIDE,
    s$INDEL_LEN_COL, s$PRIMER_ID
  )

  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(
      s$VAR_TYPE_TABLE,
      varTypeCols
    ),
    params = list(
      varID, class, changeString, typeString,
      pos, indelLen, primerID
    )
  )

  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}


#' Update Mismatch Count in a Database
#'
#'
#' @param mismatchCount New Mismatch Count
#' @param varID Var Id to edit
#' @param DB Database filepath
#'
#' @return NULL
#'
#' @noRd
updateDBMismatchCount <- function(mismatchCount, varID, DB) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Add Mismatch Count
  rs <- DBI::dbSendStatement(
    conn = con,
    genUpdate2Val(
      s$VAR_INFO_TABLE,
      s$MISMATCH_NUM,
      s$VAR_ID
    ),
    params = list(mismatchCount, varID)
  )

  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}


#' Add a row to the ID_Map Table in a Database
#'
#' @param accID Accession ID to add
#' @param varID VarID to add
#' @param primerID PrimerID associated with varID
#' @param assayID AssayID associated with Primer
#' @param DB_file Database filepath
#'
#' @return NULL
#'
#' @noRd
add2IdMapTable <- function(accID, varID, primerID, assayID, DB_file) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_file)
  varListCols <- c(s$ACC_ID, s$VAR_ID, s$PRIMER_ID, s$ASSAY_ID)
  # Add entry to Variation List table in DB
  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(
      s$ID_MAP_TABLE,
      varListCols
    ),
    params = list(accID, varID, primerID, assayID)
  )

  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}

insertToHighSeqsDB <- function(DB_file, accessionNum, percentN) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_file)
  highNCols <- c(s$ACC_ID, s$SEQ_PER_N)
  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(
      s$HIGH_N_TABLE,
      highNCols
    ),
    params = list(accessionNum, percentN)
  )
  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}

insert2SeqMeta <- function(DB_file, accessionNum, continent, country, locDetails) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_file)
  metaCols <- c(s$ACC_ID, s$CONTINENT, s$COUNTRY, s$LOC_DETAILS)
  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(
      s$SEQ_META_TABLE,
      metaCols
    ),
    params = list(accessionNum, continent, country, locDetails)
  )
  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}

# Removed subDate param (causing errors with db insertion)
insert2SeqInfo <- function(DB_file, accessionNum, location, colDate, colDay,
                           FASTA_header, pullDate, percentN, clade, lineage) { # add check for if percentN is not a number
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_file)
  # Add data to Seq_Info Table
  seqInfoCols <- c(
    s$ACC_ID, s$LOC_ACC,
    s$COLLECTION_DATE, s$COLLECTION_DAY,
    s$FASTA_HEADER, s$PULL_DATE,
    s$PER_N, s$LINEAGE,
    s$CLADE
  )

  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(
      s$SEQ_INFO_TABLE,
      seqInfoCols
    ),
    params = list(
      accessionNum, location, colDate, colDay,
      FASTA_header, pullDate, percentN, lineage,
      clade
    )
  )
  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}

insert2SeqInfoNoMeta <- function(DB_file, accessionNum, location, colDate,
                                 colDay, FASTA_header, pullDate, percentN) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_file)
  # Add data to Seq_Info Table
  seqInfoCols <- c(
    s$ACC_ID, s$LOC_ACC,
    s$COLLECTION_DATE, s$COLLECTION_DAY,
    s$FASTA_HEADER, s$PULL_DATE, s$PER_N
  )

  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(
      s$SEQ_INFO_TABLE,
      seqInfoCols
    ),
    params = list(
      accessionNum, location, colDate, colDay,
      FASTA_header, pullDate, percentN
    )
  )
  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}

addFoundVariant2DB <- function(DBfile, variant, PrimerID, AssayID) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DBfile)
  # If not present Add New Variant
  varInfoCols <- c(s$VAR_SEQ_COL, s$PRIMER_ID, s$ASSAY_ID)

  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(
      s$VAR_INFO_TABLE,
      varInfoCols
    ),
    params = list(variant, PrimerID, AssayID)
  )

  DBI::dbClearResult(rs)

  query <- paste0("select last_insert_rowid();")
  varID <- DBI::dbGetQuery(conn = con, query)
  DBI::dbDisconnect(con)
  return(varID[[1]])
}

# when a bad amplicon is found in comit_classify; accID, assay, OAL, IAL, GC
addBadAmplicon2DB <- function(accID, assayID, OAL, IAL, GC, outOfOrder = FALSE, duplicated = FALSE) { # accID, assay, OAL, IAL, GC
  con <- DBI::dbConnect(RSQLite::SQLite(), cache$db)

  insCols <- c(s$ACC_ID, s$ASSAY_ID, s$OUTER_AMP_LEN, s$INNER_AMP_LEN, s$INNER_GC, s$OUT_OF_ORDER, s$DUPLICATED)

  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(s$AMPLICON_TABLE, insCols),
    params = list(accID, assayID, OAL, IAL, GC, outOfOrder, duplicated)
  )

  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}


# Align Script--------------------------------------------------------------------

# Check if exists, if exists will update counts (freq and indel on seq_info) and add to ID_Map table

isVariantInDB <- function(DBfile, variant, currPrimerID, accID) {
  # Check if already present in DB
  con <- DBI::dbConnect(RSQLite::SQLite(), DBfile)
  varInfoCols <- c(s$VAR_ID, s$INDELS_PRESENT_VAR)
  whereCols <- c(s$VAR_SEQ_COL, s$PRIMER_ID)

  variantExists <- DBI::dbGetQuery(con, genSelect(
    s$VAR_INFO_TABLE,
    varInfoCols,
    whereCols
  ),
  params = list(variant, currPrimerID)
  )

  # Variant Already present
  if (nrow(variantExists) != 0) {
    varID <- variantExists$Var_ID[1]
    # Add entry to Variation List table in DB
    varListCols <- c(s$ACC_ID, s$VAR_ID, s$PRIMER_ID, s$ASSAY_ID)

    rs <- DBI::dbSendStatement(
      conn = con, genInsQuery(
        s$ID_MAP_TABLE,
        varListCols
      ),
      params = list(accID, varID, currPrimerID, getAssayFromPrimer(currPrimerID))
    )

    DBI::dbClearResult(rs)

    if (variantExists$Has_Indel[1] == 1) {
      # Update indel info for sequence
      rs <- DBI::dbSendStatement(
        conn = con,
        genUpdateAugVal(
          s$SEQ_INFO_TABLE,
          s$INDELS_PRESENT_VAR,
          s$ACC_ID
        ),
        params = list(accID)
      )

      DBI::dbClearResult(rs)
    }

    DBI::dbDisconnect(con)

    # Update global counts for Seq Info table
    updateVarCounts(getAssayFromPrimer(currPrimerID), isAmbiguity = FALSE, isWithin10bp = isVarId3Prime(varID, DBfile))

    return(TRUE)
  }

  DBI::dbDisconnect(con)
  # NOT FOUND IN DATABASE
  return(FALSE)
}


addAmbiguity2DB <- function(DBfile, currAccession, currPrimerID, primerString, ambiguityChars) {
  N_Percent <- round((stringr::str_count(primerString, "N") / nchar(primerString) * 100))

  con <- DBI::dbConnect(RSQLite::SQLite(), DBfile)
  # Add to Ambiguity Table
  ambSeqsCols <- c(s$ACC_ID, s$AMB_CHARS_COL, s$VAR_SEQ_COL, s$PRIMER_ID, s$SEQ_PER_N)

  rs <- DBI::dbSendStatement(
    conn = con, genInsQuery(
      s$AMBIGUOUS_SEQS_TABLE,
      ambSeqsCols
    ),
    params = list(currAccession, toString(ambiguityChars), primerString, currPrimerID, N_Percent)
  )

  DBI::dbClearResult(rs)

  # Fix Ambiguity Count
  rs <- DBI::dbSendStatement(
    conn = con,
    genUpdateAugVal(
      s$SEQ_INFO_TABLE,
      s$HAS_AMB_COL,
      s$ACC_ID
    ),
    params = list(currAccession)
  )

  DBI::dbClearResult(rs)

  DBI::dbDisconnect(con)
}

updateDBvarAndAssayCount <- function(currAccession, varCounts, assayCounts, assayAmb, assay10bp, DBfile) { # test if db is populated
  con <- DBI::dbConnect(RSQLite::SQLite(), DBfile)

  # Assay Counts for mutations
  rs <- DBI::dbSendStatement(
    conn = con,
    genUpdate2Val(
      s$SEQ_INFO_TABLE,
      s$ASSAYS_AFFECTED,
      s$ACC_ID
    ),
    params = list(assayCounts, currAccession)
  )

  DBI::dbClearResult(rs)

  # VarCount
  rs <- DBI::dbSendStatement(
    conn = con,
    genUpdateAugValCustom(
      s$SEQ_INFO_TABLE,
      s$SEQ_VAR_COUNT,
      s$ACC_ID
    ),
    params = list(varCounts, currAccession)
  )

  DBI::dbClearResult(rs)

  # Assays Count with mutations or ambiguities
  rs <- DBI::dbSendStatement(
    conn = con,
    genUpdate2Val(
      s$SEQ_INFO_TABLE,
      s$BP10_AA,
      s$ACC_ID
    ),
    params = list(assay10bp, currAccession)
  )

  DBI::dbClearResult(rs)

  # Assays Count with mutations in a primer on the 3 primer end
  rs <- DBI::dbSendStatement(
    conn = con,
    genUpdate2Val(
      s$SEQ_INFO_TABLE,
      s$WITH_AMBS_AA,
      s$ACC_ID
    ),
    params = list(assayAmb, currAccession)
  )

  DBI::dbClearResult(rs)

  DBI::dbDisconnect(con)
}

updateDBIndels <- function(varID, accID, indelNum, DB) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Update indel info for var into
  rs <- DBI::dbSendStatement(
    conn = con,
    genUpdate2Val(
      s$VAR_INFO_TABLE,
      s$INDELS_PRESENT_VAR,
      s$VAR_ID
    ),
    params = list(indelNum, varID)
  )

  DBI::dbClearResult(rs)

  rs <- DBI::dbSendStatement(
    conn = con,
    genUpdateAugVal(
      s$SEQ_INFO_TABLE,
      s$INDELS_PRESENT_VAR,
      s$ACC_ID
    ),
    params = list(accID)
  )

  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}



# Retrieval Functions-----------------------------------------
isVarId3Prime <- function(varId, DB_File) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_File)
  # Add entry to Variation List table in DB
  rs <- DBI::dbGetQuery(
    conn = con,
    genSelect(
      s$VAR_TYPE_TABLE,
      s$POS_FROM_SIDE,
      s$VAR_ID
    ),
    params = list(varId)
  )

  pos <- min(rs$Pos_from_3P)
  DBI::dbDisconnect(con)
  if (pos <= 10) {
    return(TRUE)
  } else {
    (
      return(FALSE)
    )
  }
}

getAssayFromPrimer <- function(primerID) {
  assayCol <- cache$assayKey %>% dplyr::filter(get(s$PRIMER_ID) == primerID)
  return(assayCol$Assay_ID[1])
}

isAccInDB <- function(accID, DB_File) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_File)
  # Check if already in database
  rs <- DBI::dbGetQuery(con, paste0(
    "SELECT COUNT(1) FROM '", s$SEQ_INFO_TABLE,
    "' WHERE \"", s$ACC_ID, "\" = ?;"
  ),
  params = list(accID)
  )
  DBI::dbDisconnect(con)

  return(rs$`COUNT(1)` > 0)
}

isAccInMetaDB <- function(accID, DB_File) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_File)
  # Check if already in database
  rs <- DBI::dbGetQuery(con, paste0(
    "SELECT COUNT(1) FROM '", s$SEQ_META_TABLE,
    "' WHERE \"", s$ACC_ID, "\" = ?;"
  ),
  params = list(accID)
  )
  DBI::dbDisconnect(con)

  return(rs$`COUNT(1)` > 0)
}

isAccInHighNDB <- function(accID, DB_File) {
  # Check if already in database
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_File)
  rs <- DBI::dbGetQuery(con, paste0(
    "SELECT COUNT(1) FROM '", s$HIGH_N_TABLE,
    "' WHERE \"", s$ACC_ID, "\" = ?;"
  ),
  params = list(accID)
  )
  DBI::dbDisconnect(con)

  return(rs$`COUNT(1)` > 0)
}
