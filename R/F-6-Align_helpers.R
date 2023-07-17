# Alignment Mismatch, Ambiguity, Indel, and Perfect Match classification-------------------------------------------------
sortAlignmentPrimer <- function(currAccession, currPrimerID, currPrimerSeq, variantPrimerSeq, DATABASE_FILE, isCospot) {
  hasIndel <- FALSE
  # isSame
  if (stringi::stri_cmp_eq(currPrimerSeq, variantPrimerSeq)) {
    # Exact match, don't need to log
    return(FALSE)
  }

  currAssayID <- getAssayFromPrimer(currPrimerID)

  # Check for Ambiguities
  hasAmbiguity <- FALSE
  diffChars <- Reduce(setdiff, strsplit(c(variantPrimerSeq, currPrimerSeq), split = ""))
  for (c in diffChars) {
    if ((c != "A") & (c != "T") & (c != "G") & (c != "C") & (c != "-")) {
      hasAmbiguity <- TRUE
    }
  }
  if (hasAmbiguity) {
    addAmbiguity2DB(DATABASE_FILE, currAccession, currPrimerID, variantPrimerSeq, diffChars)
    updateVarCounts(currAssayID, isAmbiguity = TRUE)
    return(FALSE)
  } else {
    # Add to Variant Info
    # Check if already present
    if (!isVariantInDB(DATABASE_FILE, variantPrimerSeq, currPrimerID, currAccession)) {
      varID <- addFoundVariant2DB(DATABASE_FILE, variantPrimerSeq, currPrimerID, currAssayID)
      hasIndel <- classifyAlignmentVar(currPrimerSeq, variantPrimerSeq, currPrimerID, varID, currAccession, DATABASE_FILE, isCospot)
      add2IdMapTable(currAccession, varID, currPrimerID, currAssayID, DATABASE_FILE)
      updateVarCounts(currAssayID, isAmbiguity = FALSE, isWithin10bp = isVarId3Prime(varID, DATABASE_FILE))
    }
  }
  return(hasIndel)
}


classifyAlignmentVar <- function(origPrimer, variant, currPrimerID, varID, accID, DB_file, isCospot) {
  origPList <- unlist(strsplit(origPrimer, ""))
  varList <- unlist(strsplit(variant, ""))

  mismatchCount <- 0
  deletionStart <- -1
  insStart <- -1
  insertion <- c()
  indelCount <- 0

  for (i in 1:length(origPList)) {
    if (origPList[i] == "-") {
      # Add to insertion
      insertion <- c(insertion, varList[i])
      # Insertion
      if (insStart == -1) {
        insStart <- i
      }
      if ((i == length(varList)) | (origPList[i + 1] != "-")) {
        insEnd <- i
        insString <- paste(insertion, collapse = "")

        InsertionInfo2DB(insStart, insEnd, insString, nchar(origPrimer), varID, currPrimerID, DB_file, isCospot)

        indelCount <- indelCount + 1
        insStart <- -1
        insertion <- c()
      }
      next
    }
    if (origPList[i] != varList[i]) {
      varChar <- varList[i]
      # Mismatch
      if (varChar != "-") {
        # Will return with position on variant (not original primer)
        classifyMismatch(origPList[i], varChar, i, nchar(origPrimer), varID, currPrimerID, isCospot)
        mismatchCount <- mismatchCount + 1
      } else {
        # DELETION
        if (deletionStart == -1) {
          deletionStart <- i
        }
        if ((i == length(varList)) | (varList[i + 1] != "-")) {
          deletionEnd <- i

          DeletionInfo2DB(deletionStart, deletionEnd, nchar(origPrimer), varID, currPrimerID, DB_file, isCospot)
          indelCount <- indelCount + 1

          deletionStart <- -1
        }
      }
    }
  }

  # Update Mismatch and Indel count for Variations
  updateDBMismatchCount(mismatchCount, varID, cache$db)
  if (indelCount > 0) {
    updateDBIndels(varID, accID, indelCount, DB_file)
    return(TRUE)
  }
  return(FALSE)
}

InsertionInfo2DB <- function(start, end, insertionString, primerLen, varID, primerID, DBfile, isCospot = FALSE) {
  # Get Direction
  primer <- NULL
  if (isCospot) {
    primer <- dplyr::filter(cache$cospotPrimers, get(s$PRIMER_ID) == primerID)
  } else {
    primer <- dplyr::filter(cache$primerDF, get(s$PRIMER_ID) == primerID)
  }

  # Get Direction
  direction <- primer[[s$DIRECTION_COL]][1]
  mType <- insertionString
  pos <- start
  if (direction == "Forward") {
    pos <- primerLen - end + 1
  }

  length <- nchar(insertionString)

  insertToVariantTypeTableDB(DBfile, varID, "Insertion", mType, pos, primerID, indelLen = length)
}

DeletionInfo2DB <- function(start, end, primerLen, varID, primerID, DB_file, isCospot = FALSE) {
  # Get Direction
  primer <- NULL
  if (isCospot) {
    primer <- dplyr::filter(cache$cospotPrimers, get(s$PRIMER_ID) == primerID)
  } else {
    primer <- dplyr::filter(cache$primerDF, get(s$PRIMER_ID) == primerID)
  }

  # Get Direction
  direction <- primer[[s$DIRECTION_COL]][1]
  mType <- paste("Deletion from", start, "to", end)
  pos <- start
  if (direction == "Forward") {
    start <- primerLen - start + 1
    end <- primerLen - end + 1
    mType <- paste("Deletion from", end, "to", start)
    pos <- end
    length <- start - end + 1
  } else {
    length <- end - start + 1
  }

  insertToVariantTypeTableDB(DB_file, varID, "Deletion", mType, pos, primerID, indelLen = length)
}
