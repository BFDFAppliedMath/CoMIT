#' Check a sequence for variations on cospoted primers and log in database
#'
#' @param possVarList Vector of possible mismatches
#' @param cSeq1 Searchable Sequence of cospot 1
#' @param cSeq2 Searchable Sequence of cospot 2
#' @param co_ID_1 Identifier for cospot 1
#' @param co_ID_2 Identifier for cospot 2
#' @param seq2search String Sequence to search for primers
#'
#' @return VarID of classified variant or -1 if no matching variant found
#' #checkPossCospotVarList({vector: mismatchList}, "ACTTGTGCCC", "ACTTGTGCGC", 10, 11, {character vector: dnaSequence})
#' @noRd

checkPossCospotVarList <- function(possVarList, cSeq1, cSeq2, co_ID_1, co_ID_2, assayID, seq2search, wholeLeftWindow, wholeRightWindow) {
  varID <- -1
  for (i in 1:length(possVarList)) {
    if (grepl(possVarList[i], seq2search, fixed = TRUE)) {
      variant <- possVarList[[i]][1]
      co1List <- unlist(strsplit(cSeq1, ""))
      co2List <- unlist(strsplit(cSeq2, ""))
      varList <- unlist(strsplit(variant, ""))

      diffCount1 <- 0
      diffCount2 <- 0

      # Which primer is closer
      for (i in 1:nchar(variant)) {
        if (co1List[i] != varList[i]) {
          diffCount1 <- diffCount1 + 1
        }

        if (co2List[i] != varList[i]) {
          diffCount2 <- diffCount2 + 1
        }
      }
      if (diffCount1 <= diffCount2) {
        # findPrimerAddLoc(co_ID_1, assayID, cSeq1, seq2search, wholeLeftWindow, wholeRightWindow)

        # If found Classify
        varID <- classifyVariant(cSeq1, variant, co_ID_1, isCospot = TRUE)

        # Terminate iteration through list if found
        break
      } else {
        # findPrimerAddLoc(co_ID_2, assayID, cSeq2, seq2search, wholeLeftWindow, wholeRightWindow)

        varID <- classifyVariant(cSeq2, variant, co_ID_2, isCospot = TRUE)

        # Terminate iteration through list if found
        break
      }
    }
  }

  return(varID)
}

#' Search for variation not already present in the database
#'
#' @param cSeq1 Searchable Sequence of cospot 1
#' @param cSeq2 Searchable Sequence of cospot 2
#' @param co_ID_1 Database Identifier for cospot 1
#' @param co_ID_2 Database Identifier for cospot 2
#' @param seq2search String Sequence to search for primers
#'
#' @return VarID number recieved from the database or -1 if no matching variant found
#'
#' @noRd
#'
findNovelCospotVar <- function(cSeq1, cSeq2, co_ID_1, co_ID_2, assayID, seq2search, wholeLeftWindow, wholeRightWindow) {
  # Get VARIANTS OF BOTH COMBINE LIST AND DEDUPLICATE
  possibleVars1 <- readRDS(paste0(cache$mismatchLibrary, "/", cSeq1, "_1.rds"))
  possibleVars2 <- readRDS(paste0(cache$mismatchLibrary, "/", cSeq2, "_1.rds"))
  possibleVars <- unique(c(possibleVars1, possibleVars2))


  misMatchNum <- 0

  varID <- checkPossCospotVarList(possibleVars, cSeq1, cSeq2, co_ID_1, co_ID_2, assayID, seq2search, wholeLeftWindow, wholeRightWindow)

  if (varID == -1) {
    # Read in saved list of possible mismatches and dedup
    possibleVars1 <- readRDS(paste0(cache$mismatchLibrary, "/", cSeq1, "_2.rds"))
    possibleVars2 <- readRDS(paste0(cache$mismatchLibrary, "/", cSeq2, "_2.rds"))
    possibleVars <- unique(c(possibleVars1, possibleVars2))

    varID <- checkPossCospotVarList(possibleVars, cSeq1, cSeq2, co_ID_1, co_ID_2, assayID, seq2search, wholeLeftWindow, wholeRightWindow)
    if (varID != -1) {
      misMatchNum <- 2
    }
  } else {
    misMatchNum <- 1
  }

  if (varID != -1) {
    updateDBMismatchCount(misMatchNum, varID, cache$db)
  }


  return(varID)
}

getDominantCospotList <- function(cospotKey) {
  domList <- c()
  seenC <- c()
  for (cRow in 1:nrow(cospotKey)) {
    if (cospotKey$Spot[cRow] %in% seenC) {
      next
    }
    dom <- min(c(cospotKey$Spot[cRow], cospotKey$Cospot[cRow]))
    domList <- c(domList, dom)
    seenC <- c(seenC, cospotKey$Spot[cRow], cospotKey$Cospot[cRow])
  }
  return(domList)
}
