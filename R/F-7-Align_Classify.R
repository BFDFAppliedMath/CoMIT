#' Get Masked Alignment List from supplied primers
#'
#' @param seq Character vector, The current sequence being analyzed
#' @param accID, Accession Number
#' @param refSeq Reference Sequence to align to
#' @param primersNotClassified Vector of Primers still needed
#'
#' @return List containing masked alignment sections for each primer
#' @noRd
#'
getMaskedAlignmentList <- function(seq, accID, refSeq, primersNotClassified) { # TODO: we're passing primers we've found perfect matches for for this sequence as primersnotClassified, is that right??

  # Pairwise Align and mask

  # Run Alignment-------------------------------
  seq <- gsub("\\s", "", seq)
  seq <- gsub("U", "T", seq)
  pair <- c(refSeq[[1]][1], seq)
  names(pair) <- c("Reference", accID)
  dnaSet <- Biostrings::DNAStringSet(pair)
  names(dnaSet) <- names(pair)

  dna1 <- DECIPHER::RemoveGaps(dnaSet, removeGaps = "all", processors = NULL)

  dna1 <- DECIPHER::OrientNucleotides(dna1, reference = 1, verbose = FALSE, processors = NULL)
  alignmentO <- DECIPHER::AlignSeqs(dna1, processors = NULL, iterations = 1, verbose = FALSE)

  # Get consensus sequence positions----------------------------------------
  MSA <- Biostrings::DNAMultipleAlignment(alignmentO) # convert alignmentO to DNAMultipleAlignment object
  alignedRef <- toString(alignmentO$Reference)

  # Create Ref Map-------------
  insV <- unlist(strsplit(alignedRef, ""))
  refMap <- rep(-1, nchar(refSeq)) # map from reference seq to alignment seq
  insPos <- 1
  refPos <- 1

  for (i in insV) {
    if (i != "-") {
      refMap[refPos] <- insPos
      refPos <- refPos + 1
    }
    # Increment up one
    insPos <- insPos + 1
  }



  # Create Seq Map (for amplicon information)-------------
  insV <- unlist(strsplit(toString(alignmentO[[2]]), "")) # alignment
  varSeqMap <- rep(-1, length(insV)) # map from alignment to variant seq
  alignPos <- 1
  varSeqPos <- 1

  for (i in insV) {
    if (i != "-") {
      varSeqMap[alignPos] <- varSeqPos
      varSeqPos <- varSeqPos + 1
    }
    # Increment up one
    alignPos <- alignPos + 1
  }





  # Find Locations-------------
  primerStart <- list()
  primerWidth <- list()
  iStart <- c()
  iWidth <- c()
  iNames <- c()

  unclassifiedPrimerLocs <- data.frame(primerID = integer(), assayID = integer(), startIdx = integer(), endIdx = integer())


  for (notfoundPrimerID in primersNotClassified) {
    primerLoc <- cache$primerLocs[which(cache$primerLocs$Primer_ID == notfoundPrimerID), ] # get original ref seq location
    # use map to get aligned positions
    startR <- refMap[primerLoc$Start[1]]
    endR <- refMap[primerLoc$End[1]]

    iStart <- c(iStart, startR)
    iWidth <- c(iWidth, (endR - startR + 1))
    iNames <- c(iNames, notfoundPrimerID)

    startV <- varSeqMap[startR] # start pos of variant seq
    endV <- varSeqMap[endR]

    unclassifiedPrimerLocs <- unclassifiedPrimerLocs %>% rbind(data.frame(
      Primer_ID = notfoundPrimerID,
      Assay_ID = getAssayFromPrimer(notfoundPrimerID),
      Acc_Primer_Start = startV,
      Acc_Primer_End = endV
    ))
  }

  primer_iRanges <- IRanges::IRanges(start = iStart, width = iWidth)
  names(primer_iRanges) <- iNames

  primerAlignment <- list()
  for (primerRow in 1:length(iStart)) {
    align_mask <- MSA
    Biostrings::colmask(align_mask, invert = TRUE) <- primer_iRanges[primerRow]
    primerAlignment[[names(primer_iRanges)[primerRow]]] <- methods::as(align_mask, "DNAStringSet")
  }

  return(list(alignment = primerAlignment, locs = unclassifiedPrimerLocs))
}


#' Classify Variant from a masked alignment
#'
#' @param primerAlignment masked alignment
#' @param currPrimerID ID for the primer to classify
#' @param accID Accession Number for the sequence
#' @param DATABASE_FILE Database to add to
#' @param cospotList List of CoSpot Primers
#' @param saveIndelAligns Boolean that indicates whether the aligments should be saved to the server for later viewing
#' @return a list containing the closest cospot ID for the variant, the variant seq, and whether or not it contains an indel
#' @noRd
#'
classifyFromAlignment <- function(primerAlignment, currPrimerID, accID, DATABASE_FILE, cospotList = NULL, saveIndelAligns = FALSE) {
  # Classify Primers-----------------------------------------------

  # Get Variant String
  alignment <- primerAlignment[[toString(currPrimerID)]]
  varPrimerSeq <- alignment[[toString(accID)]]
  variantPrimerString <- toString(varPrimerSeq)

  # Get Orig
  refPrimerSeq <- alignment[["Reference"]]
  refPrimerString <- toString(refPrimerSeq)

  hasIndels <- FALSE

  if ((!is.null(cospotList)) & (currPrimerID %in% cospotList)) {
    # Before sorting test # of char difference between cospot 1 and 2 if cospot and sort whichever is closer.
    coSpotMatchID <- getClosestCoSpot(currPrimerID, variantPrimerString)
    if (coSpotMatchID[["indel"]] == FALSE) {
      hasIndels <- sortAlignmentPrimer(accID, coSpotMatchID[["id"]], coSpotMatchID[["seq"]], variantPrimerString, DATABASE_FILE, isCospot = TRUE)
    } else {
      # Indels will always use the dominant cospot
      hasIndels <- sortAlignmentPrimer(accID, currPrimerID, refPrimerString, variantPrimerString, DATABASE_FILE, isCospot = TRUE)
    }
  } else {
    hasIndels <- sortAlignmentPrimer(accID, currPrimerID, refPrimerString, variantPrimerString, DATABASE_FILE, isCospot = FALSE)
  }
}

getClosestCoSpot <- function(primerID, variantSeq) {
  co_1_ID <- primerID
  co1_row <- cache$cospotPrimers[which(cache$cospotPrimers[[s$PRIMER_ID]] == co_1_ID), ]
  co_1_seq <- co1_row[[s$SEARCHABLE_SEQ]][1]

  currAssayID <- co1_row[[s$ASSAY_ID]][1]

  # Cospot
  cospotPair <- cache$cospotMap[which(cache$cospotMap[[s$SPOT1]] == co_1_ID), ]
  co_2_ID <- cospotPair[[s$SPOT2]][1]
  co2_row <- cache$cospotPrimers[which(cache$cospotPrimers[[s$PRIMER_ID]] == co_2_ID), ]
  co_2_seq <- co2_row[[s$SEARCHABLE_SEQ]][1]

  co1List <- unlist(strsplit(co_1_seq, ""))
  co2List <- unlist(strsplit(co_2_seq, ""))
  varList <- unlist(strsplit(variantSeq, ""))

  # If they are different lengths it will have to be the one that matches the reference sequence
  if (length(co1List) != length(varList)) {
    return(list("id" = min(c(co_1_ID, co_2_ID)), "indel" = TRUE))
  }
  diffCount1 <- 0
  diffCount2 <- 0

  # Which primer is closer
  for (i in 1:nchar(variantSeq)) {
    if (co1List[i] != varList[i]) {
      diffCount1 <- diffCount1 + 1
    }

    if (co2List[i] != varList[i]) {
      diffCount2 <- diffCount2 + 1
    }
  }
  if (diffCount2 < diffCount1) {
    closestCospotID <- co_2_ID
    seq <- co_2_seq
  } else {
    closestCospotID <- co_1_ID
    seq <- co_1_seq
  }
  return(list("id" = closestCospotID, "seq" = seq, "indel" = FALSE))
}
