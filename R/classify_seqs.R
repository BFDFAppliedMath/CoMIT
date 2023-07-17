#' COMIT Variant Classifier
#'
#' @param fileName File path to a FASTA containing sequences to sort
#' @param pullDate Date as a character vector that the FASTA was pulled from
#' GISAID,FORMAT: "YEAR-MONTH-DAY" i.e. "2021-03-05"
#' @param DATABASE_FILE File path to a CoMIT database
#' @param cladeDF Dataframe that can contain submission dates, lineages, and/or clades
#' @param refSeqFile A reference sequence for alignments. Default is the COVID Wuhan Sequence.
#' @param removeHighN Defaults to TRUE, will remove any sequences with more
#' than 1 percent N's
#' @param N_PER_FILTER Percentage of Ns to filer out
#' @param silently Run with out printing progress
#' @param leftSearchDist Integer, Left distance for buffer window region search
#' @param rightSearchDist Integer, Right distance for buffer window region search
#'
#' @return A vector with run log information
#'
comit_classify <- function(fileName, pullDate, DATABASE_FILE, cladeDF = NA,
                           refSeqFile = NULL, removeHighN = TRUE,
                           N_PER_FILTER = 1, silently = FALSE,
                           leftSearchDist = 400, rightSearchDist = 150) {
  cache$db <- DATABASE_FILE

  start_time <- Sys.time()
  print(start_time)

  # GLOBAL VARS------------------
  currSeq <- NULL
  accessionNum <- NULL
  currPrimerID <- NULL
  currVarId <- -1

  clearStrings()
  setStrings2AssayDB()
  setGISAIDstrings()
  setCacheConsts()

  # Build global data cache for use
  # cache$ primerDF, primerLocs, primerVars, cospotMap, assayKey, cospotPrimers, assayInfo
  # primerDF only has primerIDs for primers that aren't cospots
  libResult <- getDBLibrary(DATABASE_FILE)

  if ("error" %in% names(libResult)) {
    return(libResult$error)
  }

  # add left and right search distances to cache
  cache$leftSearchDist <- leftSearchDist
  cache$rightSearchDist <- rightSearchDist

  # Primer Lists for alignment classification

  # domCospots are the smaller primerID numbers in every pair
  domCospots <- getDominantCospotList(cache$cospotMap)

  primerList <- cache$primerDF[[s$PRIMER_ID]]
  primerList <- c(primerList, domCospots)

  # primerDF doesn't have cospot primers!
  # primerList has non-cospot primers + dominant cospot primers

  primerOrder <- 1:nrow(cache$primerDF)

  # Generate Primer Variant Files (Writes files)--------------------------------
  # Make mismatch directory
  cache$mismatchLibrary <- tempdir()

  generatePrimerVarLibrary(cache$primerDF)
  generatePrimerVarLibrary(cache$cospotPrimers)


  # FASTA Extract-----------------------------------
  if (is.null(refSeqFile)) {
    base <- system.file("cov_data", package = "CoMIT")
    refSeqFile <- dir(base, "refseq.fasta", full.names = TRUE)
  }

  fileResult <- getSeqs(fileName)
  # Stores seqs in cache$dnaSeqs
  refSeq <- seqinr::read.fasta(refSeqFile, as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)


  if ("error" %in% names(fileResult)) {
    return(fileResult$error)
  }

  # Get refSeqPrimerLocs table for amplicon tracking:
  # con <- DBI::dbConnect(RSQLite::SQLite(), DATABASE_FILE)
  # refSeqPrimerLocs <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_LOCS_TABLE))
  # DBI::dbDisconnect(con)
  # not needed--already in cache


  names <- attributes(cache$dnaSeqs)$names

  # Data Loading---------------

  counter <- 1
  nameNumber <- 1

  naAccCount <- 0
  hostError <- 0
  repeatedAccession <- 0

  if (!silently) {
    print("Running Seqs")
  }

  # For each Accession Loop----------------------------------------------
  for (accession in cache$dnaSeqs) {

    if (!silently & nameNumber %% 1000 == 0) {
      print(nameNumber)
    }

    # Parse Name for Location, Submission Date, and ID
    name <- names[nameNumber]
    nameNumber <- nameNumber + 1

    lineInfo <- parseFastaLine(name)
    if (!lineInfo$continue) {
      if (lineInfo$Fail_Type == "Host Error") {
        hostError <- hostError + 1
      } else if (lineInfo$Fail_Type == "NA Accession") {
        naAccCount <- naAccCount + 1
      }
      next
    }

    accessionNum <- lineInfo[["accNum"]]

    # Check if already in database
    if (isAccInDB(accessionNum, DATABASE_FILE)) {
      repeatedAccession <- repeatedAccession + 1
      # Next sequence
      next
    }

    counter <- counter + 1
    currSeq <- toupper(accession[1])


    N_Count <- round(stringr::str_count(currSeq, "N") / nchar(currSeq), 4) * 100
    if (removeHighN & N_Count >= N_PER_FILTER) {
      # If not already present, then add
      if (!isAccInHighNDB(accessionNum, DATABASE_FILE)) {
        insertToHighSeqsDB(DATABASE_FILE, accessionNum, N_Count)
      }
      # Go to next sequence
      next
    }

    # Add Accession Data to DATABASE---------------------
    if (is.data.frame(cladeDF)) {
      # Get Metadata
      accClade <- cladeDF[which(cladeDF[[s$METADATA_ACC_ID]] == accessionNum), ]
      metaVector <- getMetaData(accClade)

      insert2SeqInfo(
        DATABASE_FILE, accessionNum, lineInfo$location,
        lineInfo$colDate, lineInfo$colDay, name, pullDate,
        N_Count, metaVector[["clade"]], metaVector[["lineage"]]
      )
    } else {
      insert2SeqInfoNoMeta(
        DATABASE_FILE, accessionNum, lineInfo$location,
        lineInfo$colDate, lineInfo$colDay, name, pullDate, N_Count
      )
    }

    # Number of Variations on the Sequence
    cache$accessionVarCount <- 0

    # Create Assay List
    cache$AssayVarCounts <- NULL
    # Set Counts to 0
    for (id in unique(cache$assayKey[[s$ASSAY_ID]])) {
      cache$AssayVarCounts[[id]] <- 0
    }
    # Create tracking lists
    cache$AssayVarCountsWithAmb <- rlang::duplicate(cache$AssayVarCounts, shallow = FALSE)
    cache$AssayVarCounts10bp <- rlang::duplicate(cache$AssayVarCounts, shallow = FALSE)

    classifiedPrimerList <- c()
    accMaskedAlign <- NULL

    # # make new primer idx df
    # cache$primerIdxdf <- data.frame(Primer_ID=integer(), Assay_ID=integer(), Acc_Primer_Start=integer(), Acc_Primer_End=integer())


    # For each Primer
    for (primerRow in primerOrder) {
      currPrimerID <- cache$primerDF[[s$PRIMER_ID]][primerRow]
      currAssayID <- cache$primerDF[[s$ASSAY_ID]][primerRow]
      primerSeq <- cache$primerDF[[s$SEARCHABLE_SEQ]][primerRow]

      primerSubSeqRes <- getPrimerSubSeq(currSeq, currPrimerID)

      currSearchSeq4Primer <- primerSubSeqRes$subSeq
      wholeLeftWindow <- primerSubSeqRes$wholeLeftWindow
      wholeRightWindow <- primerSubSeqRes$wholeRightWindow

      # If primer found in sequence (perfect match)

      if (any(grep(primerSeq, currSearchSeq4Primer, fixed = TRUE))) {
        classifiedPrimerList <- c(classifiedPrimerList, currPrimerID)
        next # primer in seq
      }

      # # Version with primer idx df:
      # if(findPrimerAddLoc(currPrimerID, currAssayID, primerSeq, currSearchSeq4Primer, wholeLeftWindow, wholeRightWindow) != -1) {
      #   classifiedPrimerList <- c(classifiedPrimerList, currPrimerID)
      #   next #primer in seq
      # }

      currVarId <- -1

      # VARIANT Identification------------------
      currVarId <- isKnownVariant(currPrimerID, currAssayID, primerSeq, currSearchSeq4Primer, wholeLeftWindow, wholeRightWindow)
      if (currVarId == -1) {
        currVarId <- findNovelVar(primerSeq, currPrimerID, currAssayID, currSearchSeq4Primer, wholeLeftWindow, wholeRightWindow) # checks for all possible 1 and 2 mismatch versions of this primer in db

        if (currVarId == -1) {
          if (is.null(accMaskedAlign)) {
            # Align for indels, ambiguities, and >2 mismatches
            res <- getMaskedAlignmentList(
              currSeq, accessionNum, refSeq,
              setdiff(primerList, classifiedPrimerList)
            )
            accMaskedAlign <- res$alignment
            primerLocdf <- res$locs
            rm(res)
          }
          # Leverage Alignment
          classifyFromAlignment(accMaskedAlign, currPrimerID, accessionNum, DATABASE_FILE,
            cospotList = domCospots
          )

          # # add primer idx info for current primer
          # cache$primerIdxdf <- cache$primerIdxdf %>% rbind(dplyr::filter(primerLocdf, get(s$PRIMER_ID)==currPrimerID))

          next # primer in seq
        }
      }

      add2IdMapTable(accessionNum, currVarId, currPrimerID, getAssayFromPrimer(currPrimerID), DATABASE_FILE)

      # print(paste(currPrimerID, currVarId))
      updateVarCounts(getAssayFromPrimer(currPrimerID), isWithin10bp = isVarId3Prime(currVarId, DATABASE_FILE))

      classifiedPrimerList <- c(classifiedPrimerList, currPrimerID)
    } # end of primer loop



    cospotSeenList <- c()

    # COSPOTS
    for (primer in domCospots) {
      co_1_ID <- primer
      primerRow <- which(cache$cospotPrimers[[s$PRIMER_ID]] == co_1_ID, arr.ind = TRUE)
      co_1_seq <- cache$cospotPrimers[[s$SEARCHABLE_SEQ]][primerRow]

      currAssayID <- cache$cospotPrimers[[s$ASSAY_ID]][primerRow]

      # Cospot
      cospotPair <- cache$cospotMap[which(cache$cospotMap[[s$SPOT1]] == co_1_ID), ]
      co_2_ID <- cospotPair[[s$SPOT2]][1]
      coRow <- which(cache$cospotPrimers[[s$PRIMER_ID]] == co_2_ID, arr.ind = TRUE)


      # coDF <- filter(cache$cospotPrimers, Primer_ID == co_2_ID)

      co_2_seq <- cache$cospotPrimers[[s$SEARCHABLE_SEQ]][coRow]

      # primerSubSeqRes <- getPrimerSubSeq(currSeq, co_1_ID, rightSearchDist)
      primerSubSeqRes <- getPrimerSubSeq(currSeq, co_1_ID)

      currSearchSeq4Primer <- primerSubSeqRes$subSeq
      wholeLeftWindow <- primerSubSeqRes$wholeLeftWindow
      wholeRightWindow <- primerSubSeqRes$wholeRightWindow

      # If primer found in sequence (perfect match)
      if (any(grep(co_1_seq, currSearchSeq4Primer, fixed = TRUE))) {
        # if(findPrimerAddLoc(co_1_ID, currAssayID, co_1_seq, currSearchSeq4Primer, wholeLeftWindow, wholeRightWindow) != -1) {
        # Add to classified list
        classifiedPrimerList <- c(classifiedPrimerList, co_1_ID, co_2_ID)
        next # primer in seq
      }

      if (any(grep(co_2_seq, currSearchSeq4Primer, fixed = TRUE))) {
        # if(findPrimerAddLoc(co_2_ID, currAssayID, co_2_seq, currSearchSeq4Primer, wholeLeftWindow, wholeRightWindow) != -1) {
        # Add to classified list
        classifiedPrimerList <- c(classifiedPrimerList, co_1_ID, co_2_ID)
        next # primer in seq
      }

      currVarId <- -1

      # VARIANT Identification------------------
      co_var_1_ID <- -1
      co_var_2_ID <- -1

      co_var_1_ID <- isKnownVariant(co_1_ID, currAssayID, co_1_seq, currSearchSeq4Primer, wholeLeftWindow, wholeRightWindow)
      if (co_var_1_ID != -1) {
        currVarId <- co_var_1_ID
        classifiedPrimerList <- c(classifiedPrimerList, co_1_ID, co_2_ID)
      } else {
        co_var_2_ID <- isKnownVariant(co_2_ID, currAssayID, co_2_seq, currSearchSeq4Primer, wholeLeftWindow, wholeRightWindow)
        if (co_var_2_ID != -1) {
          currVarId <- co_var_2_ID
          classifiedPrimerList <- c(classifiedPrimerList, co_1_ID, co_2_ID)
        }
      }

      # If not found in variant list
      if (currVarId == -1) {
        currVarId <- findNovelCospotVar(co_1_seq, co_2_seq, co_1_ID, co_2_ID, currAssayID, currSearchSeq4Primer, wholeLeftWindow, wholeRightWindow)

        if (currVarId == -1) {
          if (is.null(accMaskedAlign)) {
            # Align for indels, ambiguities, and >2 mismatches
            res <- getMaskedAlignmentList(
              currSeq, accessionNum, refSeq,
              setdiff(primerList, classifiedPrimerList)
            )
            accMaskedAlign <- res$alignment
            primerLocdf <- res$locs
            # print("here's primerLocdf after cospot alignment:")
            rm(res)
          }
          # Leverage Alignment
          classifyFromAlignment(accMaskedAlign, co_1_ID, accessionNum, DATABASE_FILE,
            cospotList = domCospots
          )

          # # add primer idx info for current primer
          # cache$primerIdxdf <- cache$primerIdxdf %>% rbind(dplyr::filter(primerLocdf, get(s$PRIMER_ID)==currPrimerID))

          # Move to next primerID
          next # dominant primer in seq
        } else {
          classifiedPrimerList <- c(classifiedPrimerList, co_1_ID, co_2_ID)
        }
      }

      primerR <- cache$primerVars %>% dplyr::filter(get(s$VAR_ID) == currVarId)
      coSpotPrimerUsed <- primerR[[s$PRIMER_ID]][1]
      add2IdMapTable(accessionNum, currVarId, coSpotPrimerUsed, getAssayFromPrimer(coSpotPrimerUsed), DATABASE_FILE)

      # Variant and Assay Count
      updateVarCounts(getAssayFromPrimer(coSpotPrimerUsed), isAmbiguity = FALSE, isWithin10bp = isVarId3Prime(currVarId, DATABASE_FILE))
    } # end dominant primer loop



    # Total assay var
    assayVarCount <- 0
    AVCwithAmb <- 0
    AVC10bp <- 0

    for (id in unique(cache$assayKey[[s$ASSAY_ID]])) {
      if (cache$AssayVarCounts[[id]] >= 1) {
        assayVarCount <- assayVarCount + 1
      }
      if (cache$AssayVarCountsWithAmb[[id]] >= 1) {
        AVCwithAmb <- AVCwithAmb + 1
      }
      if (cache$AssayVarCounts10bp[[id]] >= 1) {
        AVC10bp <- AVC10bp + 1
      }
    }

    updateDBvarAndAssayCount(accessionNum, cache$accessionVarCount, assayVarCount, AVCwithAmb, AVC10bp, DATABASE_FILE)

    # analyze_amplicons(accessionNum, currSeq, refSeq, primerList)
  } # end of looping accessions

  cache$dnaSeqs <- NULL

  # Build Database Run output files--------------------------------------------

  fasta_p <- paste("Ran file", fileName)

  # SEQUENCES
  fasta_seqs <- paste("Total Sequences in FASTA:", length(names))
  total_seqs <- paste("Total Sequences Run:", counter - 1)
  cache$TotalSeqsRun <- cache$TotalSeqsRun + counter - 1
  seq_sem_accession <- paste("Number of Seqs without accession number:", naAccCount)
  host_error_state <- paste("Number of Seqs with wrong host:", hostError)

  # TIME
  end_time <- Sys.time()
  time_statment <- paste("Start:", start_time, "End:", end_time)
  totalRunTime <- round(difftime(end_time, start_time, units = "mins")[[1]], 3)
  total_time <- paste("Total Run Time:", totalRunTime, "mins")

  # PRINT
  print(total_seqs)
  print(total_time)
  classifyRate <- round((counter - 1) / totalRunTime, 2)
  print(paste("Sequence Classification Rate:", classifyRate, "seqs/min"))

  fileBuilder <- c(
    fasta_p, fasta_seqs, total_seqs, " ", seq_sem_accession, host_error_state, " ",
    time_statment, total_time, " ", " "
  )

  return(fileBuilder)
}
