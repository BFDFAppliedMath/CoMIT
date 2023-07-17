# Create Test fasta

create_test_fasta_and_DB <- function() {
  # Build DB
  base <- system.file("sql", package = "CoMIT")
  schemaFile <- dir(base, "base_comit_schema.sql", full.names = TRUE)

  synDir <- tempdir()

  DB <- paste0(synDir, "/FullFake_V2.DB")

  status <- buildDBFromSchema(schemaFile, DB,
    force = TRUE, silently = TRUE
  )

  set.seed(5)
  DNA <- c("A", "G", "T", "C")
  chain <- sample(DNA, 500, replace = T)
  refseq <- paste(chain, collapse = "")


  refseqFile <- paste0(synDir, "/gen_ref_seq.fasta")
  seqinr::write.fasta(refseq, "Reference", refseqFile)

  buildPrimers(refseq, DB)


  # Pull from Database------------------------------
  con <- DBI::dbConnect(RSQLite::SQLite(), DB)

  primerDF <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_TABLE))
  primerLocs <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_LOCS_TABLE))

  DBI::dbDisconnect(con)

  mutDNA <- refseq

  seqNum <- 1
  seqList <- c(refseq)
  mutationTracker <- tibble::tibble(
    Seq_Num = rep(1, 4),
    Primer_Num = 1:4,
    M_V = rep("M", 4),
    Mismatch = NA,
    Position = NA,
    Original = NA
  )

  # MUTATE SEQS-----------------------------------------
  seqNum <- seqNum + 1
  mutDNA <- refseq
  # Beginning single
  # print("Beginning Single")
  for (i in 1:(nrow(primerLocs) - 1)) {
    mV <- mutateAtPrimer(i, mutDNA, primerLocs, primerDF, mutLoc = "S")
    mutDNA <- unname(mV["mSeq"])
    mutationTracker <- rbind(mutationTracker, c(seqNum, i, "V", mV["MM"], mV["mPos"], mV["Orig"]))
  }
  seqList <- c(seqList, mutDNA)

  seqNum <- seqNum + 1
  mutDNA <- refseq
  # End single
  # print("End Single")
  for (i in 1:(nrow(primerLocs) - 1)) {
    mV <- mutateAtPrimer(i, mutDNA, primerLocs, primerDF, mutLoc = "E")
    mutDNA <- unname(mV["mSeq"])
    mutationTracker <- rbind(mutationTracker, c(seqNum, i, "V", mV["MM"], mV["mPos"], mV["Orig"]))
  }
  seqList <- c(seqList, mutDNA)

  seqNum <- seqNum + 1
  mutDNA <- refseq
  # Middle
  # print("Random Single")
  for (i in 1:(nrow(primerLocs) - 1)) {
    mV <- mutateAtPrimer(i, mutDNA, primerLocs, primerDF)
    mutDNA <- unname(mV["mSeq"])
    mutationTracker <- rbind(mutationTracker, c(seqNum, i, "V", mV["MM"], mV["mPos"], mV["Orig"]))
  }
  seqList <- c(seqList, mutDNA)

  # Random Primer Mutation Loop
  for (j in 1:25) {
    seqNum <- seqNum + 1
    mutDNA <- refseq

    for (i in 1:(nrow(primerLocs) - 1)) {
      isCospot <- FALSE

      # If primer 4
      if (i == 4) {
        # Randomization to choose cospot 1/3 of the time
        n <- sample(1:3, 1)
        if (n == 1) {
          b <- primerLocs[[s$LOC_START]][5]
          e <- primerLocs[[s$LOC_END]][5]
          substr(mutDNA, b, e) <- primerDF[[s$SEARCHABLE_SEQ]][5]
          isCospot <- TRUE
          i <- 5
        }
      }

      # Randomization to mutate 25% of the time
      n <- sample(1:4, 1)
      if (n == 1) {
        mV <- mutateAtPrimer(i, mutDNA, primerLocs, primerDF)
        mutDNA <- unname(mV["mSeq"])

        if (isCospot) {
          mutationTracker <- rbind(
            mutationTracker,
            c(seqNum, i, "VC", mV["MM"], mV["mPos"], mV["Orig"])
          )
        } else {
          mutationTracker <- rbind(
            mutationTracker,
            c(seqNum, i, "V", mV["MM"], mV["mPos"], mV["Orig"])
          )
        }
      } else {
        if (isCospot) {
          mutationTracker <- rbind(
            mutationTracker,
            c(seqNum, i, "MC", NA, NA, NA)
          )
        } else {
          mutationTracker <- rbind(
            mutationTracker,
            c(seqNum, i, "M", NA, NA, NA)
          )
        }
      }
    }
    seqList <- c(seqList, mutDNA)
  }


  # Containing ambiguities and double mismatches
  # print("Random Double Am")
  for (j in 1:25) {
    seqNum <- seqNum + 1
    mutDNA <- refseq

    for (i in 1:(nrow(primerLocs) - 1)) {
      isCospot <- FALSE
      firstAmb <- FALSE
      # if primer 4
      if (i == 4) {
        # Randomization to choose cospot 1/3 of the time
        n <- sample(1:3, 1)
        if (n == 1) {
          b <- primerLocs[[s$LOC_START]][5]
          e <- primerLocs[[s$LOC_END]][5]
          substr(mutDNA, b, e) <- primerDF[[s$SEARCHABLE_SEQ]][5]
          isCospot <- TRUE
          i <- 5
        }
      }

      # Randomization to mutate 25% of the time
      n <- sample(1:4, 1)
      if (n == 1) {
        # Randomization for 2 mismatches 50% of the time
        n <- sample(1:2, 1)
        if (n == 1) {
          # Randomization for ambiguities 33% of the time
          n <- sample(1:3, 1)
          if (n == 1) {
            mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = TRUE)
            hasAmbiguity <- TRUE
          } else {
            mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF)
            hasAmbiguity <- FALSE
          }
          mutDNA <- unname(mV["mSeq"])
        } else {
          # TWO MISMATCHES

          # Randomization for ambiguities 25% of the time
          n <- sample(1:4, 1)
          if (n == 1) {
            mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = TRUE)
            firstAmb <- TRUE
            hasAmbiguity <- TRUE
          } else {
            mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = FALSE)
            hasAmbiguity <- FALSE
          }
          mutDNA <- unname(mV["mSeq"])

          if (hasAmbiguity) {
            varType <- "A"
          } else {
            if (isCospot) {
              varType <- "VC"
            } else {
              varType <- "V"
            }
          }
          mutationTracker <- rbind(
            mutationTracker,
            c(seqNum, i, varType, mV["MM"], mV["mPos"], mV["Orig"])
          )

          # Mutation 2
          # Randomization for ambiguities 25% of the time
          n <- sample(1:4, 1)
          if (n == 1) {
            mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = TRUE)
            hasAmbiguity <- TRUE
          } else {
            mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = FALSE)
            hasAmbiguity <- FALSE
          }
          mutDNA <- unname(mV["mSeq"])
        }



        if (hasAmbiguity | firstAmb) {
          varType <- "A"
        } else {
          if (isCospot) {
            varType <- "VC"
          } else {
            varType <- "V"
          }
        }
        mutationTracker <- rbind(
          mutationTracker,
          c(seqNum, i, varType, mV["MM"], mV["mPos"], mV["Orig"])
        )
      } else {
        if (isCospot) {
          mutationTracker <- rbind(
            mutationTracker,
            c(seqNum, i, "MC", NA, NA, NA)
          )
        } else {
          mutationTracker <- rbind(
            mutationTracker,
            c(seqNum, i, "M", NA, NA, NA)
          )
        }
      }
    }
    seqList <- c(seqList, mutDNA)
  }

  seqNames <- generateSeqNames(length(seqList) + 1)
  fastaFile <- paste0(synDir, "/gen_data.fasta")
  seqinr::write.fasta(as.list(seqList), seqNames, fastaFile)

  return(list(
    "FASTA" = fastaFile,
    "MutDF" = mutationTracker,
    "DatabaseFile" = DB,
    "SeqNames" = seqNames,
    "refFile" = refseqFile
  ))
}
