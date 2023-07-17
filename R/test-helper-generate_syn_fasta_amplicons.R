# Create Test fasta

create_test_fasta_and_DB_amplicons <- function() {
  # cache$TOO_BIG_INNER_AMPLICON = 50 # idk...is there a set inner/outer amplicon length that messes up Tm that applies to all sequences? or is it a proportional thing??
  # cache$TOO_BIG_OUTER_AMPLICON = 200
  #
  # cache$TOO_SMALL_INNER_AMPLICON = 5
  # cache$TOO_SMALL_OUTER_AMPLICON = 50
  #
  # cache$TOO_MUCH_GC = .9 # 90%??
  # cache$TOO_LITTLE_GC = .1 # 10%??

  # Build DB
  base <- system.file("sql", package = "CoMIT")
  schemaFile <- dir(base, "base_comit_schema.sql", full.names = TRUE)

  synDir <- tempdir()

  DB <- paste0(synDir, "/FullFake_AmpliconChanges.DB")

  status <- buildDBFromSchema(schemaFile, DB,
    force = TRUE, silently = TRUE
  )

  set.seed(5)
  DNA <- c("A", "G", "T", "C")
  chain <- sample(DNA, 1200, replace = T)
  DNAs <- paste(chain, collapse = "")


  refseqFile <- paste0(synDir, "/gen_ref_seq.fasta")
  seqinr::write.fasta(DNAs, "Reference", refseqFile)

  buildPrimers(DNAs, DB, multipleAssays = TRUE)

  # addAmpliconInfoToDB(refseqFile, DB) # TODO: first it needs other assay info...

  # Pull from Database------------------------------
  con <- DBI::dbConnect(RSQLite::SQLite(), DB)

  primerDF <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_TABLE))
  primerLocs <- DBI::dbGetQuery(conn = con, genSelect(s$PRIMER_LOCS_TABLE))

  DBI::dbDisconnect(con)

  mutDNA <- DNAs # set original mutDNA to ref seq

  seqNum <- 1
  seqList <- c(DNAs) # All perfect matches (original seqList is set to ref seq)
  mutationTracker <- tibble::tibble(
    Seq_Num = 1,
    og_inner_length = NA,
    mut_inner_length = NA,
    og_outer_length = NA,
    mut_outer_length = NA,
    og_GC = NA,
    mut_GC = NA
  )
  # Primer_Num = 1:9,
  # M_V = rep("M", 9),
  # Mismatch = NA,
  # Position = NA,
  # Original = NA)


  # mutationTracker["mutSeq"]  entire mutated sequence
  # mutationInfo["ogInnerLength"],
  # mutationInfo["mutInnerLength"],
  # mutationInfo["ogOuterLength"],
  # mutationInfo["mutOuterLength"],
  # mutationInfo["ogGC"],
  # mutationInfo["mutGC"]

  # MUTATE SEQS-----------------------------------------

  # increase GC content in inner amplicon to > threshold
  # insertions inside inner amplicon that make inner length > max threshold
  #  + increase GC
  # insertions inside outer amplicon that make outer length > max threshold
  #  + increase GC
  # insertions inside both amplicons
  #  + increase GC
  # deletions inside inner amplicon that make inner length < min threshold
  #  + increase GC
  # deletions inside outer amplicon that make outer length < min threshold
  #  + increase GC
  # deletions inside both amplicons
  #  + increase GC

  # some combos




  # just GC ----------------------------------------

  # next seq
  seqNum <- seqNum + 1
  mutDNA <- DNAs
  # increase GC content in inner amplicon to > max threshold
  mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID = 1, indel = NA, indelRegionNums = NA, changeGC = "inc")
  mutationTracker <- rbind(mutationTracker, c(
    seqNum, mutationInfo["ogInnerLength"],
    mutationInfo["mutInnerLength"],
    mutationInfo["ogOuterLength"],
    mutationInfo["mutOuterLength"],
    mutationInfo["ogGC"],
    mutationInfo["mutGC"]
  ))
  mutDNA <- unname(mutationInfo["mutSeq"])
  seqList <- c(seqList, mutDNA)

  # mutDNA <- unname(mV["mutSeq"]) # set mutDNA to the new mutated seq
  # mutationTracker  <- rbind(mutationTracker, c(seqNum, i, "V", mV["MM"], mV["mPos"], mV["Orig"])) # r bind every element of the returned vector except for the mutated seq
  # seqList = c(seqList, mutDNA) # add the mutated seq to the seqList



  # next seq
  seqNum <- seqNum + 1
  mutDNA <- DNAs
  # decrease GC content in inner amplicon to < min threshold
  mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID = 2, indel = NA, indelRegionNums = NA, changeGC = "dec")
  mutationTracker <- rbind(mutationTracker, c(
    seqNum, mutationInfo["ogInnerLength"],
    mutationInfo["mutInnerLength"],
    mutationInfo["ogOuterLength"],
    mutationInfo["mutOuterLength"],
    mutationInfo["ogGC"],
    mutationInfo["mutGC"]
  ))
  mutDNA <- unname(mutationInfo["mutSeq"])
  seqList <- c(seqList, mutDNA)



  # insertions -----------------------------------------

  # next seq
  seqNum <- seqNum + 1
  mutDNA <- DNAs
  # insertions inside inner amplicon that make inner length > max threshold
  mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID = 1, indel = "ins", indelRegionNums = 2, changeGC = NA)
  mutationTracker <- rbind(mutationTracker, c(
    seqNum, mutationInfo["ogInnerLength"],
    mutationInfo["mutInnerLength"],
    mutationInfo["ogOuterLength"],
    mutationInfo["mutOuterLength"],
    mutationInfo["ogGC"],
    mutationInfo["mutGC"]
  ))
  mutDNA <- unname(mutationInfo["mutSeq"])
  seqList <- c(seqList, mutDNA)



  # # next seq
  # seqNum = seqNum + 1
  # mutDNA = DNAs
  # # insertions inside inner amplicon that make inner length > max threshold + increase inner GC content
  # mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID=2, indel="ins", indelRegionNums=2, changeGC="inc")
  # mutationTracker <- rbind(mutationTracker, c(seqNum, mutationInfo["ogInnerLength"],
  #                                             mutationInfo["mutInnerLength"],
  #                                             mutationInfo["ogOuterLength"],
  #                                             mutationInfo["mutOuterLength"],
  #                                             mutationInfo["ogGC"],
  #                                             mutationInfo["mutGC"]))
  # mutDNA <- unname(mutationTracker["mutSeq"])
  # seqList <- c(seqList, mutDNA)



  # next seq
  seqNum <- seqNum + 1
  mutDNA <- DNAs
  # insertions inside outer amplicon that make outer length > max threshold
  mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID = 1, indel = "ins", indelRegionNums = 1, changeGC = NA)
  mutationTracker <- rbind(mutationTracker, c(
    seqNum, mutationInfo["ogInnerLength"],
    mutationInfo["mutInnerLength"],
    mutationInfo["ogOuterLength"],
    mutationInfo["mutOuterLength"],
    mutationInfo["ogGC"],
    mutationInfo["mutGC"]
  ))
  mutDNA <- unname(mutationInfo["mutSeq"])
  seqList <- c(seqList, mutDNA)



  # # next seq
  # seqNum = seqNum + 1
  # mutDNA = DNAs
  # # insertions inside outer amplicon that make outer length > max threshold + increase inner GC content
  # mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID=2, indel="ins", indelRegionNums=1, changeGC="inc")
  # mutationTracker <- rbind(mutationTracker, c(seqNum, mutationInfo["ogInnerLength"],
  #                                             mutationInfo["mutInnerLength"],
  #                                             mutationInfo["ogOuterLength"],
  #                                             mutationInfo["mutOuterLength"],
  #                                             mutationInfo["ogGC"],
  #                                             mutationInfo["mutGC"]))
  # mutDNA <- unname(mutationTracker["mutSeq"])
  # seqList <- c(seqList, mutDNA)



  # # next seq
  # seqNum = seqNum + 1
  # mutDNA = DNAs
  # # insertions inside both amplicons
  # mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID=1, indel="ins", indelRegionNums=c(2,3), changeGC=NA)
  # mutationTracker <- rbind(mutationTracker, c(seqNum, mutationInfo["ogInnerLength"],
  #                                             mutationInfo["mutInnerLength"],
  #                                             mutationInfo["ogOuterLength"],
  #                                             mutationInfo["mutOuterLength"],
  #                                             mutationInfo["ogGC"],
  #                                             mutationInfo["mutGC"]))
  # mutDNA <- unname(mutationTracker["mutSeq"])
  # seqList <- c(seqList, mutDNA)



  # # next seq
  # seqNum = seqNum + 1
  # mutDNA = DNAs
  # # insertions inside both amplicons + increase inner GC content
  # mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID=2, indel="ins", indelRegionNums=c(2,3), changeGC="inc")
  # mutationTracker <- rbind(mutationTracker, c(seqNum, mutationInfo["ogInnerLength"],
  #                                             mutationInfo["mutInnerLength"],
  #                                             mutationInfo["ogOuterLength"],
  #                                             mutationInfo["mutOuterLength"],
  #                                             mutationInfo["ogGC"],
  #                                             mutationInfo["mutGC"]))
  # mutDNA <- unname(mutationTracker["mutSeq"])
  # seqList <- c(seqList, mutDNA)




  # deletions--------------------------------------------------

  # next seq
  seqNum <- seqNum + 1
  mutDNA <- DNAs
  # deletions inside inner amplicon that make inner length < min threshold
  mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID = 1, indel = "del", indelRegionNums = 2, changeGC = NA)
  mutationTracker <- rbind(mutationTracker, c(
    seqNum, mutationInfo["ogInnerLength"],
    mutationInfo["mutInnerLength"],
    mutationInfo["ogOuterLength"],
    mutationInfo["mutOuterLength"],
    mutationInfo["ogGC"],
    mutationInfo["mutGC"]
  ))
  mutDNA <- unname(mutationInfo["mutSeq"])
  seqList <- c(seqList, mutDNA)



  # next seq
  seqNum <- seqNum + 1
  mutDNA <- DNAs
  # deletions inside outer amplicon that make outer length < min threshold
  mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID = 1, indel = "del", indelRegionNums = 3, changeGC = NA)
  mutationTracker <- rbind(mutationTracker, c(
    seqNum, mutationInfo["ogInnerLength"],
    mutationInfo["mutInnerLength"],
    mutationInfo["ogOuterLength"],
    mutationInfo["mutOuterLength"],
    mutationInfo["ogGC"],
    mutationInfo["mutGC"]
  ))
  mutDNA <- unname(mutationInfo["mutSeq"])
  seqList <- c(seqList, mutDNA)



  # # next seq
  # seqNum = seqNum + 1
  # mutDNA = DNAs
  # # deletions inside inner amplicon that make inner length < min threshold + increase inner GC content
  # mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID=2, indel="del", indelRegionNums=2, changeGC="inc")
  # mutationTracker <- rbind(mutationTracker, c(seqNum, mutationInfo["ogInnerLength"],
  #                                             mutationInfo["mutInnerLength"],
  #                                             mutationInfo["ogOuterLength"],
  #                                             mutationInfo["mutOuterLength"],
  #                                             mutationInfo["ogGC"],
  #                                             mutationInfo["mutGC"]))
  # mutDNA <- unname(mutationTracker["mutSeq"])
  # seqList <- c(seqList, mutDNA)





  # GC content decrease


  # # next seq
  # seqNum = seqNum + 1
  # mutDNA = DNAs
  # # insertions inside both amplicons + decrease inner GC content
  # mutationInfo <- mutateAtAmplicon(mutDNA, primerLocs, primerDF, assayID=1, indel="ins", indelRegionNums=c(2,3), changeGC="dec")
  # mutationTracker <- rbind(mutationTracker, c(seqNum, mutationInfo["ogInnerLength"],
  #                                             mutationInfo["mutInnerLength"],
  #                                             mutationInfo["ogOuterLength"],
  #                                             mutationInfo["mutOuterLength"],
  #                                             mutationInfo["ogGC"],
  #                                             mutationInfo["mutGC"]))
  # mutDNA <- unname(mutationTracker["mutSeq"])
  # seqList <- c(seqList, mutDNA)












  # #same mutations at each primer
  #
  # #Beginning single
  # # print("Beginning Single")
  # for(i in 1:(nrow(primerLocs) - 1)){
  #   mV <- mutateAtPrimer(i, mutDNA, primerLocs, primerDF, mutLoc = "S")
  #   mutDNA <- unname(mV["mSeq"])
  #   mutationTracker  <- rbind(mutationTracker, c(seqNum, i, "V", mV["MM"], mV["mPos"], mV["Orig"]))
  # }
  # seqList = c(seqList, mutDNA)
  #
  #
  # # next seq
  # seqNum = seqNum + 1
  # mutDNA = DNAs
  # #End single
  # # print("End Single")
  # for(i in 1:(nrow(primerLocs) - 1)){
  #   mV <- mutateAtPrimer(i, mutDNA, primerLocs, primerDF, mutLoc = "E")
  #   mutDNA <- unname(mV["mSeq"])
  #   mutationTracker  <- rbind(mutationTracker, c(seqNum, i, "V", mV["MM"], mV["mPos"], mV["Orig"]))
  #
  # }
  # seqList = c(seqList, mutDNA)
  #
  #
  # # next seq
  # seqNum = seqNum + 1
  # mutDNA = DNAs
  # #Middle
  # # print("Random Single")
  # for(i in 1:(nrow(primerLocs) - 1)){
  #   mV <- mutateAtPrimer(i, mutDNA, primerLocs, primerDF)
  #   mutDNA <- unname(mV["mSeq"])
  #   mutationTracker  <- rbind(mutationTracker, c(seqNum, i, "V", mV["MM"], mV["mPos"], mV["Orig"]))
  # }
  # seqList = c(seqList, mutDNA)


  # #Random Primer Mutation Loop
  # for (j in 1:25){
  #   seqNum = seqNum + 1
  #   mutDNA = DNAs
  #
  #   for(i in 1:(nrow(primerLocs) - 1)){
  #     isCospot = FALSE
  #
  #     #If primer 4
  #     if (i == 4){
  #       #Randomization to choose cospot 1/3 of the time
  #       n = sample(1:3, 1)
  #       if (n == 1){
  #         b = primerLocs[[s$LOC_START]][5]
  #         e = primerLocs[[s$LOC_END]][5]
  #         substr(mutDNA, b, e) <- primerDF[[s$SEARCHABLE_SEQ]][5]
  #         isCospot = TRUE
  #         i = 5
  #       }
  #     }
  #
  #     #Randomization to mutate 25% of the time
  #     n = sample(1:4, 1)
  #     if (n == 1){
  #       mV <- mutateAtPrimer(i, mutDNA, primerLocs, primerDF)
  #       mutDNA <- unname(mV["mSeq"])
  #
  #       if (isCospot){
  #         mutationTracker  <- rbind(mutationTracker,
  #                                   c(seqNum, i, "VC",  mV["MM"], mV["mPos"], mV["Orig"]))
  #       } else {
  #         mutationTracker  <- rbind(mutationTracker,
  #                                   c(seqNum, i, "V",  mV["MM"], mV["mPos"], mV["Orig"]))
  #       }
  #
  #     } else{
  #       if (isCospot){
  #         mutationTracker  <- rbind(mutationTracker,
  #                                   c(seqNum, i, "MC", NA, NA, NA))
  #       } else{
  #         mutationTracker  <- rbind(mutationTracker,
  #                                   c(seqNum, i, "M", NA, NA, NA))
  #       }
  #
  #     }
  #
  #   }
  #   seqList = c(seqList, mutDNA)
  # }


  # #Containing ambiguities and double mismatches
  # # print("Random Double Am")
  # for (j in 1:25){
  #   seqNum = seqNum + 1
  #   mutDNA = DNAs
  #
  #   for(i in 1:(nrow(primerLocs) - 1)){
  #
  #     isCospot = FALSE
  #     firstAmb = FALSE
  #     #if primer 4
  #     if (i == 4){
  #       #Randomization to choose cospot 1/3 of the time
  #       n = sample(1:3, 1)
  #       if (n == 1){
  #         b = primerLocs[[s$LOC_START]][5]
  #         e = primerLocs[[s$LOC_END]][5]
  #         substr(mutDNA, b, e) <- primerDF[[s$SEARCHABLE_SEQ]][5]
  #         isCospot = TRUE
  #         i = 5
  #       }
  #     }
  #
  #     #Randomization to mutate 25% of the time
  #     n = sample(1:4, 1)
  #     if (n == 1){
  #
  #       #Randomization for 2 mismatches 50% of the time
  #       n = sample(1:2, 1)
  #       if (n == 1){
  #         #Randomization for ambiguities 33% of the time
  #         n = sample(1:3, 1)
  #         if (n == 1){
  #           mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = TRUE)
  #           hasAmbiguity = TRUE
  #         } else{
  #           mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF)
  #           hasAmbiguity = FALSE
  #         }
  #         mutDNA <- unname(mV["mSeq"])
  #       } else {
  #         #TWO MISMATCHES
  #
  #         #Randomization for ambiguities 25% of the time
  #         n = sample(1:4, 1)
  #         if (n == 1){
  #           mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = TRUE)
  #           firstAmb = TRUE
  #           hasAmbiguity = TRUE
  #         } else{
  #           mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = FALSE)
  #           hasAmbiguity = FALSE
  #         }
  #         mutDNA <- unname(mV["mSeq"])
  #
  #         if(hasAmbiguity){
  #           varType <-"A"
  #         } else{
  #           if (isCospot){
  #             varType <-"VC"
  #           } else {
  #             varType <-"V"
  #           }
  #         }
  #         mutationTracker  <- rbind(mutationTracker,
  #                                   c(seqNum, i, varType,  mV["MM"], mV["mPos"], mV["Orig"]))
  #
  #         #Mutation 2
  #         #Randomization for ambiguities 25% of the time
  #         n = sample(1:4, 1)
  #         if (n == 1){
  #           mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = TRUE)
  #           hasAmbiguity = TRUE
  #         } else{
  #           mV <- extremeMutateAtPrimer(i, mutDNA, primerLocs, primerDF, ambiguity = FALSE)
  #           hasAmbiguity = FALSE
  #         }
  #         mutDNA <- unname(mV["mSeq"])
  #       }
  #
  #
  #
  #       if(hasAmbiguity | firstAmb){
  #         varType <-"A"
  #
  #       } else{
  #         if (isCospot){
  #           varType <-"VC"
  #         } else {
  #           varType <-"V"
  #         }
  #       }
  #       mutationTracker  <- rbind(mutationTracker,
  #                                 c(seqNum, i, varType,  mV["MM"], mV["mPos"], mV["Orig"]))
  #
  #     } else{
  #       if (isCospot){
  #         mutationTracker  <- rbind(mutationTracker,
  #                                   c(seqNum, i, "MC", NA, NA, NA))
  #       } else{
  #         mutationTracker  <- rbind(mutationTracker,
  #                                   c(seqNum, i, "M", NA, NA, NA))
  #       }
  #
  #     }
  #
  #   }
  #   seqList = c(seqList, mutDNA)
  # }

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
