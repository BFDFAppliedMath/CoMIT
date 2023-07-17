insertIntoSeq <- function(seq, pos, insertion) {
  newSeq <- gsub(
    paste0("^(.{", pos, "})(.*)$"),
    paste0("\\1", insertion, "\\2"),
    seq
  )

  return(newSeq)
}



# for changes to combo of locations, make indels starting on right and update lengths based on previous adds/deletes

# mutDNA, primerLocs, primerDF, assayID=2, indel=NA, indelRegionNums=NA, changeGC="dec"


# indel can be NA, "ins", or "del"
# indelRegionNums can be NA or any combo of c(1,2,3); 1 is between OF and IF, 2 is between IF and IR, 3 is between IR and OR
# changeGC can be NA, "inc" or "dec" (always in inner amplicon)
mutateAtAmplicon <- function(ogDNA, primerLocs, primerDF, assayID, indel = NA, indelRegionNums = NA, changeGC = NA) { # assayID=1, indel="ins", indelRegionNums=2, changeGC=NA

  # 1. make df with OFstart, OFend, IFstart, IFend, IRstart, IRend, ORstart, ORend for this assay
  assayDF <- primerDF %>%
    dplyr::filter(get(s$ASSAY_ID) == assayID) %>%
    dplyr::left_join(primerLocs, by = "Primer_ID") %>%
    dplyr::select(s$PRIMER_ID, s$REACTION_COL, s$DIRECTION_COL, s$LOC_START, s$LOC_END)

  OF_start <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Outer", get(s$DIRECTION_COL) == "Forward") %>%
    dplyr::filter(dplyr::row_number() == 1) %>% # first cospot is dominant cospot in Covid
    dplyr::pull(get(s$LOC_START))

  IF_start <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Inner", get(s$DIRECTION_COL) == "Forward") %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(get(s$LOC_START))

  IF_end <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Inner", get(s$DIRECTION_COL) == "Forward") %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(get(s$LOC_END))

  IR_start <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Inner", get(s$DIRECTION_COL) == "Reverse") %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(get(s$LOC_START))

  IR_end <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Inner", get(s$DIRECTION_COL) == "Reverse") %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(get(s$LOC_END))

  OR_end <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Outer", get(s$DIRECTION_COL) == "Reverse") %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(get(s$LOC_END))

  ogInnerLength <- IR_end - IF_start + 1
  og_inner_seq <- substr(ogDNA, IF_start, IR_end)

  ogOuterLength <- OR_end - OF_start + 1
  og_outer_seq <- substr(ogDNA, OF_start, OR_end)



  innerStart <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Inner" & get(s$DIRECTION_COL) == "Forward") %>%
    dplyr::pull(get(s$LOC_END)) + 1
  innerEnd <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Inner" & get(s$DIRECTION_COL) == "Reverse") %>%
    dplyr::pull(get(s$LOC_START)) - 1
  outerStart <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Outer" & get(s$DIRECTION_COL) == "Forward") %>%
    dplyr::pull(get(s$LOC_END)) + 1
  outerEnd <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Outer" & get(s$DIRECTION_COL) == "Reverse" & get(s$PRIMER_ID) != 5) %>%
    dplyr::pull(get(s$LOC_START)) - 1 # First is dominant cospot

  # ogInnerLengthIndel <- innerEnd - innerStart + 1
  # og_inner_seq_indel <- substr(ogDNA, innerStart, innerEnd)

  # ogOuterLengthIndel <- outerEnd - outerStart + 1
  # og_outer_seq_indel <- substr(ogDNA, outerStart, outerEnd)

  # for outer amplicon changes:
  region1End <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Inner" & get(s$DIRECTION_COL) == "Forward") %>%
    dplyr::pull(get(s$LOC_START)) - 1 # bp before start of IF
  region3Start <- assayDF %>%
    dplyr::filter(get(s$REACTION_COL) == "Inner" & get(s$DIRECTION_COL) == "Reverse") %>%
    dplyr::pull(get(s$LOC_END)) + 1 # bp after end of IR



  # 2. Return these things:
  # mutationTracker["mutSeq"]  entire mutated sequence
  # mutationInfo["ogInnerLength"],
  # mutationInfo["mutInnerLength"],
  # mutationInfo["ogOuterLength"],
  # mutationInfo["mutOuterLength"],
  # mutationInfo["ogGC"],
  # mutationInfo["mutGC"]

  # indels ---------------------------------

  if (!is.na(indelRegionNums)) {
    if (indelRegionNums == 2) { # inner amplicon

      if (indel == "ins") {
        # ins_len <- cache$TOO_BIG_INNER_AMPLICON - ogInnerLengthIndel

        too_big_len <- 2 * ogInnerLength
        ins_len <- too_big_len - ogInnerLength

        ins_start <- sample(innerStart:innerEnd, size = 1)
        mutInnerLength <- ogInnerLength + ins_len

        ins_seq <- sample(Biostrings::DNA_ALPHABET[1:4], size = ins_len, replace = TRUE)
        ins_seq <- paste(ins_seq, collapse = "")

        # update seq to have the insertion
        mutSeq <- ogDNA
        mutSeq <- insertIntoSeq(mutSeq, ins_start, ins_seq)

        # get mutated amplicon seq
        mutAmplicon <- substr(mutSeq, IF_start, IR_end + ins_len)

        return(c(
          "mutSeq" = mutSeq,
          "ogInnerLength" = ogInnerLength,
          "mutInnerLength" = mutInnerLength,
          "ogOuterLength" = NA,
          "mutOuterLength" = NA,
          "ogGC" = NA,
          "mutGC" = NA
        ))
      } else if (indel == "del") {
        # del_len <- ogInnerLength - cache$TOO_SMALL_INNER_AMPLICON

        too_small_len <- ceiling(0.8 * ogInnerLength)
        del_len <- ogInnerLength - too_small_len

        # get a deletion start position that wouldn't delete the Inner Reverse primer
        # del_start <- innerEnd + 1
        # while(del_start + del_len - 1 > innerEnd){
        #   del_start <- sample(innerStart:innerEnd, size=1)
        # }

        del_start <- sample(innerStart:(innerEnd - (del_len + 1)), size = 1)

        mutInnerLength <- ogInnerLength - del_len

        # update seq to have the deletion
        mutSeq <- ogDNA
        del_seq <- substr(mutSeq, del_start, del_start + del_len - 1)
        mutSeq <- stringr::str_remove(mutSeq, del_seq)

        # get mutated amplicon seq
        mutAmplicon <- substr(mutSeq, IF_start, IF_end - del_len)

        return(c(
          "mutSeq" = mutSeq,
          "ogInnerLength" = ogInnerLength,
          "mutInnerLength" = mutInnerLength,
          "ogOuterLength" = NA,
          "mutOuterLength" = NA,
          "ogGC" = NA,
          "mutGC" = NA
        ))
      }
    } else if (indelRegionNums == 1) { # indel in left part of outer amplicon
      if (indel == "ins") {
        # ins_len <- cache$TOO_BIG_OUTER_AMPLICON - ogOuterLength

        too_big_len <- 2 * ogOuterLength
        ins_len <- too_big_len - ogOuterLength

        ins_start <- sample(outerStart:region1End, size = 1)
        mutOuterLength <- ogOuterLength + ins_len

        ins_seq <- sample(Biostrings::DNA_ALPHABET[1:4], size = ins_len, replace = TRUE)
        ins_seq <- paste(ins_seq, collapse = "")

        # update seq to have the insertion
        mutSeq <- ogDNA
        mutSeq <- insertIntoSeq(mutSeq, ins_start, ins_seq)

        # get mutated amplicon seq
        mutAmplicon <- substr(mutSeq, OF_start, OR_end + ins_len)

        return(c(
          "mutSeq" = mutSeq,
          "ogInnerLength" = NA,
          "mutInnerLength" = NA,
          "ogOuterLength" = ogOuterLength,
          "mutOuterLength" = mutOuterLength,
          "ogGC" = NA,
          "mutGC" = NA
        ))
      } else if (indel == "del") {
        # del_len <- ogOuterLength - cache$TOO_SMALL_INNER_AMPLICON

        too_small_len <- ceiling(0.8 * ogOuterLength)
        del_len <- ogOuterLength - too_small_len

        # get a deletion start position that wouldn't delete the Inner Forward primer
        # del_start <- region1End + 1
        # while(del_start + del_len - 1 > region1End){
        #   del_start <- sample(outerStart:region1End, size=1)
        # }

        del_start <- sample(outerStart:(region1End - (del_len + 1)), size = 1)

        mutOuterLength <- ogOuterLength - del_len

        # update seq to have the deletion
        mutSeq <- ogDNA
        del_seq <- substr(mutSeq, del_start, del_start + del_len - 1)
        mutSeq <- stringr::str_remove(mutSeq, del_seq)

        # get mutated amplicon seq
        mutAmplicon <- substr(mutSeq, OF_start, OR_end - del_len)

        return(c(
          "mutSeq" = mutSeq,
          "ogInnerLength" = NA,
          "mutInnerLength" = NA,
          "ogOuterLength" = ogOuterLength,
          "mutOuterLength" = mutOuterLength,
          "ogGC" = NA,
          "mutGC" = NA
        ))
      }
    } else if (indelRegionNums == 3) { # indel in right part of outer amplicon
      if (indel == "ins") {
        # ins_len <- cache$TOO_BIG_OUTER_AMPLICON - ogOuterLength

        too_big_len <- 2 * ogOuterLength
        ins_len <- too_big_len - ogOuterLength

        ins_start <- sample(region3Start:outerEnd, size = 1)
        mutOuterLength <- ogOuterLength + ins_len

        ins_seq <- sample(Biostrings::DNA_ALPHABET[1:4], size = ins_len, replace = TRUE)
        ins_seq <- paste(ins_seq, collapse = "")

        # update seq to have the insertion
        mutSeq <- ogDNA
        mutSeq <- insertIntoSeq(mutSeq, ins_start, ins_seq)

        # get mutated amplicon seq
        mutAmplison <- substr(mutSeq, OF_start, OR_end + ins_len)

        return(c(
          "mutSeq" = mutSeq,
          "ogInnerLength" = NA,
          "mutInnerLength" = NA,
          "ogOuterLength" = ogOuterLength,
          "mutOuterLength" = mutOuterLength,
          "ogGC" = NA,
          "mutGC" = NA
        ))
      } else if (indel == "del") {
        # del_len <- ogOuterLength - cache$TOO_SMALL_INNER_AMPLICON

        too_small_len <- ceiling(0.8 * ogOuterLength)
        del_len <- ogOuterLength - too_small_len

        # get a deletion start position that wouldn't delete the Outer Reverse primer
        # del_start <- outerEnd + 1
        # while(del_start + del_len - 1 > outerEnd){
        #   del_start <- sample(region3Start:outerEnd, size=1)
        # }

        del_start <- sample(region3Start:(outerEnd - (del_len + 1)), size = 1)


        mutOuterLength <- ogOuterLength - del_len

        # update seq to have the deletion
        mutSeq <- ogDNA
        del_seq <- substr(mutSeq, del_start, del_start + del_len - 1)
        mutSeq <- stringr::str_remove(mutSeq, del_seq)

        # get mutated amplicon seq
        mutAmplicon <- substr(mutSeq, OF_start, OR_end - del_len)

        return(c(
          "mutSeq" = mutSeq,
          "ogInnerLength" = NA,
          "mutInnerLength" = NA,
          "ogOuterLength" = ogOuterLength,
          "mutOuterLength" = mutOuterLength,
          "ogGC" = NA,
          "mutGC" = NA
        ))
      }
    }
  } # end of indels


  # GC content change -------------------------------------

  if (!is.na(changeGC)) {
    num_g <- stringr::str_count(og_inner_seq, "G")
    num_c <- stringr::str_count(og_inner_seq, "C")
    total_og_num_GC <- num_g + num_c
    ogGC <- (total_og_num_GC) / ogInnerLength

    if (changeGC == "inc") {
      # mutGC <- cache$TOO_MUCH_GC

      mutGC <- ogGC + 0.2 # abs(new_GC - ref_GC) >= .5

      total_mut_num_GC <- mutGC * ogInnerLength
      num_GC_to_add <- ceiling(total_mut_num_GC - total_og_num_GC)

      gc_seq <- rep("G", times = num_GC_to_add)
      gc_seq <- paste(gc_seq, collapse = "")

      # replace the beginning of the inner amplicon sequence with G's

      # initialize mutated sequence to original sequence
      mutSeq <- ogDNA

      # sequence we will replace with G's
      og_seq <- substr(mutSeq, innerStart, innerStart + num_GC_to_add - 1)

      # entire mutated sequence
      mutSeq <- sub(og_seq, gc_seq, mutSeq)

      # mutated sequence at inner amplicon
      mutAmplicon <- substr(mutSeq, IF_start, IR_end)

      return(c(
        "mutSeq" = mutSeq,
        "ogInnerLength" = NA,
        "mutInnerLength" = NA,
        "ogOuterLength" = NA,
        "mutOuterLength" = NA,
        "ogGC" = ogGC,
        "mutGC" = mutGC
      ))
    } else if (changeGC == "dec") {
      # mutGC <- cache$TOO_LITTLE_GC

      mutGC <- ogGC - 0.2

      total_mut_num_GC <- mutGC * ogInnerLength
      num_GC_to_subtract <- ceiling(total_og_num_GC - total_mut_num_GC)

      mutSeq <- ogDNA
      og_inner_seq <- substr(ogDNA, innerStart, innerEnd)

      # replace number of G/C with A/T needed to get to too little GC content
      mutAmplicon <- og_inner_seq
      for (i in 1:num_GC_to_subtract) {
        random_replacement_nucleotide <- sample(c("A", "T"), size = 1)
        mutAmplicon <- sub("[GC]", random_replacement_nucleotide, mutAmplicon) # this is actually just from end of IF to beginning of IR, it's hard to get actual inner amplicon and we don't need it
      }

      # replace the inner amplicon seq with mutAmplicon
      mutSeq <- sub(og_inner_seq, mutAmplicon, mutSeq)

      return(c(
        "mutSeq" = mutSeq,
        "ogInnerLength" = NA,
        "mutInnerLength" = NA,
        "ogOuterLength" = NA,
        "mutOuterLength" = NA,
        "ogGC" = ogGC,
        "mutGC" = mutGC
      ))
    }
  }




















  # # other way to do mutations (based on region, not amplicon, so doensn't work with test-helper-generate_syn_fasta_amplicons)
  #
  # if(!is.na(indelRegionNums)){
  #
  #   # make dataframe with start and end positions of regions to add indels to
  #
  #   regionsdf <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("regionNum", "regionStart", "regionEnd"))))
  #
  #   for(region in seq_along(indelRegionNums)){
  #     if(region==1){
  #       primerIDs <- primerDF %>% dplyr::filter(get(s$ASSAY_ID)==assayID, (get(s$DIRECTION_COL)=="Forward" & get(s$REACTION_COL)=="Outer") | (get(s$DIRECTION_COL)=="Forward" & get(s$REACTION_COL)=="Inner")) %>%
  #         dplyr::pull(get(s$PRIMER_ID))
  #
  #     } else if(region==2){
  #       primerIDs <- primerDF %>% dplyr::filter(get(s$ASSAY_ID)==assayID, (get(s$DIRECTION_COL)=="Forward" & get(s$REACTION_COL)=="Inner") | (get(s$DIRECTION_COL)=="Reverse" & get(s$REACTION_COL)=="Inner")) %>%
  #         dplyr::pull(get(s$PRIMER_ID))
  #
  #     } else if(region==3){
  #       primerIDs <- primerDF %>% dplyr::filter(get(s$ASSAY_ID)==assayID, (get(s$DIRECTION_COL)=="Reverse" & get(s$REACTION_COL)=="Inner") | (get(s$DIRECTION_COL)=="Reverse" & get(s$REACTION_COL)=="Outer")) %>%
  #         dplyr::pull(get(s$PRIMER_ID))
  #     }
  #
  #     locs <- primerLocs %>% dplyr::filter(get(s$PRIMER_ID) %in% primerIDs)
  #     region_start <- locs[1, ] %>% dplyr::pull(get(s$LOC_END))
  #     region_end <- locs[2, ] %>% dplyr::pull(Start)
  #
  #     regionsdf <- regionsdf %>% rbind(data.frame(regionNum = region, regionStart = region_start, regionEnd = region_end))
  #   }
  #
  #   # add columns to df to be returned
  #   regionsdf <- regionsdf %>%
  #     dplyr::mutate(ogRegionLen = regionEnd - regionStart, mutRegionLen = NA, indelStart = NA, indelEnd = NA, indelSeq = NA, ogRegionSeq = NA, mutRegionSeq = NA, mutSeq = NA)
  #
  #   ogDNA <- mutDNA
  #
  #   for(i in nrow(regionsdf)){
  #
  #     if(indel=="ins"){
  #       # set insertion length to too big length (so that Tm is affected)
  #       ins_len = cache$TOO_BIG_INNER_AMPLICON - regionsdf[i, ]$ogRegionLen
  #
  #       # generate random insertion position in region (position is in reference to entire seq):
  #       ins_start <- sample(regionsdf[i, ]$regionStart:regionsdf[i, ]$regionEnd, size=1)
  #
  #       # add ins_start and ins_end to df
  #       regionsdf[i, ]$indelStart <- ins_start
  #       regionsdf[i, ]$indelEnd <- ins_start + ins_len
  #       regionsdf[i, ]$mutRegionLen <- regionsdf[i, ]$ogRegionLen + ins_len
  #
  #       # generate random sequence for insertion and add to df
  #       ins_seq <- sample(Biostrings::DNA_ALPHABET[1:4], size=ins_len, replace=TRUE)
  #       ins_seq <- paste(ins_seq, collapse="")
  #       regionsdf[i, ]$indelSeq <- ins_seq
  #
  #       # get ogRegionSeq and add to df
  #       ogRegionSeq <- substr(ogDNA, regionsdf[i, ]$regionStart, regionsdf[i, ]$regionEnd)
  #       regionsdf[i, ]$ogRegionSeq <- ogRegionSeq
  #
  #       # update region seq to have the insertion and add to df
  #       mutDNASeq <- ogDNA
  #       mutDNASeq <- insertIntoSeq(mutDNASeq, ins_start, ins_seq)
  #       regionsdf[i, ]$mutSeq <- mutDNASeq
  #       mutRegionSeq <- substr(mutDNASeq, regionsdf[i, ]$regionStart, regionsdf[i, ]$regionEnd + ins_len - 1) # idk why the -1 is needed
  #       regionsdf[i, ]$mutRegionSeq <- mutRegionSeq
  #     }
  #
  #     else if(indel=="del"){
  #       # set deletion length to too small length (so that Tm is affected)
  #       del_len = regionsdf[i, ]$ogRegionLen - cache$TOO_SMALL_INNER_AMPLICON
  #
  #       # etc
  #
  #     }
  #
  #
  #   }
  #
  #
  #
  # }
  #
  # # make row to be added to mutations df (1 row per mutated seq)
}

# mutateAtAmpliconWithCospots <- function() ????



# mutLoc can be NA, "S", or "E"
mutateAtPrimer <- function(primerNum, mutDNA, primerLocs, primerDF, mutLoc = NA) {
  b <- primerLocs[[s$LOC_START]][primerNum]
  e <- primerLocs[[s$LOC_END]][primerNum]
  lenP <- e - b + 1

  # Where to mutate
  mLoc <- sample(b:e, 1)

  if (!is.na(mutLoc) & mutLoc == "S") {
    mLoc <- b
  } else if (!is.na(mutLoc) & mutLoc == "E") {
    mLoc <- e
  }


  mutations <- c("A", "T", "G", "C")
  # get mutation
  o <- substr(mutDNA, mLoc, mLoc)
  m <- sample(mutations, 1) # m is the actual mutation nucleotide
  while (o == m) { # keep getting random sample until the mutation is different from the original nucleotide
    m <- sample(mutations, 1)
  }

  substr(mutDNA, mLoc, mLoc) <- m # update mutDNA to have the mutation

  pos <- mLoc - b + 1 # pos of mutation
  primer <- dplyr::filter(primerDF, get(s$PRIMER_ID) == primerNum)
  direction <- primer[[s$DIRECTION_COL]][1]
  if (direction == "Forward") {
    pos <- e - mLoc + 1
  }

  # MM is mutation sequence
  # mSeq is the whole sequence with the mutation
  # mPos is the start position of mutation
  # Orig is the original sequence at the mutation location
  return(c("MM" = m, "mSeq" = mutDNA, "mPos" = pos, "Orig" = o))
}


# createDeletionAtPrimer <- function(primerNum, mutDNA, primerLocs, primerDF, mutLoc = NA, size = 1){
#   b = primerLocs[[s$LOC_START]][primerNum]
#   e = primerLocs[[s$LOC_END]][primerNum]
#   lenP = e-s + 1
#
#   #Where to mutate
#   #Random location
#   mLoc = sample(s:e, 1)
#   deletionStart = mLoc
#   deletionEnd = mLoc + size
#   if(deletionEnd > e){
#     deletionEnd = e
#   }
#   if (!is.na(mutLoc) & mutLoc == "S") {
#     mLoc = b
#   } else if (!is.na(mutLoc) & mutLoc == "E") {
#     mLoc = e
#   }
#
#
#
#   if(size > lenP){
#     size <- lenP
#   }
#
#   substr(mutDNA, ends[1], ends[2]) <- paste(rep("-", (ends[2] - ends[1] + 1)), collapse = "")
#
#   #Get whole primer sequence with mutation for Var_Info Table
#
#   pos = mLoc - b + 1
#   primer = dplyr::filter(primerDF, get(s$PRIMER_ID) == primerNum)
#   direction = primer[[s$DIRECTION_COL]][1]
#   if(direction == "Forward"){
#     pos = e - mLoc + 1
#   }
#   d = paste0("Deletion from ", delStart, " to ", delEnd)
#
#   return(c("MM" = d, "mSeq" = mutDNA, "mPos" = pos, "Orig" = o))
# }


extremeMutateAtPrimer <- function(primerNum, mutDNA, primerLocs, primerDF, ambiguity = FALSE, prevMLoc = NA) {
  b <- primerLocs[[s$LOC_START]][primerNum]
  e <- primerLocs[[s$LOC_END]][primerNum]
  lenP <- e - b + 1

  # Set to true for first iteration
  sameMutationPos <- TRUE

  while (sameMutationPos) {
    mLoc <- sample(b:e, 1)

    pos <- mLoc - b + 1
    primer <- dplyr::filter(primerDF, get(s$PRIMER_ID) == primerNum)
    direction <- primer[[s$DIRECTION_COL]][1]
    if (direction == "Forward") {
      pos <- e - mLoc + 1
    }

    if (is.na(prevMLoc) | prevMLoc != pos) {
      sameMutationPos <- FALSE
    }
  }

  mutations <- c("A", "T", "G", "C")
  if (ambiguity == TRUE) {
    mutations <- c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N")
  }

  # get mutation
  o <- substr(mutDNA, mLoc, mLoc)
  m <- sample(mutations, 1)
  while (o == m) {
    m <- sample(mutations, 1)
  }

  substr(mutDNA, mLoc, mLoc) <- m

  pos <- mLoc - b + 1
  primer <- dplyr::filter(primerDF, get(s$PRIMER_ID) == primerNum)
  direction <- primer[[s$DIRECTION_COL]][1]
  if (direction == "Forward") {
    pos <- e - mLoc + 1
  }

  return(c("MM" = m, "mSeq" = mutDNA, "mPos" = pos, "Orig" = o))
}

insert_Primer_for_test <- function(assayNum, direction, place, seq, cospotStat, DB) {
  searchable <- seq
  # if(direction == "Reverse"){
  #   searchable <- c2s(rev(comp(s2c(searchable), ambiguous = TRUE, forceToLower = FALSE)))
  # }

  insCols <- c(
    s$ASSAY_ID,
    s$DIRECTION_COL, s$REACTION_COL, s$P_SEQ_COL,
    s$S_SEQ_COL, s$IS_COSPOT_STRING
  )

  primerUS <- genInsQuery(s$PRIMER_TABLE, insCols)

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  rs <- DBI::dbSendStatement(
    conn = con, primerUS,
    params = list(assayNum, direction, place, seq, searchable, cospotStat)
  )
  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}

insert_loc_for_test <- function(Primer_ID, Start, End, DB) {
  primerlocQ <- genInsQuery(s$PRIMER_LOCS_TABLE, c(s$PRIMER_ID, s$LOC_START, s$LOC_END))

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  rs <- DBI::dbSendStatement(
    conn = con, primerlocQ,
    params = list(Primer_ID, Start, End)
  )
  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}

insert_cospot_for_test <- function(c1, c2, DB) {
  cospotQ <- genInsQuery(s$COSPOT_KEY, c(s$SPOT1, s$SPOT2))

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  rs <- DBI::dbSendStatement(conn = con, cospotQ, params = list(c1, c2))
  DBI::dbClearResult(rs)
  rs <- DBI::dbSendStatement(conn = con, cospotQ, params = list(c2, c1))
  DBI::dbClearResult(rs)
  DBI::dbDisconnect(con)
}

buildPrimers <- function(DNAstring, DB, multipleAssays = FALSE) {
  assayNum <- 1
  primer_starts <- c(20, 60, 100, 150, 150)
  primer_ends <- c(40, 85, 117, 175, 175)
  dirs <- c("Forward", "Forward", "Reverse", "Reverse", "Reverse")
  reactions <- c("Outer", "Inner", "Inner", "Outer", "Outer")
  cospotStat <- c(rep(0, 3), rep(1, 2))

  for (i in seq_along(primer_starts)) {
    primerSeq <- substring(DNAstring, primer_starts[i], primer_ends[i])
    if (i == 5) {
      substr(primerSeq, 6, 6) <- "C"
    }
    insert_Primer_for_test(assayNum, dirs[i], reactions[i], primerSeq, cospotStat[i], DB)
    insert_loc_for_test(i, primer_starts[i], primer_ends[i], DB)
  }

  insert_cospot_for_test(4, 5, DB)

  if (isTRUE(multipleAssays)) {
    assayNum <- 2
    primer_starts <- c(220, 260, 300, 350)
    primer_ends <- c(240, 285, 317, 375)
    dirs <- c("Forward", "Forward", "Reverse", "Reverse")
    reactions <- c("Outer", "Inner", "Inner", "Outer")
    cospotStat <- c(rep(0, 4), rep(1, 2))

    for (i in seq_along(primer_starts)) {
      primerSeq <- substring(DNAstring, primer_starts[i], primer_ends[i])
      if (i == 5) {
        substr(primerSeq, 6, 6) <- "C"
      }
      insert_Primer_for_test(assayNum, dirs[i], reactions[i], primerSeq, cospotStat[i], DB)
      insert_loc_for_test(i + 5, primer_starts[i], primer_ends[i], DB)
    }
  }
}

generateSeqNames <- function(seqRange) {
  seqNames <- c()
  for (i in 1:seqRange) {
    locations <- c("USA", "Ireland", "Ghana", "China", "Russia", "Mexico", "Brazil", "Canada", "South Africa")

    ss <- paste0("hCoV-19/", sample(locations, 1)[1], "/coded-info/2021|EPI_ISL_", i, "|", sample(seq(as.Date("2020/01/01"), as.Date(Sys.Date()), by = "day"), 1)[1])

    seqNames <- c(seqNames, ss)
  }

  return(seqNames)
}
