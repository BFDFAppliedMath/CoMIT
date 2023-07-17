compKey2VarTypes <- function(DB, mutKey) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  varList <- dplyr::tbl(con, s$ID_MAP_TABLE) %>% dplyr::collect()
  varType <- dplyr::tbl(con, s$VAR_TYPE_TABLE) %>% dplyr::collect()
  ambTable <- dplyr::tbl(con, s$AMBIGUOUS_SEQS_TABLE) %>% dplyr::collect()

  DBI::dbDisconnect(con)

  # Check individual Mismatches--------------------------------------------------------
  # Iterate through VarList (Catches extra classifications in DB)
  for (row in 1:nrow(varList)) {
    a <- varList[[s$ACC_ID]][row]
    varId <- varList[[s$VAR_ID]][row]
    currPrimer <- varList[[s$PRIMER_ID]][row]

    varKey <- mutKey %>%
      dplyr::filter(get("Accession_ID") == a) %>%
      dplyr::filter(get("Primer_ID") == currPrimer)
    testVar <- varType %>% dplyr::filter(get(s$VAR_ID) == varId)

    # Equal number of mismatches for a variant
    if (nrow(testVar) != nrow(varKey)) {
      print(testVar)
      print(varKey)
    }
    testthat::expect_equal(nrow(testVar), nrow(varKey))


    for (i in 1:nrow(testVar)) {
      varMatch <- FALSE
      for (j in 1:nrow(varKey)) {
        # Check Mismatch Type and Position Match
        if (testVar[[s$VAR_TYPE]][i] == "Mismatch") {
          if ((varKey$Position[j] == testVar[[s$POS_FROM_SIDE]][i]) &
            (varKey[["Seq Change"]][j] == testVar[[s$SEQ_CHANGE]][i]) &
            (varKey[["Primer Template Mismatch"]][j] == testVar[[s$MISMATCH_TYPE]][i])) {
            varMatch <- TRUE
          }
        } else if (testVar[[s$VAR_TYPE]][i] == "Deletion") {
          if ((varKey$Position[j] == testVar[[s$POS_FROM_SIDE]][i]) &
            (varKey[["Seq Change"]][j] == testVar[[s$SEQ_CHANGE]][i])) {
            varMatch <- TRUE
          }
        } else if (testVar[[s$VAR_TYPE]][i] == "Insertion") {
          if ((varKey$Position[j] == testVar[[s$POS_FROM_SIDE]][i]) &
            (varKey[["Seq Change"]][j] == testVar[[s$SEQ_CHANGE]][i])) {
            varMatch <- TRUE
          }
        }
      }
      if (!varMatch) {
        print(a)
        # print("DB Rows")
        # print(varInfo %>% dplyr::filter(get(s$VAR_ID) == varId) %>%
        #         dplyr::pull(get(s$VAR_SEQ_COL)))
        # print("Key")
        # print(varKey %>% dplyr::pull(Mut_Seq))
        print("BD Rows")
        print(testVar)
        print("Key")
        print(varKey)
      }
      testthat::expect_true(varMatch)
    }
  }

  # Iterate through the Key and check for differences (Catches missed classification in DB)
  for (row in 1:nrow(mutKey)) {
    seq <- mutKey$Accession_ID[row]
    currPrimer <- mutKey$Primer_ID[row]
    isMut <- mutKey$M_V[row]

    if ((isMut == "M")) {
      next
    }

    if (isMut == "V") {
      varTotal <- varList %>%
        dplyr::filter(get(s$ACC_ID) == seq) %>%
        dplyr::filter(get(s$PRIMER_ID) == currPrimer)

      # Should only have one matching row
      testthat::expect_equal(nrow(varTotal), 1)

      varId <- varTotal[[s$VAR_ID]][1]

      varKey <- mutKey[row, ]
      testVar <- varType %>% dplyr::filter(get(s$VAR_ID) == varId)

      # Should be a variant type description
      testthat::expect_false(nrow(testVar) == 0)

      varMatch <- FALSE
      for (i in 1:nrow(testVar)) {
        if ((varKey$Position[1] == testVar[[s$POS_FROM_SIDE]][i]) &
          (varKey[["Seq Change"]][1] == testVar[[s$SEQ_CHANGE]][i]) &
          (varKey[["Primer Template Mismatch"]][1] == testVar[[s$MISMATCH_TYPE]][i])) {
          varMatch <- TRUE
        }
      }
      if (!varMatch) {
        print("Key")
        print(varKey)
        print("BD Rows")
        print(testVar)
      }
      testthat::expect_true(varMatch)
    }

    if ((isMut == "D") | (isMut == "I")) {
      varTotal <- varList %>%
        dplyr::filter(get(s$ACC_ID) == seq) %>%
        dplyr::filter(get(s$PRIMER_ID) == currPrimer)

      # Should only have one matching row
      testthat::expect_equal(nrow(varTotal), 1)

      varId <- varTotal[[s$VAR_ID]][1]

      varKey <- mutKey[row, ]
      testVar <- varType %>% dplyr::filter(get(s$VAR_ID) == varId)

      # Should be a variant type description
      testthat::expect_false(nrow(testVar) == 0)

      varMatch <- FALSE
      for (i in 1:nrow(testVar)) {
        if ((varKey$Position[1] == testVar[[s$POS_FROM_SIDE]][i]) &
          (varKey[["Seq Change"]][1] == testVar[[s$SEQ_CHANGE]][i])) {
          varMatch <- TRUE
        }
      }
      if (!varMatch) {
        print("Key")
        print(varKey)
        print("BD Rows")
        print(testVar)
      }
      testthat::expect_true(varMatch)
    }

    if (isMut == "A") {
      # DataSorter will default to most prevalent for unclassified
      if (currPrimer == 13) {
        currPrimer <- 12
      }
      ambPrimer <- ambTable %>%
        dplyr::filter(get(s$ACC_ID) == seq) %>%
        dplyr::filter(get(s$PRIMER_ID) == currPrimer)
      if (nrow(ambPrimer) != 1) {
        print(paste0(seq, " Primer: ", currPrimer))
        print(ambPrimer)
      } else {
        varKey <- mutKey %>%
          dplyr::filter(get("Accession_ID") == seq) %>%
          dplyr::filter(get("Primer_ID") == currPrimer)
        testthat::expect_equal(ambPrimer[[s$AMB_CHARS_COL]][1], varKey$`Seq Change`[1])
      }
      testthat::expect_equal(nrow(ambPrimer), 1)
    }
  }
}


compKey2SeqCounts <- function(DB, mutKey) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  seqInfo <- dplyr::tbl(con, s$SEQ_INFO_TABLE) %>% dplyr::collect()
  DBI::dbDisconnect(con)

  for (seqR in 1:nrow(seqInfo)) {
    seqID <- seqInfo[[s$ACC_ID]][seqR]
    # Get Key
    varCount <- mutKey %>%
      dplyr::filter(get("Accession_ID") == seqID) %>%
      dplyr::filter(get("M_V") == "V" | get("M_V") == "I" | get("M_V") == "D")

    varCountSingle <- length(unique(varCount$Primer_ID))

    # Get Var Count in DB
    vC <- seqInfo$Var_Count[seqR]

    if (vC != varCountSingle) {
      print(seqInfo[seqR, c("Accession_ID", "Var_Count")])
      print(varCount)
    }
    testthat::expect_equal(vC, varCountSingle)

    # Assays Affected
    aK <- length(unique(varCount$Assay_ID))

    aC <- seqInfo[[s$ASSAYS_AFFECTED]][seqR]

    if (aK != aC) {
      print(seqInfo[seqR, c("Accession_ID", "Assays_Affected")])
      print(varCount)
    }
    testthat::expect_equal(aC, aK)

    # Assays Affected WC (with Ambiguities)
    varCount <- mutKey %>%
      dplyr::filter(get("Accession_ID") == seqID) %>%
      dplyr::filter(get("M_V") != "M")

    aK <- length(unique(varCount$Assay_ID))

    aC <- seqInfo[[s$WITH_AMBS_AA]][seqR]

    if (aK != aC) {
      print(seqInfo[seqR, c("Accession_ID", "WC_Assays_Affected")])
      print(varCount)
    }
    testthat::expect_equal(aC, aK)

    # Assays Affected HR within 10bp
    varCount <- mutKey %>%
      dplyr::filter(get("Accession_ID") == seqID) %>%
      dplyr::filter(get("M_V") == "V" | get("M_V") == "I" | get("M_V") == "D") %>%
      dplyr::filter(Position <= 10)

    aK <- length(unique(varCount$Assay_ID))

    aC <- seqInfo[[s$BP10_AA]][seqR]

    if (aK != aC) {
      print(seqInfo[seqR, c("Accession_ID", "HR_Assays_Affected")])
      print(varCount)
    }
    testthat::expect_equal(aC, aK)
  }
}


checkIndelsInVarInfo <- function(DB, mutKey) {
  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  varInfo <- dplyr::tbl(con, s$VAR_INFO_TABLE) %>% dplyr::collect()
  primers <- dplyr::tbl(con, s$PRIMER_TABLE) %>% dplyr::collect()
  DBI::dbDisconnect(con)

  # Deltions
  indels <- mutKey %>% dplyr::filter(get("M_V") == "D" | get("M_V") == "I")

  for (i in 1:nrow(indels)) {
    seq <- indels$Mut_Seq[i]


    # Check Var Info
    varDel <- varInfo %>% dplyr::filter(get(s$PRIMER_ID) == indels$Primer_ID[i] &
      get(s$VAR_SEQ_COL) == seq)
    testthat::expect_equal(nrow(varDel), 1)

    indel_present <- FALSE
    if (varDel[[s$INDELS_PRESENT_VAR]] >= 1) {
      indel_present <- TRUE
    }
    testthat::expect_true(indel_present)
  }
}
