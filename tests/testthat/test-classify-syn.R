# TESTS for Biological COVID sequences ALL 7 ASSAYS
#  Data Sorter Biological Indel and Reversed Comp Seq Data Test

testEnv <- new.env()

rm(list = ls(envir = testEnv), envir = testEnv)

setStrings2AssayDB()

v <- create_test_fasta_and_DB()

testEnv$fastaFile <- v$FASTA
testEnv$TESTDB <- v$DatabaseFile
testEnv$mutationKey <- v$MutDF
testEnv$refSeqF <- v$refFile

withr::defer(file.remove(testEnv$TESTDB))

# Tests---------------------------------------------------------
test_that("Bad DB", {
  # File doesn't exist
  response <- comit_classify(testEnv$fastaFile, toString(Sys.Date()), "FAKE")
  expect_true(grepl("ERROR", response)[1])


  # File isn't a database
  response <- comit_classify(testEnv$fastaFile, toString(Sys.Date()), "tests/testthat.R")
  expect_true(grepl("ERROR", response)[1])

  # Incorrect tables
  response <- comit_classify(testEnv$fastaFile, toString(Sys.Date()), "tests/Test_DB/empty.db")
  expect_true(grepl("ERROR", response)[1])
})

test_that("Bad FASTA", {
  DB <- testEnv$TESTDB

  # File doesn't exist
  response <- comit_classify("FAKE", toString(Sys.Date()), DB)
  expect_true(grepl("ERROR", response)[1])


  # File isn't a fasta
  response <- comit_classify("tests/testthat.R", toString(Sys.Date()), DB)
  expect_true(grepl("ERROR", response)[1])
})

test_that("Run File", {
  # Function Throws no Errors (Following functions test for correct outputs)
  expect_error(comit_classify(testEnv$fastaFile,
    toString(Sys.Date()),
    testEnv$TESTDB,
    refSeqFile = testEnv$refSeqF
  ), NA)
})

test_that("Primer and Seq Table Counts", {
  DB <- testEnv$TESTDB
  mutKey <- testEnv$mutationKey

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  primer_db <- dplyr::tbl(con, s$PRIMER_TABLE) %>% dplyr::collect()
  seqInfo <- dplyr::tbl(con, s$SEQ_INFO_TABLE) %>% dplyr::collect()

  DBI::dbDisconnect(con)


  # Check that primers were added
  expect_equal(nrow(primer_db), 5)

  # Check number of Sequences
  totalseqs <- mutKey %>%
    dplyr::group_by(Seq_Num) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    nrow()
  expect_equal(nrow(seqInfo), totalseqs)
})

test_that("Variant Tables Counts", {
  DB <- testEnv$TESTDB
  mutKey <- testEnv$mutationKey

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  varInfo <- dplyr::tbl(con, s$VAR_INFO_TABLE) %>% dplyr::collect()
  varList <- dplyr::tbl(con, s$ID_MAP_TABLE) %>% dplyr::collect()
  DBI::dbDisconnect(con)


  # Same 2 Mismatch Numbers
  vars <- mutKey %>%
    dplyr::filter(M_V == "V" | M_V == "VC") %>%
    dplyr::group_by(Seq_Num, Primer_Num) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n == 2)

  doubleM <- varInfo %>% dplyr::filter(get(s$MISMATCH_NUM) == 2)

  expect_equal(nrow(doubleM), nrow(vars))

  # Check Total Number of Entries in the Var List Table
  vars <- mutKey %>% dplyr::filter(M_V == "V" | M_V == "VC")

  # Check if the varList is the correct length
  expect_equal(nrow(varList), nrow(vars) - nrow(doubleM))
})

test_that("Primer Match Check", {
  DB <- testEnv$TESTDB
  mutKey <- testEnv$mutationKey

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  varList <- dplyr::tbl(con, s$ID_MAP_TABLE) %>% dplyr::collect()
  DBI::dbDisconnect(con)


  # Invidivual Row check------------------------------------------------------------------------
  # Check Matches (if in the variation list then not classified as match)
  matches <- mutKey %>%
    dplyr::filter(M_V == "M" | M_V == "MC") %>%
    dplyr::select(Seq_Num, Primer_Num)
  for (row in 1:nrow(matches)) {
    varCount <- varList %>%
      dplyr::filter(get(s$ACC_ID) == matches$Seq_Num[row]) %>%
      dplyr::filter(get(s$PRIMER_ID) == matches$Primer_Num[row]) %>%
      dplyr::summarise(n = dplyr::n())
    expect_equal(varCount$n[1], 0)
  }
})

# Accessions are different, can't use function, maybe revise?
test_that("Variant Classification", {
  # Check Variant List and Variant Type Tables match Key
  DB <- testEnv$TESTDB
  mutKey <- testEnv$mutationKey

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  varList <- dplyr::tbl(con, s$ID_MAP_TABLE) %>% dplyr::collect()
  varType <- dplyr::tbl(con, s$VAR_TYPE_TABLE) %>% dplyr::collect()
  ambTable <- dplyr::tbl(con, s$AMBIGUOUS_SEQS_TABLE) %>% dplyr::collect()


  DBI::dbDisconnect(con)

  testCount <- 0
  # Check individual Mismatches--------------------------------------------------------
  # Iterate through VarList (Catches extra classifications in DB)
  for (row in 1:nrow(varList)) {
    a <- varList[[s$ACC_ID]][row]
    varId <- varList[[s$VAR_ID]][row]
    currPrimer <- varList[[s$PRIMER_ID]][row]
    aNum <- stringr::str_extract(a, "\\d+")

    varKey <- mutKey %>%
      dplyr::filter(Seq_Num == aNum) %>%
      dplyr::filter(Primer_Num == currPrimer)
    testVar <- varType %>% dplyr::filter(get(s$VAR_ID) == varId)

    # Equal number of mismatches for a variant
    expect_equal(nrow(testVar), nrow(varKey))
    if (nrow(testVar) != nrow(varKey)) {
      print(paste("Test:", nrow(testVar), "Key:", nrow(varKey)))
      print(testVar)
      print(varKey)
    }
    testCount <- testCount + 1

    for (i in 1:nrow(testVar)) {
      varMatch <- FALSE
      for (j in 1:nrow(varKey)) {
        keySeqChange <- paste0(varKey$Original[j], "->", varKey$Mismatch[j])

        # print("error here")
        # print(varKey)$Position[j]
        # print(testVar[[s$POS_FROM_SIDE]][i])
        # print(keySeqChange)
        # print(testVar[[s$SEQ_CHANGE]][i])

        if ((varKey$Position[j] == testVar[[s$POS_FROM_SIDE]][i]) &
          (keySeqChange == testVar[[s$SEQ_CHANGE]][i])) {
          varMatch <- TRUE
        }
      }
      if (!varMatch) {
        print("Key")
        print(varKey)
        print("BD Rows")
        print(testVar)
      }
      expect_true(varMatch)
      testCount <- testCount + 1
    }
  }

  # Iterate through the Key and check for differences (Catches missed classification in DB)
  for (row in 1:nrow(mutKey)) {
    seq <- mutKey$Seq_Num[row]
    currPrimer <- mutKey$Primer_Num[row]
    isMut <- mutKey$M_V[row]

    if ((isMut == "M") | (isMut == "MC")) {
      next
    }

    if ((isMut == "V") | (isMut == "VC")) {
      varTotal <- varList %>%
        dplyr::filter(get(s$ACC_ID) == paste0("EPI_ISL_", seq)) %>%
        dplyr::filter(get(s$PRIMER_ID) == currPrimer)

      # Should only have one matching row
      expect_equal(nrow(varTotal), 1)
      testCount <- testCount + 1

      varId <- varTotal[[s$VAR_ID]][1]

      varKey <- mutKey[row, ]
      testVar <- varType %>% dplyr::filter(get(s$VAR_ID) == varId)

      expect_false(nrow(testVar) == 0)
      testCount <- testCount + 1

      varMatch <- FALSE
      for (i in 1:nrow(testVar)) {
        keySeqChange <- paste0(varKey$Original[1], "->", varKey$Mismatch[1])
        if ((varKey$Position[1] == testVar[[s$POS_FROM_SIDE]][i]) &
          (keySeqChange == testVar[[s$SEQ_CHANGE]][i])) {
          varMatch <- TRUE
        }
      }
      if (!varMatch) {
        print("Key")
        print(varKey)
        print("BD Rows")
        print(testVar)
      }
      expect_true(varMatch)
      testCount <- testCount + 1
    }

    if (isMut == "A") {
      ambPrimer <- ambTable %>%
        dplyr::filter(get(s$ACC_ID) == paste0("EPI_ISL_", seq)) %>%
        dplyr::filter(get(s$PRIMER_ID) == currPrimer)

      # print(paste("acc", seq, "primer", currPrimer))

      containsAmbiguity <- FALSE

      for (aRow in 1:nrow(ambPrimer)) {
        if (!(ambPrimer[[s$AMB_CHARS_COL]][aRow] %in% c("A", "T", "G", "C"))) {
          containsAmbiguity <- TRUE
        }
      }
      if (!containsAmbiguity) {
        print(ambPrimer)
      }
      expect_true(containsAmbiguity)
      testCount <- testCount + 1
    }
  }
  # print(testCount)
})


# Different Seq names so can't use function, Maybe fix this?
test_that("Accession Variant Counts", {
  # Check Var Count and Assay Count on Seq_Info
  DB <- testEnv$TESTDB
  mutKey <- testEnv$mutationKey

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  seqInfo <- dplyr::tbl(con, s$SEQ_INFO_TABLE) %>% dplyr::collect()
  DBI::dbDisconnect(con)


  for (seqR in 1:nrow(seqInfo)) {
    seqID <- seqInfo[[s$ACC_ID]][seqR]
    # Get Key
    varCount <- mutKey %>%
      dplyr::filter(Seq_Num == seqR) %>%
      dplyr::filter(M_V == "V" | M_V == "VC")
    doubleCount <- varCount %>%
      dplyr::group_by(Seq_Num, Primer_Num) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    doubleCount <- doubleCount %>% dplyr::filter(n > 1)

    # Get Var Count in DB
    vC <- seqInfo[[s$SEQ_VAR_COUNT]][seqR]
    expect_equal(vC, nrow(varCount) - nrow(doubleCount))

    # Assays
    # Assays
    aK <- 1
    if (vC == 0) {
      aK <- 0
    }

    aC <- seqInfo[[s$ASSAYS_AFFECTED]][seqR]

    if (aK != aC) {
      print(seqInfo[seqR, ])
      print(varCount)
    }
    expect_equal(aC, aK)
  }
})
