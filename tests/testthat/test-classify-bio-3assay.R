# TESTS for Biological COVID sequences 3 ASSAYS ONLY

testEnv <- new.env()

clearTestEnv <- function() {
  rm(list = ls(envir = testEnv), envir = testEnv)
}
clearTestEnv()

setStrings2AssayDB()

testEnv$TESTDB <- build_COVID_DB("COVID_BIO_TEST", tempdir(),
                                 include_all_assays = FALSE,
                                 force = TRUE, silently = TRUE
)

testEnv$testCaseDir <- "ref" # System run tests
# testEnv$testCaseDir = "tests/testthat/ref" #Troubleshooting


testEnv$NUMBER_OF_PRIMERS <- 13

testEnv$mutationKey <- readr::read_csv(paste0(testEnv$testCaseDir, "/classify_bio_3A_KEY.csv"))

testEnv$mutationKey <- testEnv$mutationKey %>% dplyr::rename_at("M/V", ~"M_V")

withr::defer(file.remove(testEnv$TESTDB))




# Tests---------------------------------------------------------
test_that("Run File", {
  covidTestSeqs <- paste0(testEnv$testCaseDir, "/COVID_Bio_test_seqs.fasta")

  myvars <- readr::cols_only(`Accession ID` = "c", Clade = "c", Lineage = "c")
  cladeDF <- readr::read_tsv(paste0(testEnv$testCaseDir, "/COVID_Bio_test_meta.tsv"), col_types = myvars)


  expect_error(comit_classify(covidTestSeqs, toString(Sys.Date()), testEnv$TESTDB, cladeDF), NA)
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
  expect_equal(nrow(primer_db), testEnv$NUMBER_OF_PRIMERS)

  # Check number of Sequences
  totalseqs <- mutKey %>%
    dplyr::group_by(Accession_ID) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    nrow()
  expect_equal(nrow(seqInfo), totalseqs)
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
    dplyr::filter(M_V == "M") %>%
    dplyr::select(s$ACC_ID, s$PRIMER_ID)
  for (row in 1:nrow(matches)) {
    varCount <- varList %>%
      dplyr::filter(Accession_ID == matches$Accession_ID[row]) %>%
      dplyr::filter(Primer_ID == matches$Primer_ID[row]) %>%
      dplyr::summarise(n = dplyr::n())
    expect_equal(varCount$n[1], 0)
  }
})

test_that("Variant List and Variant Type Tables Check", {
  # Function compares all variations tracked in DB to key
  compKey2VarTypes(testEnv$TESTDB, testEnv$mutationKey)
})

test_that("Check Var Count and Assay Count on Seq_Info", {
  compKey2SeqCounts(testEnv$TESTDB, testEnv$mutationKey)
})

# Builds Var Sequence with Deletion, can't use function
test_that("Variant Info Check Indels", {
  DB <- testEnv$TESTDB
  mutKey <- testEnv$mutationKey

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  varInfo <- dplyr::tbl(con, s$VAR_INFO_TABLE) %>% dplyr::collect()
  primers <- dplyr::tbl(con, s$PRIMER_TABLE) %>% dplyr::collect()
  DBI::dbDisconnect(con)

  # Deltions
  dels <- mutKey %>% dplyr::filter(M_V == "D")

  for (i in 1:nrow(dels)) {
    ends <- as.numeric(unlist(stringr::str_extract_all(dels$`Seq Change`[i], "\\d+")))
    delPrimer <- primers %>%
      dplyr::filter(Primer_ID == dels$Primer_ID[i]) %>%
      dplyr::pull(Searchable_Seq)
    substr(delPrimer, ends[1], ends[2]) <- paste(rep("-", (ends[2] - ends[1] + 1)), collapse = "")

    # Check Var Info
    varDel <- varInfo %>% dplyr::filter(get(s$PRIMER_ID) == dels$Primer_ID[i] &
      get(s$VAR_SEQ_COL) == delPrimer)

    expect_equal(nrow(varDel), 1)

    expect_equal(varDel[[s$INDELS_PRESENT_VAR]], 1)
  }
})

test_that("Check GISAID NAMES/Seq info", {
  DB <- testEnv$TESTDB
  mutKey <- testEnv$mutationKey

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  seqInfo <- dplyr::tbl(con, s$SEQ_INFO_TABLE) %>% dplyr::collect()
  DBI::dbDisconnect(con)

  high_seqs_files <- "ref/gisaid_non_human.fasta"
  comit_classify(high_seqs_files, toString(Sys.Date()), DB)

  nonHumans <- c("EPI_ISL_2663240", "EPI_ISL_2521769", "EPI_ISL_2803214")

  seqInfoCount <- seqInfo %>% dplyr::filter(Accession_ID %in% nonHumans)

  expect_equal(nrow(seqInfoCount), 0)
})

test_that("High N Seqs", {
  DB <- testEnv$TESTDB
  mutKey <- testEnv$mutationKey

  high_seqs_files <- paste0(testEnv$testCaseDir, "/gisaid_high_n.fasta")
  comit_classify(high_seqs_files, toString(Sys.Date()), DB)


  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  high_n <- dplyr::tbl(con, s$HIGH_N_TABLE) %>% dplyr::collect()
  seqInfo <- dplyr::tbl(con, s$SEQ_INFO_TABLE) %>% dplyr::collect()

  DBI::dbDisconnect(con)

  expect_equal(nrow(high_n), 4)

  highNs <- c("EPI_ISL_416737", "EPI_ISL_413885", "EPI_ISL_413900", "EPI_ISL_416735")

  seqInfoCount <- seqInfo %>% dplyr::filter(get(s$ACC_ID) %in% highNs)

  expect_equal(nrow(seqInfoCount), 0)
})
