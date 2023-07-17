# TESTS for Biological COVID sequences ALL 7 ASSAYS
#  Data Sorter Biological Indel and Reversed Comp Seq Data Test

testEnv <- new.env()

rm(list = ls(envir = testEnv), envir = testEnv)


setStrings2AssayDB()

testEnv$TESTDB <- build_COVID_DB("COVID_BIO_TEST_Indels", tempdir(),
                                 include_all_assays = TRUE,
                                 force = TRUE, silently = TRUE
)


testEnv$testCaseDir <- "ref" # System run tests
# testEnv$testCaseDir = "tests/testthat/ref" #Troubleshooting

testEnv$NUMBER_OF_PRIMERS <- 30

testEnv$mutationKey <- readr::read_csv(paste0(
  testEnv$testCaseDir,
  "/classify_bio_indel_7A_KEY.csv"
))


testEnv$mutationKey <- testEnv$mutationKey %>% dplyr::rename_at("M/V", ~"M_V")

withr::defer(file.remove(testEnv$TESTDB))


# Tests---------------------------------------------------------
test_that("Successful Run", {
  covidTestSeqs <- paste0(testEnv$testCaseDir, "/indel_test.fasta")

  reversed_seq <- paste0(testEnv$testCaseDir, "/reversed_test.fasta")

  myvars <- readr::cols_only(`Accession ID` = "c", Clade = "c", Lineage = "c")
  cladeDF <- readr::read_tsv(paste0(testEnv$testCaseDir, "/indel_test_meta.tsv"),
    col_types = myvars
  )


  expect_error(comit_classify(covidTestSeqs, toString(Sys.Date()), testEnv$TESTDB, cladeDF), NA)

  expect_error(comit_classify(reversed_seq, toString(Sys.Date()), testEnv$TESTDB), NA)
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

test_that("Variant List and Variant Type Tables Check", {
  # Function compares all variations tracked in DB to key
  compKey2VarTypes(testEnv$TESTDB, testEnv$mutationKey)
})

test_that("Check Var Count and Assay Count on Seq_Info", {
  compKey2SeqCounts(testEnv$TESTDB, testEnv$mutationKey)
})

test_that("Variant Info Indels Check", {
  checkIndelsInVarInfo(testEnv$TESTDB, testEnv$mutationKey)
})
