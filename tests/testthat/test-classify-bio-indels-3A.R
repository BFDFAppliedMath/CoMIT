# Test with only 3A, uses same reference sequence


testEnv <- new.env()

clearTestEnv <- function() {
  rm(list = ls(envir = testEnv), envir = testEnv)
}
clearTestEnv()

setStrings2AssayDB()

testEnv$TESTDB <- build_COVID_DB("COVID_BIO_TEST_Indels", tempdir(),
                                 include_all_assays = FALSE,
                                 force = TRUE, silently = TRUE
)

testEnv$NUMBER_OF_PRIMERS <- 13

testEnv$mutationKey <- readr::read_csv("ref/classify-bio-indels-3A_KEY.csv")
# FOr troubleshootin tests
# testEnv$mutationKey <- readr::read_csv("tests/testthat/ref/classify-bio-indels-3A_KEY.csv")

testEnv$mutationKey <- testEnv$mutationKey %>% dplyr::rename_at("M/V", ~"M_V")

withr::defer(file.remove(testEnv$TESTDB))

# Tests---------------------------------------------------------
test_that("Successful Run", {
  covidTestSeqs <- "ref/indel_test.fasta"
  # covidTestSeqs <- "tests/testthat/ref/indel_test.fasta"

  myvars <- readr::cols_only(`Accession ID` = "c", Clade = "c", Lineage = "c")
  cladeDF <- readr::read_tsv("ref/indel_test_meta.tsv", col_types = myvars)
  # cladeDF <- readr::read_tsv("tests/testthat/ref/indel_test_meta.tsv", col_types = myvars)

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
