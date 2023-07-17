# Synthetic Data with Indels and more mismatches

testEnv <- new.env()
rm(list = ls(envir = testEnv), envir = testEnv)

setStrings2AssayDB()

testEnv$testCaseDir <- "ref" # Sustem run tests
# testEnv$testCaseDir = "tests/testthat/ref" #Troubleshooting

TEST_DB_NAME <- "SYN_TEST_Indels"

primerCSV <- paste0(testEnv$testCaseDir, "/Primer_List_Syn_DNA.csv")

testEnv$synRef <- paste0(testEnv$testCaseDir, "/syn_ref.fasta")

testEnv$TESTDB <- build_Synthetic_DB(TEST_DB_NAME, tempdir(), primerCSV, testEnv$synRef)

testEnv$NUMBER_OF_PRIMERS <- 5

# Because of random order of classification order matters or sometimes you get a mismatch instead of a indel
set.seed(12)


testEnv$mutationKey <- readr::read_csv(paste0(
  testEnv$testCaseDir,
  "/classify_syn_indels_KEY.csv"
))
testEnv$mutationKey <- testEnv$mutationKey %>% dplyr::rename_at("M/V", ~"M_V")



# Tests---------------------------------------------------------

test_that("Successful Run", {
  synTestSeqs <- paste0(testEnv$testCaseDir, "/syn_indel.fasta")

  expect_error(comit_classify(synTestSeqs, toString(Sys.Date()), testEnv$TESTDB,
    refSeqFile = testEnv$synRef
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
