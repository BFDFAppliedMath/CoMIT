testEnv <- new.env()

clearTestEnv <- function() {
  rm(list = ls(envir = testEnv), envir = testEnv)
}
clearTestEnv()

setStrings2AssayDB()

testEnv$DB <- build_COVID_DB("db_manip_test", tempdir(),
                             include_all_assays = TRUE,
                             force = TRUE, silently = TRUE
)

withr::defer(file.remove(testEnv$DB))


# Tests---------------------------------------------------------
test_that("Successful Seq Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  AccId <- "EPI_TEST"
  LOC <- "Outer Space"
  COL_D <- "2050-10"
  DAY <- "31"
  FASTA_H <- AccId
  PULL_D <- "2050-12-31"
  P_N <- .3
  CLADE <- "T-Rex"
  LINEAGE <- "Coolness"

  insert2SeqInfo(DB_FILE, AccId, LOC, COL_D, DAY, FASTA_H, PULL_D, P_N, CLADE, LINEAGE)

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM ", s$SEQ_INFO_TABLE, " WHERE Accession_ID = '", AccId, "'")

  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)

  testDF <- DBI::dbGetQuery(conn = con, TestQuery)

  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Accession_ID[1], AccId)
  expect_equal(testDF$Location[1], LOC)
  expect_equal(testDF$Collection_Date[1], COL_D)
  expect_equal(testDF$Collection_Day[1], DAY)
  expect_equal(testDF$Lineage[1], LINEAGE)
  expect_equal(testDF$Clade[1], CLADE)
  expect_equal(testDF$FASTA_Header[1], FASTA_H)
  expect_equal(testDF$Pull_Date[1], PULL_D)
  expect_equal(testDF$Percent_N[1], P_N)
  expect_equal(testDF$Has_Ambiguities[1], 0)
  expect_equal(testDF$Has_Indel[1], 0)
})

test_that("Successful Seq Submission No Meta", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  AccId <- "EPI_TEST_2"
  LOC <- "Outer Space"
  COL_D <- "2050-10"
  DAY <- "31"
  FASTA_H <- AccId
  PULL_D <- "2050-12-31"
  P_N <- .3

  insert2SeqInfoNoMeta(DB_FILE, AccId, LOC, COL_D, DAY, FASTA_H, PULL_D, P_N)

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM ", s$SEQ_INFO_TABLE, " WHERE Accession_ID = '", AccId, "'")

  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)

  testDF <- DBI::dbGetQuery(conn = con, TestQuery)

  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Accession_ID[1], AccId)
  expect_equal(testDF$Location[1], LOC)
  expect_equal(testDF$Collection_Date[1], COL_D)
  expect_equal(testDF$Collection_Day[1], DAY)
  expect_equal(testDF$FASTA_Header[1], FASTA_H)
  expect_equal(testDF$Pull_Date[1], PULL_D)
  expect_equal(testDF$Percent_N[1], P_N)
  expect_equal(testDF$Has_Ambiguities[1], 0)
  expect_equal(testDF$Has_Indel[1], 0)
})

test_that("Failed Seq Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  AccId <- "EPI_TEST_Error"
  LOC <- "Tatooine"
  COL_D <- "2050-10"
  DAY <- "31"
  FASTA_H <- AccId
  PULL_D <- "2050-12-31"
  P_N <- .9999
  CLADE <- "Ewok"
  LINEAGE <- "Smile"

  # Wrong DB
  expect_error(insert2SeqInfo(tempfile(pattern = "seq_info"), AccId, LOC, COL_D, DAY, FASTA_H, PULL_D, P_N, CLADE, LINEAGE))

  # Duplicate Accession (Inserted prev test)
  expect_error(insert2SeqInfo(DB_FILE, "EPI_TEST", LOC, COL_D, DAY, FASTA_H, PULL_D, P_N, CLADE, LINEAGE))

  # String for an int (Should pass, add checks in the future?)
  insert2SeqInfo(DB_FILE, AccId, LOC, COL_D, DAY, FASTA_H, PULL_D, "One", CLADE, LINEAGE)

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM '", s$SEQ_INFO_TABLE, "' WHERE Accession_ID = '", AccId, "'")

  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)

  testDF <- DBI::dbGetQuery(conn = con, TestQuery)

  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(nrow(testDF), 1)
})


test_that("Successful High N Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  H_AccId <- "EPI_High_N"
  H_P_N <- 56.3

  insertToHighSeqsDB(DB_FILE, H_AccId, H_P_N)

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM '", s$HIGH_N_TABLE, "' WHERE Accession_ID = '", H_AccId, "'")
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)
  testDF <- DBI::dbGetQuery(conn = con, TestQuery)
  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Accession_ID[1], H_AccId)
  expect_equal(testDF$Seq_Percent_N[1], H_P_N)
})

test_that("Failed High N Submission", {
  # Test Set up-----------------------------
  DB_FILE <- "tests/COVID_TEST.db"

  H_AccId <- "EPI_High_N"
  H_P_N <- 56.3

  expect_error(insertToHighSeqsDB(tempfile(pattern = "high_n"), H_AccId, H_P_N))
})



test_that("Successful Variant Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  varSeq <- "AAAAVARIANT"
  P_ID <- 6
  Assay_ID <- 2

  Var_ID <- addFoundVariant2DB(DB_FILE, varSeq, P_ID, Assay_ID)

  expect_false(is.na(Var_ID))

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM '", s$VAR_INFO_TABLE, "' WHERE Var_ID = '", Var_ID, "'")
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)
  testDF <- DBI::dbGetQuery(conn = con, TestQuery)
  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Var_Seq[1], varSeq)
  expect_equal(testDF$Primer_ID[1], P_ID)
  expect_equal(testDF$Assay_ID[1], Assay_ID)
})

test_that("Failed Variant Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  varSeq <- "AAAAVARIANT"
  P_ID <- 6
  Assay_ID <- 2

  # Invalid DB
  expect_error(addFoundVariant2DB(tempfile(pattern = "add_variant"), varSeq, P_ID, Assay_ID))
})

test_that("Successful Variant AccessionID Map Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  AccId <- "EPI_Darth_Vader"
  varID <- 351
  P_ID <- 5
  Assay_ID <- 9

  add2IdMapTable(AccId, varID, P_ID, Assay_ID, DB_FILE)

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM '", s$ID_MAP_TABLE, "' WHERE Var_ID = '", varID, "' AND Accession_ID = '", AccId, "'")
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)
  testDF <- DBI::dbGetQuery(conn = con, TestQuery)
  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Accession_ID[1], AccId)
  expect_equal(testDF$Var_ID[1], varID)
  expect_equal(testDF$Primer_ID[1], P_ID)
  expect_equal(testDF$Assay_ID[1], Assay_ID)
})

test_that("Successful Variant AccessionID Map Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  AccId <- "EPI_Darth_Vader"
  varID <- 351
  P_ID <- 5
  Assay_ID <- 9

  expect_error(add2IdMapTable(AccId, varID, P_ID, Assay_ID))
  expect_error(add2IdMapTable(AccId, varID, P_ID, Assay_ID, tempfile(pattern = "var_accId_map")))
})

test_that("Successful Variant Type Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  varID <- 351
  class <- "Mismatch"
  changeString <- "A->G"
  pos <- 13
  P_ID <- 5
  mismatchType <- "G/T"

  insertToVariantTypeTableDB(DB_FILE, varID, class, changeString, pos, P_ID, mismatchType)

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM '", s$VAR_TYPE_TABLE, "' WHERE Var_ID = '", varID, "'")
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)
  testDF <- DBI::dbGetQuery(conn = con, TestQuery)
  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  # expect_equal(testDF$Var_ID[1], varID)
  expect_equal(testDF[[s$VAR_TYPE]][1], class)
  expect_equal(testDF$Seq_Change[1], changeString)
  expect_equal(testDF$Mismatch_Type[1], mismatchType)
  expect_equal(testDF$Pos_from_3P[1], pos)
  expect_equal(testDF$Primer_ID[1], P_ID)
})

test_that("Failed Variant Type Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  varID <- 351
  class <- "Mismatch"
  changeString <- "A->G"
  pos <- 13
  P_ID <- 5
  mismatchType <- "G/T"

  expect_error(insertToVariantTypeTableDB(tempfile(pattern = "variant_type"), varID, class, changeString, pos, P_ID, mismatchType))
})

test_that("Successful Ambiguity Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  AccId <- "EPI_Bad_Read"
  P_ID <- 9
  primerString <- "AASeagullsStopItNoW"
  ambiguityChars <- "W"

  addAmbiguity2DB(DB_FILE, AccId, P_ID, primerString, ambiguityChars)

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM '", s$AMBIGUOUS_SEQS_TABLE, "' WHERE Accession_ID = '", AccId, "'")
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)
  testDF <- DBI::dbGetQuery(conn = con, TestQuery)
  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Accession_ID[1], AccId)
  expect_equal(testDF$Var_Seq[1], primerString)
  expect_equal(testDF$Ambiguity_Chars[1], ambiguityChars)
  expect_equal(testDF$Seq_Percent_N[1], 5)
  expect_equal(testDF$Primer_ID[1], P_ID)


  AccId <- "EPI_Bad_Read_2"
  P_ID <- 49
  primerString <- "BushesOfLove"
  ambiguityChars <- c("B", "L")

  addAmbiguity2DB(DB_FILE, AccId, P_ID, primerString, ambiguityChars)

  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM '", s$AMBIGUOUS_SEQS_TABLE, "' WHERE Accession_ID = '", AccId, "'")
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)
  testDF <- DBI::dbGetQuery(conn = con, TestQuery)
  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Accession_ID[1], AccId)
  expect_equal(testDF$Var_Seq[1], primerString)
  expect_equal(testDF$Ambiguity_Chars[1], "B, L")
  expect_equal(testDF$Seq_Percent_N[1], 0)
  expect_equal(testDF$Primer_ID[1], P_ID)
})

test_that("Failed Ambiguity Submission", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  AccId <- "EPI_Bad_Read"
  P_ID <- 9
  primerString <- "AASeagullsStopItNoW"
  ambiguityChars <- "W"

  expect_error(addAmbiguity2DB(tempfile(pattern = "ambiguity"), AccId, P_ID, primerString, ambiguityChars))
})
