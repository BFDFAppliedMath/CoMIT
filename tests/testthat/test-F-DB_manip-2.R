testEnv <- new.env()

clearTestEnv <- function() {
  rm(list = ls(envir = testEnv), envir = testEnv)
}
clearTestEnv()

setStrings2AssayDB()

testEnv$DB <- build_COVID_DB("db_manip_test2", tempdir(),
                             include_all_assays = TRUE,
                             force = TRUE, silently = TRUE
)

withr::defer(file.remove(testEnv$DB))


# Tests---------------------------------------------------------

test_that("Success Assay from Primer", {
  # Pull from Database------------------------------
  assayKeyQuery <- paste0(
    "SELECT ", s$PRIMER_ID, ", ", s$ASSAY_ID,
    ' FROM "', s$PRIMER_TABLE, '"'
  )
  con <- DBI::dbConnect(RSQLite::SQLite(), testEnv$DB)
  testEnv$assayKey <- DBI::dbGetQuery(conn = con, assayKeyQuery)
  DBI::dbDisconnect(con)

  # Build library for cache variable
  getDBLibrary(testEnv$DB)

  expect_equal(getAssayFromPrimer(3), 1)
  expect_equal(getAssayFromPrimer(11), 3)
})

test_that("Failed Assay from Primer", {
  # Build library for cache variable
  getDBLibrary(testEnv$DB)

  expect_true(is.na(getAssayFromPrimer(40)))
})

test_that("Success Update Mismatch Count", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  varSeq <- "AAAAVARIANT"
  P_ID <- 6
  Assay_ID <- 2

  Var_ID <- addFoundVariant2DB(DB_FILE, varSeq, P_ID, Assay_ID)

  updateDBMismatchCount(5, Var_ID, DB_FILE)

  # Pull from Database------------------------------
  TestQuery <- paste0(
    "SELECT ", s$MISMATCH_NUM, " FROM '", s$VAR_INFO_TABLE,
    "' WHERE ", s$VAR_ID, " = '", Var_ID, "'"
  )
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)
  testDF <- DBI::dbGetQuery(conn = con, TestQuery)
  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Num_Mismatches[1], 5)
})

test_that("Failed Update Mismatch Count", {
  # Test Set up-----------------------------
  varSeq <- "AAAAVARIANT"
  P_ID <- 6
  Assay_ID <- 2

  Var_ID <- addFoundVariant2DB(testEnv$DB, varSeq, P_ID, Assay_ID)

  expect_error(updateDBMismatchCount(5, Var_ID, "not"))
})

test_that("Success Update Var and Assay Count", {
  # Test Set up-----------------------------
  DB_FILE <- testEnv$DB

  AccId <- "EPI_TEST"
  LOC <- "Outer Space"
  COL_D <- "2050-10"
  COL_DAY <- "31"
  name <- AccId
  PULL_D <- "2050-12-31"
  N_C <- .3
  CLADE <- "T-Rex"
  LINEAGE <- "Coolness"
  SUB_D <- "2050-12-25"

  insert2SeqInfo(DB_FILE, AccId, LOC, COL_D, COL_DAY, AccId, PULL_D, N_C, CLADE, LINEAGE)

  varCount <- 4
  assayCount <- 3
  WCcount <- 4
  HR_count <- 2


  updateDBvarAndAssayCount(AccId, varCount, assayCount, WCcount, HR_count, DB_FILE)


  # Pull from Database------------------------------
  TestQuery <- paste0("SELECT * FROM '", s$SEQ_INFO_TABLE, "' WHERE Accession_ID = '", AccId, "'")
  con <- DBI::dbConnect(RSQLite::SQLite(), DB_FILE)
  testDF <- DBI::dbGetQuery(conn = con, TestQuery)
  DBI::dbDisconnect(con)

  # Test Functions------------------------------------
  expect_equal(testDF$Assays_Affected[1], assayCount)
  expect_equal(testDF$Var_Count[1], varCount)
  expect_equal(testDF$HR_Assays_Affected[1], HR_count)
  expect_equal(testDF$WC_Assays_Affected[1], WCcount)
})

test_that("Failed Update Var and Assay Count", {
  # Test Set up-----------------------------
  AccId <- "EPI_TEST"

  varCount <- 4
  assayCount <- 3
  WCcount <- 4
  HR_count <- 2

  expect_error(updateDBvarAndAssayCount(AccId, varCount, assayCount, WCcount, HR_count, "not"))
})
