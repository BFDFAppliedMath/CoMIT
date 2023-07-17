testEnv <- new.env()

clearTestEnv <- function() {
  rm(list = ls(envir = testEnv), envir = testEnv)
}
clearTestEnv()

setStrings2AssayDB()

testEnv$DB <- build_COVID_DB("cospot_test", tempdir(),
                             include_all_assays = TRUE,
                             force = TRUE, silently = TRUE
)

getDBLibrary(testEnv$DB)

withr::defer(file.remove(testEnv$DB))


test_that("Closest Cospot Test", {
  # 16 is "CGTGGATGAGGCTGGTTC"
  # 17 is "CCTGGATGAGGCTGGTTC"
  closest <- getClosestCoSpot(16, "CGTGGATGAGGCTGGATG")
  expect_equal(closest$id, 16)
  expect_equal(closest$seq, "CGTGGATGAGGCTGGTTC")
  expect_equal(closest$indel, FALSE)

  closest <- getClosestCoSpot(16, "CCTGGATGAGGCTGGTTG")
  expect_equal(closest$id, 17)
  expect_equal(closest$seq, "CCTGGATGAGGCTGGTTC")
  expect_equal(closest$indel, FALSE)
})


test_that("Get Dominant Cospot Test", {
  TestQuery <- paste0("SELECT * FROM ", s$COSPOT_KEY)

  con <- DBI::dbConnect(RSQLite::SQLite(), testEnv$DB)

  coDF <- DBI::dbGetQuery(conn = con, TestQuery)

  DBI::dbDisconnect(con)

  doms <- getDominantCospotList(coDF)
  # print(doms)
  expect_equal(doms, c(16, 27))
})
