# Maybe build tests for individual function parts (fill, xl construct)

test_that("Check COVID Build", {
  setStrings2AssayDB()

  dbFile <- build_COVID_DB("test", tempdir(),
                           include_all_assays = TRUE,
                           force = FALSE, silently = TRUE
  )

  expect_true(file.exists(dbFile))


  TestQuery <- paste0("SELECT * FROM ", s$PRIMER_TABLE)

  con <- DBI::dbConnect(RSQLite::SQLite(), dbFile)

  testDF <- DBI::dbGetQuery(conn = con, TestQuery)

  DBI::dbDisconnect(con)

  expect_equal(testDF[[s$PRIMER_ID]][1], 1)
  expect_equal(testDF[[s$ASSAY_ID]][10], 3)
  expect_equal(testDF[[s$REACTION_COL]][2], "Inner")
  expect_equal(testDF[[s$DIRECTION_COL]][3], "Reverse")
  expect_equal(testDF[[s$P_SEQ_COL]][3], "ATTAGGTGAATTGTCCATAC")
  expect_equal(testDF[[s$S_SEQ_COL]][3], "GTATGGACAATTCACCTAAT")

  TestQuery <- paste0("SELECT * FROM ", s$PRIMER_LOCS_TABLE)

  con <- DBI::dbConnect(RSQLite::SQLite(), dbFile)

  testDF <- DBI::dbGetQuery(conn = con, TestQuery)

  DBI::dbDisconnect(con)


  expect_equal(testDF[[s$PRIMER_ID]][5], 5)
  expect_equal(testDF[[s$LOC_START]][5], 13988)
  expect_equal(testDF[[s$LOC_END]][5], 14011)

  expect_equal(testDF[[s$PRIMER_ID]][19], 19)
  expect_equal(testDF[[s$LOC_START]][19], 27999)
  expect_equal(testDF[[s$LOC_END]][19], 28024)

  file.remove(dbFile)
})
