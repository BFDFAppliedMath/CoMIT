test_that("Database builds", {
  base <- system.file("sql", package = "CoMIT")
  schemaFile <- dir(base, "base_comit_schema.sql", f = TRUE)

  testName <- tempfile("test.db")
  expect_true(buildDBFromSchema(schemaFile, testName, silently = TRUE))
})

test_that("Update val builder", {
  expect_equal(
    genUpdate2Val("Var_Table", "MM", "ID"),
    "UPDATE \"Var_Table\" SET \"MM\"=? WHERE \"ID\"=?;"
  )
})

test_that("Insertion statement builder", {
  colsTest <- c("Cool_Col", "fun_level", "bounciness")
  expect_equal(
    genInsQuery("AwesomeTable", colsTest),
    "INSERT INTO 'AwesomeTable' ('Cool_Col', 'fun_level', 'bounciness') VALUES(?, ?, ?);"
  )
})



test_that("Select Statement Builder", {
  colsTest <- c("Cool_Col", "fun_level", "bounciness")
  whereTest1 <- c("ice", "cream")
  whereTest2 <- c("ice", "cream", "sunday")
  whereTest3 <- "bash"

  # NO WHERE STATEMENT
  expect_equal(
    genSelect("AwesomeTable", cols = colsTest),
    "SELECT \"Cool_Col\", \"fun_level\", \"bounciness\" FROM \"AwesomeTable\";"
  )

  # WHERE STATEMENT w/2
  expect_equal(
    genSelect("AwesomeTable", cols = colsTest, where = whereTest1),
    "SELECT \"Cool_Col\", \"fun_level\", \"bounciness\" FROM \"AwesomeTable\" WHERE \"ice\"=? AND \"cream\"=?;"
  )

  # WHERE STATEMENT w/3
  expect_equal(
    genSelect("AwesomeTable", cols = colsTest, where = whereTest2),
    "SELECT \"Cool_Col\", \"fun_level\", \"bounciness\" FROM \"AwesomeTable\" WHERE \"ice\"=? AND \"cream\" =? AND \"sunday\"=?;"
  )

  # WHERE STATEMENT w/3
  expect_equal(
    genSelect("AwesomeTable", cols = colsTest, where = whereTest3),
    "SELECT \"Cool_Col\", \"fun_level\", \"bounciness\" FROM \"AwesomeTable\" WHERE \"bash\"=?;"
  )
})

test_that("Aug by 1 statement builder", {
  expect_equal(
    genUpdateAugVal("Wanted_Posters", "Bounty_Price", "Sokka"),
    "UPDATE \"Wanted_Posters\" SET \"Bounty_Price\"=\"Bounty_Price\" + 1 WHERE \"Sokka\"=?;"
  )
})

test_that("Aug statement builder custom", {
  expect_equal(
    genUpdateAugValCustom("Wanted_Posters", "Bounty_Price", "Zuko"),
    "UPDATE \"Wanted_Posters\" SET \"Bounty_Price\"=\"Bounty_Price\" + ? WHERE \"Zuko\"=?;"
  )
})
