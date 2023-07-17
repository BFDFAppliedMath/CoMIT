setGISAIDstrings()

test_that("Metadata extraction", {
  metaWithSub <- tibble::tibble(!!s$GISAID_CLADE := c("G"),
    !!s$PANGO_LINEAGE := c("B.666"),
    "Submission date" = c("2020-03")
  )

  metaVector <- getMetaData(metaWithSub)
  expect_equal(metaVector[["clade"]], "G")
  expect_equal(metaVector[["lineage"]], "B.666")




  metaWithSub <- tibble::tibble(!!s$GISAID_CLADE := c("L"),
    !!s$PANGO_LINEAGE := c("B.20"),
    "Submission date" = c("2021-05-13")
  )

  metaVector <- getMetaData(metaWithSub)
  expect_equal(metaVector[["clade"]], "L")
  expect_equal(metaVector[["lineage"]], "B.20")


  # metaWithSub = tibble(!!GISAID_CLADE := c("G", "L", "O", "w"),
  #                      !!PANGO_LINEAGE := c("B.666", "B.20", "A", "M"),
  #                      "Submission date" = c("2020-03", "2021-05-13", "2070-3-9", "2111-9"))
})


test_that("FASTA Line extraction Full Date", {
  # Full Date
  name1 <- "hCoV-19/Iran/T11560/2021|EPI_ISL_2499970|2021-02-04"

  metaInfo <- parseFastaLine(name1)

  expect_equal(metaInfo$location, "Iran")
  expect_equal(metaInfo$accNum, "EPI_ISL_2499970")
  expect_equal(metaInfo$colDate, "2021-02")
  expect_equal(metaInfo$colDay, "04")
  expect_equal(metaInfo$continue, TRUE)
})

test_that("FASTA Line extraction No day", {
  # Only Month Year
  name1 <- "hCoV-19/Iran/T11560/2021|EPI_ISL_2499970|2021-02"

  metaInfo <- parseFastaLine(name1)

  expect_equal(metaInfo$location, "Iran")
  expect_equal(metaInfo$accNum, "EPI_ISL_2499970")
  expect_equal(metaInfo$colDate, "2021-02")
  expect_equal(metaInfo$colDay, NA)
  expect_equal(metaInfo$continue, TRUE)
})

test_that("FASTA Line extraction Single Digit Month", {
  # Single Digit Month
  name1 <- "hCoV-19/USA/NY-CUIMC-NP-3817/2020|EPI_ISL_2500603|2020-5-13"

  metaInfo <- parseFastaLine(name1)

  expect_equal(metaInfo$location, "USA")
  expect_equal(metaInfo$accNum, "EPI_ISL_2500603")
  expect_equal(metaInfo$colDate, "2020-05")
  expect_equal(metaInfo$colDay, "13")
  expect_equal(metaInfo$continue, TRUE)
})

test_that("FASTA Line extraction Single Digit Month and day", {
  # Single Digit Month and Day
  name1 <- "hCoV-19/USA/NY-CUIMC-NP-3817/2020|EPI_ISL_2500603|2020-5-2"

  metaInfo <- parseFastaLine(name1)

  expect_equal(metaInfo$location, "USA")
  expect_equal(metaInfo$accNum, "EPI_ISL_2500603")
  expect_equal(metaInfo$colDate, "2020-05")
  expect_equal(metaInfo$colDay, "02")
  expect_equal(metaInfo$continue, TRUE)
})


test_that("FASTA Line invalid input", {
  # Non-human
  name2 <- "hCoV-19/dog/Thailand/CU27042N/2021|EPI_ISL_2628963|2021-05-04"

  metaInfo <- parseFastaLine(name2)

  expect_equal(metaInfo$continue, FALSE)
  expect_equal(metaInfo$Fail_Type, "Host Error")

  # No EPI tag
  name3 <- "hCoV-19/Belgium/ULG-16349/2021|ISL_2499972|2021-06-04"

  metaInfo <- parseFastaLine(name3)

  expect_equal(metaInfo$continue, FALSE)
  expect_equal(metaInfo$Fail_Type, "NA Accession")

  # Just Year
  name4 <- "hCoV-19/Belgium/ULG-16349/2021|EPI_ISL_2499972|2021"

  metaInfo <- parseFastaLine(name4)

  expect_equal(metaInfo$continue, FALSE)
  expect_equal(metaInfo$Fail_Type, "Date is too short")
})

test_that("FASTA Line extraction EPI before accession", {
  name1 <- "hCoV-19/Wales/PHWC-PREPI4/2021|EPI_ISL_5290509|2021-10-02"

  metaInfo <- parseFastaLine(name1)

  expect_equal(metaInfo$location, "Wales")
  expect_equal(metaInfo$accNum, "EPI_ISL_5290509")
  expect_equal(metaInfo$colDate, "2021-10")
  expect_equal(metaInfo$colDay, "02")
  expect_equal(metaInfo$continue, TRUE)
})
