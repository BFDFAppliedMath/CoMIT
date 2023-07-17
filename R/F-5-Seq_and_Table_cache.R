# INPUT ERROR HANDLING FUNCTIONS
getDBLibrary <- function(DATABASE_FILE) {
  # QUERY GENERATION
  primerDFQuery <- genSelect(s$PRIMER_TABLE)
  cospotMapQuery <- genSelect(s$COSPOT_KEY)
  primerLocQuery <- genSelect(s$PRIMER_LOCS_TABLE)
  primerVarsQuery <- genSelect(s$VAR_INFO_TABLE,
    cols = c(
      s$VAR_ID,
      s$VAR_SEQ_COL,
      s$PRIMER_ID
    )
  )
  assayKeyQuery <- genSelect(s$PRIMER_TABLE,
    cols = c(
      s$PRIMER_ID,
      s$ASSAY_ID
    )
  )
  assayInfoQuery <- genSelect(s$ASSAY_TABLE)
  varCountQuery <- paste0(
    'SELECT "', s$VAR_ID,
    '", COUNT(*) AS Var_Count FROM "',
    s$ID_MAP_TABLE,
    '" GROUP BY "', s$VAR_ID, '"'
  )

  result <- tryCatch(
    expr = {
      # Pull from Database------------------------------
      if (file.exists(DATABASE_FILE)) {
        con <- DBI::dbConnect(RSQLite::SQLite(), DATABASE_FILE)


        cache$primerDF <- DBI::dbGetQuery(conn = con, primerDFQuery)

        cache$primerLocs <- DBI::dbGetQuery(conn = con, primerLocQuery)

        cache$primerVars <- DBI::dbGetQuery(conn = con, primerVarsQuery)

        cache$cospotMap <- DBI::dbGetQuery(conn = con, cospotMapQuery)

        cache$assayKey <- DBI::dbGetQuery(conn = con, assayKeyQuery) # from primer to assay

        cache$assayInfo <- DBI::dbGetQuery(conn = con, assayInfoQuery)

        VarCounts <- DBI::dbGetQuery(conn = con, varCountQuery)

        DBI::dbDisconnect(con)

        # Variant_Seq_Info Dataframe (Ordered with most common variation first)
        if (nrow(cache$primerVars) != 0) {
          cache$primerVars <- cache$primerVars %>%
            dplyr::filter(!grepl("-", s$VAR_SEQ_COL)) %>%
            dplyr::left_join(VarCounts, by = s$VAR_ID) %>%
            dplyr::arrange(dplyr::desc("Var_Count"))
          cache$primerVars <- cache$primerVars %>% dplyr::select(-"Var_Count")
        }

        cache$cospotPrimers <- dplyr::filter(
          cache$primerDF,
          get(s$IS_COSPOT_STRING) == 1
        )
        cache$primerDF <- dplyr::filter(
          cache$primerDF,
          get(s$IS_COSPOT_STRING) == 0
        )
      } else {
        fileBuilder <- c(
          "*****ERROR*******",
          paste(DATABASE_FILE, "does not exist"),
          "*****ERROR*******", " ", " "
        )
        # print(fileBuilder)
        return(list(error = fileBuilder))
      }
    },
    warning = function(w) {
      fileBuilder <- c(
        "*****ERROR*******",
        paste(w), "*****ERROR*******", " ", " "
      )
      print(fileBuilder)
      return(list(error = fileBuilder))
    },
    error = function(e) {
      fileBuilder <- c(
        "*****ERROR*******",
        paste(e), "*****ERROR*******", " ", " "
      )
      print(fileBuilder)
      return(list(error = fileBuilder))
    }
  )
  if ("error" %in% names(result)) {
    return(result)
  } else {
    return(NA)
  }
}

getSeqs <- function(fileName) {
  # Sequences
  result <- tryCatch(
    expr = {
      cache$dnaSeqs <- seqinr::read.fasta(fileName,
        as.string = TRUE,
        forceDNAtolower = FALSE,
        whole.header = TRUE
      )
    },
    warning = function(w) {
      if (grepl("final line", w)) {
        return(NA)
      } else {
        fileBuilder <- c(
          "*****ERROR*******",
          paste(fileName, "could not be read"),
          paste("Check FASTA File List"),
          "*****ERROR*******", " ", " "
        )
        print(w)
        return(list(error = fileBuilder))
      }
    },
    error = function(e) {
      fileBuilder <- c(
        "*****ERROR*******",
        paste(fileName, "could not be read"),
        paste("Check FASTA File List"),
        "*****ERROR*******", " ", " "
      )
      print(e)
      return(list(error = fileBuilder))
    }
  )

  if ("error" %in% names(result)) {
    return(result)
  } else {
    return(NA)
  }
}
