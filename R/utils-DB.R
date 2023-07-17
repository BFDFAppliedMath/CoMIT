#' Build a database from a given SQL schema file
#'
#' @description This is a generic database build function that can create a database file from any .sql schema. If you use a schema other than the default one, later functions in the pipeline may not work.
#' @param schemaFile a .sql schema file (if no value given, defaults to the base Covid-19 schema used in build_COVID_DB())
#' @param dbName name of database (with file path)
#' @param force if true, replace DB if already exists; if false, will not overwrite file
#' @param silently whether or not to print tables made
#'
#' @export
#'
buildDBFromSchema <- function(schemaFile = NULL, dbName, force = FALSE, silently = FALSE) {
  dbSuccessful <- FALSE

  if(is.null(schemaFile)) {
    base <- system.file("sql", package = "CoMIT")
    schemaFile <- dir(base, "base_comit_schema.sql", full.names = TRUE)
  }

  if (!file.exists(schemaFile)) {
    warning("SQL schema file does not exist. Please check your file.")
    return()
  }

  if (file.exists(dbName)) {
    if (!force) {
      warning("File exists in database write location. To overwrite file, run with force = TRUE")
      return()
    }
  }

  result <- tryCatch(
    expr = {
      queries <- readLines(schemaFile)

      queries <- stringr::str_replace_all(queries, "--.*$", " ") ### remove any commented lines
      queries <- paste(queries, collapse = "\n") ### collapse with new lines
      schemaUpdates <- unlist(stringr::str_split(queries, "(?<=;)"))

      con <- DBI::dbConnect(RSQLite::SQLite(), dbname = dbName)

      for (s in schemaUpdates) {
        if (nchar(s) > 1 & grepl("\\S", s)) { # check if valid sql statement
          DBI::dbExecute(con, s)
        }
      }

      if (!silently) {
        print("The following tables were built:")
        print(DBI::dbListTables(con))
      }

      DBI::dbDisconnect(con)

      dbSuccessful <- TRUE
    }, error = function(e) {
      print("Connection or SQL type error")
      print(e)
    }
  )
  return(dbSuccessful)
}

runSQLstring <- function(sqlstrings, dbCon, silently = FALSE) {
  dbSuccessful <- FALSE

  result <- tryCatch(
    expr = {
      queries <- stringr::str_replace_all(sqlstrings, "--.*$", " ") ### remove any commented lines
      queries <- paste(queries, collapse = "\n") ### collapse with new lines
      schemaUpdates <- unlist(stringr::str_split(queries, "(?<=;)"))

      for (s in schemaUpdates) {
        if (nchar(s) > 1) {
          rs <- DBI::dbExecute(dbCon, s)

          if (!silently) {
            print(rs)
          }
        }
      }

      dbSuccessful <- TRUE
    }, error = function(e) {
      print("Connection or SQL type error")
      print(e)
    }
  )


  return(dbSuccessful)
}

genInsQuery <- function(table_name, insVars) {
  V <- paste0(insVars, collapse = "', '")
  q <- paste0(rep("?", length(insVars)), collapse = ", ")

  return(paste0("INSERT INTO '", table_name, "' ('", V, "') VALUES(", q, ");"))
}

genUpdate2Val <- function(table, updateCol, whereCol) {
  statement <- paste0('UPDATE "', table, '" SET "', updateCol, '"=? WHERE "', whereCol, '"=?;')
  return(statement)
}

genUpdateAugVal <- function(table, updateCol, whereCol) {
  statement <- paste0('UPDATE "', table, '" SET "', updateCol, '"="', updateCol, '" + 1 WHERE "', whereCol, '"=?;')
  return(statement)
}

genUpdateAugValCustom <- function(table, updateCol, whereCol) {
  statement <- paste0('UPDATE "', table, '" SET "', updateCol, '"="', updateCol, '" + ? WHERE "', whereCol, '"=?;')
  return(statement)
}

genSelect <- function(table, cols = NULL, where = NULL) {
  if (!is.null(cols)) {
    V <- paste0(cols, collapse = '", "')
    selState <- paste0('SELECT "', V, '" FROM "', table, '"')
  } else {
    selState <- paste0('SELECT * FROM "', table, '"')
  }
  if (is.null(where)) {
    where <- ";"
  } else if (length(where) > 1) {
    W <- paste0(where[2:length(where)], collapse = '" =? AND "')
    where <- paste0(' WHERE "', where[1], '"=? AND "', W, '"=?;')
  } else {
    where <- paste0(' WHERE "', where, '"=?;')
  }
  return(paste0(selState, where))
}
