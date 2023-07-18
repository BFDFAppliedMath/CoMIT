###### Make Database Breakdown Table for COVID Inclusivity Reports ######

#' Format Database Breakdown Dataframe
#'
#' @param df the dataframe to be formatted
#' @param dateColName the user-inputted string describing the date analyzed, either
#' the supplied 3-month period or 1-month period, e.g. May - July 2022 or July 2022
#' @param total_seqs total number of sequences for the 3-month period
#' @param total_seqs_1m total number of sequences for the 1-month period
#' @param oneMonth a logical variable, whether the dataframe to be formatted describes data over a 3-month or 1-month period
#'
#' @return out_df, a dataframe with preliminary formatting applied (e.g. counts, percentages, column names) before it is combined with the other to create the gt table
#' @noRd
#'
formatDF <- function(df, dateColName, total_seqs, total_seqs_1m, oneMonth = FALSE) { #df_lastMonth, lastMonth, total_seqs, total_seqs_1m, oneMonth = TRUE
  Lineage <- Subgroup <- WHO_Designation <- LinSubgroup <- Displayed_Lineage <- Displayed_Group <- tempDate <- NULL

  dateColName <- paste0(dateColName, " # (%)")

  # Group and get counts
  lin_count <- df %>%
    dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group) %>%
    dplyr::summarise(tempDate = dplyr::n()) %>%
    dplyr::collect() %>%
    dplyr::ungroup()

  # Fix order
  lin_count <- lin_count %>% dplyr::select(c(Displayed_Group, Displayed_Lineage, WHO_Designation, Subgroup, tempDate))

  # Add "All" row for totals
  if (oneMonth) {
    lin_count <- rbind(lin_count, c(Displayed_Group = "All Sequences", Displayed_Lineage = "All", WHO_Designation = "All", Subgroup = "All", tempDate = total_seqs_1m))
  } else {
    lin_count <- rbind(lin_count, c(Displayed_Group = "All Sequences", Displayed_Lineage = "All", WHO_Designation = "All", Subgroup = "All", tempDate = total_seqs))
  }

  # change column name "tempDate" to the actual value of dateColName
  names(lin_count)[names(lin_count) == "tempDate"] <- dateColName

  out_df <- lin_count


  # Convert counts to numeric type
  out_df[, 5] <- sapply(out_df[, 5], as.numeric)

  # Add percentage column
  if (oneMonth) {
    out_df <- out_df %>%
      dplyr::mutate_at(dplyr::vars(names(out_df[5])), list("1M%" = ~ . / total_seqs_1m * 100))
  } else {
    out_df <- out_df %>%
      dplyr::mutate_at(dplyr::vars(names(out_df[5])), list("3M%" = ~ . / total_seqs * 100))
  }

  out_df <- out_df %>% dplyr::mutate_if(is.numeric, ~ round(., 2))

  # Remove subgroup col
  out_df <- out_df %>% dplyr::select(-c(Subgroup))

  # Rename some columns for gt table
  out_df <- out_df %>% dplyr::rename("rowgrp_col" = Displayed_Group, "Pangolin Lineage" = Displayed_Lineage, "WHO Label" = WHO_Designation)

  return(out_df)
}



#' Make Database Breakdown Dataframe
#'
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param interval string describing the first period analyzed, which is a 3-month period, e.g. May - July 2022
#' @param lastMonth string describing the second period analyzed, which is a 1-month period, e.g. July 2022
#' @param total_seqs total number of sequences for the 3-month period
#' @param total_seqs_1m total number of sequences for the 1-month period
#'
#' @return combined_df, a dataframe with data for both time periods and all needed formatting applied for creation of the gt table
#' @noRd
#'
makeDFdatabaseBreakdown <- function(DB, interval, lastMonth, total_seqs, total_seqs_1m) {
  Lineage <- Subgroup <- WHO_Designation <- Displayed_Lineage <- Displayed_Group <- NULL

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Table with Needed Information
  seqInfo <- dplyr::tbl(con, "Seq_Info") %>% dplyr::collect()
  DBI::dbDisconnect(con)

  # 3-month version---------------------------------
  df <- seqInfo %>% dplyr::select(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group)
  df <- formatDF(df, interval, total_seqs, total_seqs_1m)

  # 1-month version---------------------------------
  df_lastMonth <- seqInfo %>% dplyr::select(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group, s$COLLECTION_DATE)

  formatted_date <- lubridate::parse_date_time(lastMonth, orders = "%B %Y")
  formatted_date <- toString(format(formatted_date, "%Y-%m"))

  df_lastMonth <- df_lastMonth %>% dplyr::filter(get(s$COLLECTION_DATE) == formatted_date)

  df_lastMonth <- formatDF(df_lastMonth, lastMonth, total_seqs, total_seqs_1m, oneMonth = TRUE)

  # Now combine the 3-month and last month dataframes
  combined_df <- dplyr::full_join(df, df_lastMonth, keep = FALSE)

  # Last-minute formatting of combined df for gt table-------------------

  # Change Other and All to "" for gt table formatting (they will be row group names instead)
  combined_df[combined_df == "Other"] <- ""
  combined_df[combined_df == "All"] <- " "

  combined_df[is.na(combined_df)] <- 0

  return(combined_df)
}



#' Make Database Breakdown Table
#'
#' @description Make a table showing sequence counts categorized by lineage, comparing a given 3-month period to the last month of the 3-month period
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param variantFile csv file containing "Displayed_Lineage_Order" column for ordering row names/row groups in the GT table
#' @param interval string describing the first period analyzed, which is a 3-month period, e.g. May - July 2022
#' @param lastMonth string describing the second period analyzed, which is a 1-month period, e.g. July 2022
#' @param colWidth desired column width; increase if text is wrapping (applies to all columns, for both printed and saved version)
#' @param extraNote string with extra source note to be displayed at the bottom of the table
#' @param saveTable logical value indicating whether to save the tables to the specified saveFolder
#' @param saveFolder file path containing folder the finished tables will be saved in
#' @param height desired height of table .png images (when saveTable is TRUE)
#' @param width desired width of table .png images (when saveTable is TRUE)
#'
#' @export
#'
makedbBreakdownTable <- function(DB, variantFile, interval, lastMonth, colWidth, extraNote = NULL, saveTable = FALSE, saveFolder = NULL, height = NULL, width = NULL) {
  setStrings2AssayDB()

  `Pangolin Lineage` <- Lineage <- Subgroup <- WHO_Designation <- Displayed_Lineage <- Displayed_Group <- n <- `:=` <- NULL

  # get total seqs for 3 months and 1 month

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  seqInfo <- dplyr::tbl(con, "Seq_Info") %>% dplyr::collect()
  DBI::dbDisconnect(con)

  total_seqs <- seqInfo %>%
    dplyr::tally() %>%
    dplyr::pull(n)

  seqInfo <- seqInfo %>% dplyr::select(-c(Lineage))

  df_lastMonth <- seqInfo %>% dplyr::select(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group, s$COLLECTION_DATE)

  formatted_date <- lubridate::parse_date_time(lastMonth, orders = "%B %Y")
  formatted_date <- toString(format(formatted_date, "%Y-%m"))

  df_lastMonth <- df_lastMonth %>% dplyr::filter(get(s$COLLECTION_DATE) == formatted_date)

  total_seqs_1m <- df_lastMonth %>%
    dplyr::tally() %>%
    dplyr::pull(n)


  df <- makeDFdatabaseBreakdown(DB, interval, lastMonth, total_seqs, total_seqs_1m)

  col_names <- names(df)
  col_names <- col_names[-c(1, 5, 7)]

  # since pretty much everything is omicron, remove who label column
  df <- df %>% dplyr::select(-c(3))


  # Order rows based on Displayed Lineage Order column in lineage input file
  WHO_lin_des <- readr::read_csv(variantFile, show_col_types = FALSE)
  rowname_order <- WHO_lin_des$Displayed_Lineage_Order
  rowname_order <- rowname_order[!is.na(rowname_order)]
  rowname_order <- c(rowname_order, " ")

  df <- df %>%
    dplyr::mutate(`Pangolin Lineage` = forcats::fct_relevel(`Pangolin Lineage`, rowname_order)) %>%
    dplyr::arrange(`Pangolin Lineage`)


  # Create initial gt table with row groups, drop trailing zeros, and align center
  gt_tbl <- gt::gt(df, groupname_col = "rowgrp_col") %>%
    gt::fmt_number(
      columns = 3:6,
      drop_trailing_zeros = TRUE
    ) %>%
    gt::tab_style(
      style = gt::cell_text(align = "center"),
      locations = list(gt::cells_column_labels(), gt::cells_body())
    )

  # Add title
  gt_tbl <- gt_tbl %>%
    gt::tab_header(
      title = "Database Breakdown"
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bolder"),
      location = gt::cells_title()
    )

  # Format percent columns (the 3rd and 5th columns) and merge count and percent columns
  gt_tbl <- gt_tbl %>%
    gt::fmt_percent(
      scale_values = F,
      columns = c(4, 6),
      rows = gt::everything(),
      drop_trailing_zeros = TRUE
    ) %>%
    gt::cols_merge_n_pct(
      col_n = c(3),
      col_pct = c(4),
      autohide = TRUE
    ) %>%
    gt::cols_merge_n_pct(
      col_n = c(5),
      col_pct = c(6),
      autohide = TRUE
    )

  # Format count columns to be markdown format so that I can use line breaks
  gt_tbl <- gt_tbl %>% gt::cols_label(
    !!rlang::sym(col_names[3]) := gt::html(paste0(interval, "<br># (%)")),
    !!rlang::sym(col_names[4]) := gt::html(paste0(lastMonth, "<br># (%)"))
  )

  # Normalize column width
  gt_tbl <- gt_tbl %>% gt::cols_width(
    gt::everything() ~ px(colWidth)
  )

  # Add borders and boldness
  gt_tbl <- gt_tbl %>%
    gt::tab_style(
      style = gt::cell_borders(
        sides = c("left", "right"),
        color = "light grey"
      ),
      locations = list(gt::cells_body())
    ) %>%
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold"),
        gt::cell_borders(color = "light grey")
      ),
      locations = gt::cells_column_labels()
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = list(
        gt::cells_row_groups(groups = "All Sequences"),
        gt::cells_body(rows = `Pangolin Lineage` == " ")
      ) # the All row
    )

  # Add color formatting and note explaining the formatting
  gt_tbl <- gt_tbl %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#FFF2CC"), # light yellow
      locations = gt::cells_body(
        rows = abs(!!rlang::sym(col_names[3]) / total_seqs - !!rlang::sym(col_names[4]) / total_seqs_1m) >= 0.05
      )
    ) %>%
    gt::tab_source_note( # use paste0 to print >= with no space between it and the number
      source_note = paste0("Note: Yellow: ", "\U2265", "5% change between periods")
    ) %>%
    gt::tab_style(
      style = gt::cell_text(
        style = "italic",
        color = "#000000",
        size = "medium"
      ), locations = gt::cells_source_notes()
    )

  # Add extra source note to the bottom of the table if the user has specified one
  if (!is.null(extraNote)) {
    gt_tbl <- gt_tbl %>% gt::tab_source_note(
      source_note = extraNote
    )
  }


  # Add footnotes explaining recombinant lineages
  gt_tbl <- gt_tbl %>%
    gt::tab_footnote(
      footnote = "All Recombinant lineages of Omicron sublineages (e.g., XE is a recombinant of BA.1 and BA.2)",
      locations = gt::cells_body(
        columns = `Pangolin Lineage`,
        rows = `Pangolin Lineage` == "Omicron Recombinant"
      )
    ) %>%
    gt::tab_style(
      style = gt::cell_text(
        style = "italic",
        color = "#000000",
        size = "medium"
      ), locations = gt::cells_footnotes()
    )


  # Color unused cells grey
  if ("" %in% df$`Pangolin Lineage`) {
    gt_tbl <- gt_tbl %>% gt::tab_style(
      style = gt::cell_fill(color = "#D5D8DC"), # light grey
      locations = gt::cells_body(
        rows = `Pangolin Lineage` == "" | `Pangolin Lineage` == " ",
        columns = c("Pangolin Lineage")
      )
    )
  } else {
    gt_tbl <- gt_tbl %>% gt::tab_style(
      style = gt::cell_fill(color = "#D5D8DC"), # light grey
      locations = gt::cells_body(
        rows = `Pangolin Lineage` == " ",
        columns = c("Pangolin Lineage")
      )
    )
  }

  print(gt_tbl)

  if (saveTable) {
    if (!dir.exists(saveFolder)) {
      dir.create(saveFolder, recursive = TRUE)
    }

    gt_tbl %>% gt::gtsave(paste0(saveFolder, "database_breakdown.png"), vheight = height)
  }
}
