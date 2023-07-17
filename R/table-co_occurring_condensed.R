###### Make Co-Occurring Mutations Tables for COVID Inclusivity Reports ######


# ---- Helper functions ----
light_blue_builder <- function(cols, condition) {
  gt::cells_body(
    columns = !!rlang::sym(cols),
    rows = !!rlang::sym(cols) / !!rlang::sym(condition) < 0.01
  )
}


light_yellow_builder <- function(cols, condition) {
  gt::cells_body(
    columns = !!rlang::sym(cols),
    rows = !!rlang::sym(cols) / !!rlang::sym(condition) > 0.05
  )
}


#' Format Co-occurring Mutations GT Table V2
#'
#' @param in_tbl GT table object
#' @param HR aka high risk, a boolean value indicating whether or not the input table has only data associated with high-risk mutations
#' @param colWidth1 width of "Pangolin Lineage" column in px
#' @param colWidth2 width of 0, 1, 2, and 3 columns in px
#' @param colWidth3 width of "5" column in px
#' @param colWidth4 width of ">=6" column in px
#' @param colWidth5 width of "Sequences by Lineage" column in px
#'
#' @return formatted GT table
#'
format_co_gt_tbl_2 <- function(in_tbl, colWidth1 = NULL, colWidth2 = NULL, colWidth3 = NULL, colWidth4 = NULL, colWidth5 = NULL, HR = FALSE) {
  `Pangolin Lineage` <- `Sequences by Lineage` <- `Total # Sequences` <- `7` <- NULL


  # Add headers and footnote
  out_tbl <- in_tbl %>% gt::tab_spanner(
    label = paste0("# Assays Affected"),
    columns = c("0", "1", "2", "3", "4", "5", !!s$NEW_COL_NAME)
  )

  if (HR) {
    out_tbl <- out_tbl %>% gt::tab_footnote(
      footnote = "Non-ambiguous mutations under the primer binding regions that fall
    within 10bp of 3' end of the primer were considered in this analysis.",
      locations = gt::cells_title()
    )
  }

  out_tbl <- out_tbl %>%
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

  # Add color formatting and note explaining the formatting
  cols <- c("0", "1", "2", "3", "4", "5", {{ s$NEW_COL_NAME }})
  last_col <- c("Sequences by Lineage")

  out_tbl <- out_tbl %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#CDE9EB"), # light blue
      locations = lapply(cols, light_blue_builder, condition = "Sequences by Lineage")
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#FFE699"), # light yellow
      locations = lapply(cols, light_yellow_builder, condition = "Sequences by Lineage")
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#D5D8DC"), # grey
      locations = gt::cells_body(columns = `Sequences by Lineage`)
    ) %>%
    gt::tab_source_note( # use paste0 to print >= with no space between it and the number
      source_note = paste0("Color guide - Yellow: mutations occur at ", "\U2265", "5%; Blue: mutations occur at <1%")
    ) %>%
    gt::tab_style(
      style = gt::cell_text(
        style = "italic",
        color = "#000000",
        size = "medium"
      ), locations = gt::cells_source_notes()
    )

  # Add text styling and borders
  out_tbl <- out_tbl %>%
    gt::tab_style(
      style = gt::cell_text(align = "center"),
      locations = list(gt::cells_column_labels(), gt::cells_body())
    ) %>%
    gt::fmt_number(
      columns = c(cols, last_col),
      use_seps = TRUE,
      sep_mark = ",",
      drop_trailing_zeros = TRUE
    ) %>%
    gt::tab_style(
      style = gt::cell_borders(
        sides = c("left", "right"),
        color = "light grey"
      ), locations = list(gt::cells_body())
    ) %>%
    gtExtras::gt_add_divider(
      columns = c(!!s$NEW_COL_NAME),
      sides = "right",
      color = "grey",
      weight = gt::px(3)
    )

  if (!is.null(colWidth1)) {
    out_tbl <- out_tbl %>% gt::cols_width(
      c(`Pangolin Lineage`) ~ gt::px(colWidth1)
    )
  }

  if (!is.null(colWidth2)) {
    out_tbl <- out_tbl %>% gt::cols_width(
      c(`0`, `1`, `2`, `3`) ~ gt::px(colWidth2)
    )
  }

  if (!is.null(colWidth3)) {
    out_tbl <- out_tbl %>% gt::cols_width(
      c(`4`, `5`) ~ gt::px(colWidth3)
    )
  }

  if (!is.null(colWidth4)) {
    out_tbl <- out_tbl %>% gt::cols_width(
      c(!!s$NEW_COL_NAME) ~ gt::px(colWidth4)
    )
  }

  if (!is.null(colWidth5)) {
    out_tbl <- out_tbl %>% gt::cols_width(
      c(`Sequences by Lineage`) ~ gt::px(colWidth5)
    )
  }

  out_tbl <- out_tbl %>% gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = list(
      gt::cells_row_groups(groups = "Summary: All Sequences by Assays"),
      gt::cells_body(columns = gt::everything(), rows = `Pangolin Lineage` == " ")
    ) # the All row
  )

  # change height of rows
  out_tbl <- out_tbl %>%
    gt::tab_options(data_row.padding = gt::px(1))


  return(out_tbl)
}




# ---- Main functions ----


#' Get Co-Occurring Counts V2
#'
#' @description Make a dataframe with counts of sequences containing co-occurring mutations, organized by WHO designations/Pango lineages
#'
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param HR if TRUE, returns 3 prime end, i.e. high risk, version
#' @param startDate start date for subset DB, if NULL then the .1 percent analysis won't be run
#' @param endDate end date for subset DB, if NULL then the .1 percent analysis won't be run
#' @param archiveDBFile File path to archive db (used to retrieve older subset of data for comparison; older subset uses same length of time as input subset db, e.g. start_date is "2022-5" and end_date is "2022-07", then older subset start date is "2022-02 and "2022-04)
#'
#' @return df, a dataframe with mutated sequence counts co-occurring on 0 to 7 assays
#'
getCoOccCounts_2 <- function(DB, HR = FALSE, startDate = NULL, endDate = NULL, archiveDBFile = NULL) {
  `:=` <- `.` <- `0` <- `5` <- `6` <- `7` <- new_col <- Assays_Affected <- count <- All <- HR_Assays_Affected <- Subgroup <- Lineage <- Displayed_Lineage <- Displayed_Group <- WHO_Designation <- Collection_Date <- prop_change <- NULL

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  seqInfo <- dplyr::tbl(con, "Seq_Info") %>% dplyr::collect()
  assayInfo <- dplyr::tbl(con, "Assay_Info") %>%
    dplyr::select(s$ASSAY_ID, s$ASSAY_NAME) %>%
    dplyr::collect()

  DBI::dbDisconnect(con)

  total_seqs <- seqInfo %>%
    dplyr::summarise(total_seqs = dplyr::n()) %>%
    dplyr::pull(total_seqs)

  seqInfo <- seqInfo %>% dplyr::select(-c(Lineage))

  lin_df <- seqInfo %>%
    dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group) %>%
    dplyr::summarise(`Total # Sequences` = dplyr::n()) %>%
    dplyr::ungroup()


  if (isTRUE(HR)) { # mutations on the 3 prime end
    total <- seqInfo %>%
      dplyr::group_by(HR_Assays_Affected) %>%
      dplyr::summarise(All = dplyr::n()) %>%
      dplyr::rename(Assays_Affected = s$BP10_AA) %>%
      tidyr::pivot_wider(names_from = s$ASSAYS_AFFECTED, values_from = All)

    seqInfoC_HR <- seqInfo %>%
      dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group, HR_Assays_Affected) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::ungroup()

    df <- seqInfoC_HR %>%
      dplyr::rename(Assays_Affected = s$BP10_AA) %>%
      tidyr::pivot_wider(names_from = s$ASSAYS_AFFECTED, values_from = count) %>%
      dplyr::mutate_all(~ replace(., is.na(.), 0))
  } else { # all mutations
    total <- seqInfo %>%
      dplyr::group_by(Assays_Affected) %>%
      dplyr::summarise(All = dplyr::n()) %>%
      tidyr::pivot_wider(names_from = s$ASSAYS_AFFECTED, values_from = All)

    seqInfo <- seqInfo %>%
      dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group, Assays_Affected) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::ungroup()

    df <- seqInfo %>%
      tidyr::pivot_wider(names_from = s$ASSAYS_AFFECTED, values_from = count) %>%
      dplyr::mutate_all(~ replace(., is.na(.), 0))
  }

  # Order numbered columns
  non_numbered <- df %>% dplyr::select(c(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group))
  numbered <- df %>%
    dplyr::select(-c(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group)) %>%
    dplyr::select(sort(names(.)))
  df <- cbind(non_numbered, numbered)

  df$Displayed_Lineage[df$Displayed_Lineage == 0] <- NA

  df <- df %>%
    dplyr::left_join(lin_df, by = c("WHO_Designation", "Subgroup", "Displayed_Lineage", "Displayed_Group")) %>%
    rbind(c("All", " ", " ", "Summary: All Sequences by Assays", unlist(total[1, ]), total_seqs))

  # Add columns 5-7 if needed
  for (col in 0:7) {
    if (!(col %in% colnames(df))) {
      df <- df %>% tibble::add_column(!!toString(col) := 0, .before = "Total # Sequences")
    }
  }

  # Make number columns double type
  df <- df %>% dplyr::mutate(dplyr::across(`0`:`7`, as.numeric))

  # Combine 6 and 7 assay columns
  df <- df %>%
    dplyr::mutate(new_col = `6` + `7`) %>%
    dplyr::rename(!!s$NEW_COL_NAME := new_col) %>%
    dplyr::select(-c(`6`, `7`))


  # -----------------------------------------------------------------------------
  # 7-assay .1% test (for conclusion slide in SPC presentations) ---------------

  if (!is.null(startDate) && !is.null(endDate)) {
    # get frequency change for sequences with mutations affecting 5 or more assays from last period to this one and save (for conclusion slide in SPC presentation)

    total_5_or_more <- total %>% dplyr::select(c(`5`:tidyselect::last_col()))
    total_5_or_more <- rowSums(total_5_or_more)

    curr_prop <- total_5_or_more / total_seqs


    # get past prop (previous time period does not overlap current)
    monthDiff <- lubridate::interval(lubridate::ym(startDate), lubridate::ym(endDate)) %/% months(1)
    last_end_date <- lubridate::ym(startDate) - months(1)
    last_start_date <- lubridate::ym(startDate) - months(1 + monthDiff)
    last_end_date <- toString(format(last_end_date, "%Y-%m"))
    last_start_date <- toString(format(last_start_date, "%Y-%m"))

    con <- DBI::dbConnect(RSQLite::SQLite(), archiveDBFile)
    allSeqInfo <- dplyr::tbl(con, "Seq_Info") %>% dplyr::collect()
    DBI::dbDisconnect(con)

    pastSeqInfo <- allSeqInfo %>% dplyr::filter(Collection_Date <= last_end_date & Collection_Date >= last_start_date)

    pastTotalSeqs <- pastSeqInfo %>%
      dplyr::summarise(total_seqs = dplyr::n()) %>%
      dplyr::pull(total_seqs)


    if (isTRUE(HR)) {
      pastTotal <- pastSeqInfo %>%
        dplyr::group_by(HR_Assays_Affected) %>%
        dplyr::summarise(All = dplyr::n()) %>%
        dplyr::rename(Assays_Affected = s$BP10_AA) %>%
        tidyr::pivot_wider(names_from = s$ASSAYS_AFFECTED, values_from = All)
    } else {
      pastTotal <- pastSeqInfo %>%
        dplyr::group_by(Assays_Affected) %>%
        dplyr::summarise(All = dplyr::n()) %>%
        tidyr::pivot_wider(names_from = s$ASSAYS_AFFECTED, values_from = All)
    }

    past_total_5_or_more <- pastTotal %>% dplyr::select(c(`5`:tidyselect::last_col()))
    past_total_5_or_more <- rowSums(past_total_5_or_more)

    past_prop <- past_total_5_or_more / pastTotalSeqs

    prop_change <- curr_prop - past_prop
  }

  # -----------------------------------------------------------------------------

  return(list(df = df, prop_change = prop_change))
}



#' Make Co-occurring Mutations Dataframe V2
#'
#' @param in_df input dataframe, taken from getCoOccCounts function
#'
#' @return out_df, a dataframe ready to be made into a GT table
#'
makeDFcoOcc_2 <- function(in_df) {
  `0` <- `7` <- `Total # Sequences` <- `Total_%` <- `WHO_Designation` <- `WHO Variant Designation` <- Displayed_Group <- Displayed_Lineage <- WHO_Designation <- `0_%` <-
    `1` <- `1_%` <- `2` <- `2_%` <- `3` <- `3_%` <- `4` <- `4_%` <- `5` <- `5_%` <-
    `6` <- `6_%` <- `7_%` <- NULL


  out_df <- in_df %>%
    dplyr::mutate(dplyr::across(5:ncol(in_df), as.double)) %>%
    dplyr::ungroup()

  # Add percent columns
  out_df <- out_df %>%
    dplyr::mutate_at(dplyr::vars(`0`:!!s$NEW_COL_NAME), list("%" = ~ . / `Total # Sequences`))

  # Temporarily remove bottom row and round it to 1 decimal place (for percentages)
  bottom_totals <- utils::tail(out_df, n = 1) %>% dplyr::mutate_if(is.numeric, ~ round(., 3))

  out_df <- utils::head(out_df, -1)

  # Round everything else to nearest integer (for percentages)
  out_df <- out_df %>% dplyr::mutate_if(is.numeric, ~ round(., 2))

  # Get total number of sequences
  overall_total <- bottom_totals %>%
    dplyr::select(`Total # Sequences`) %>%
    dplyr::pull()

  # Add bottom row back
  out_df <- out_df %>% rbind(bottom_totals)

  # Add Total # sequences value to bottom row
  out_df$`Total # Sequences`[is.na(out_df$`Total # Sequences`)] <- overall_total

  # Add the last percent column, which is the totals for each variant divided by the overall total
  # Round the last column to 1 decimal place (for percentages)
  out_df <- out_df %>%
    dplyr::mutate(`Total_%` = `Total # Sequences` / overall_total) %>%
    dplyr::mutate(dplyr::across(`Total_%`, round, 3))

  # Change "Other" to empty string (Since this will be the only row in the row group, it doesn't need a name)
  out_df <- within(out_df, Subgroup[WHO_Designation == "Other"] <- "")

  out_df["Displayed_Lineage"][is.na(out_df["Displayed_Lineage"])] <- ""

  # Order columns correctly (count, percentage, count, percentage, etc)
  out_df <- out_df %>%
    dplyr::select(
      `Displayed_Group`, `Displayed_Lineage`, `WHO_Designation`, `0`, `0_%`, `1`, `1_%`, `2`, `2_%`, `3`, `3_%`, `4`, `4_%`, `5`,
      `5_%`, !!s$NEW_COL_NAME, !!s$NEW_PERCENT_COL_NAME, `Total # Sequences`, `Total_%`
    )

  # Make displayed lin match who designation for other and all values
  out_df <- out_df %>%
    dplyr::mutate(Displayed_Lineage = ifelse(WHO_Designation == "Other",
      WHO_Designation,
      Displayed_Lineage
    ))
  out_df <- out_df %>%
    dplyr::mutate(Displayed_Lineage = ifelse(WHO_Designation == "All",
      WHO_Designation,
      Displayed_Lineage
    ))


  # for gt table formatting
  out_df[out_df == "Other"] <- ""
  out_df[out_df == "All"] <- " "


  # Rename some columns for gt table
  out_df <- out_df %>% dplyr::rename("rowgrp_col" = Displayed_Group, "Pangolin Lineage" = Displayed_Lineage, "WHO Variant Designation" = WHO_Designation, "Sequences by Lineage" = `Total # Sequences`)


  out_df[is.na(out_df)] <- 0

  # Remove WHO Variant Designation column
  out_df <- out_df %>% dplyr::select(-`WHO Variant Designation`)

  return(out_df)
}


#' Make Tables for Co-occurring mutations (Condensed Version)
#'
#' @description Make 2 tables (one version with all mutated sequences and one version with sequences with high-risk mutations) showing mutated sequence counts categorized by number of assays affected on top and by variant on the left
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param variantFile csv file containing "Displayed_Lineage_Order" column for ordering row names/row groups in the GT table
#' @param startDate start date for subset DB, if NULL then the .1 percent analysis won't be run
#' @param endDate start date for subset DB, if NULL then the .1 percent analysis won't be run
#' @param archiveDBFile File path to archive db (used to retrieve older subset of data for comparison; older subset uses same length of time as input subset db, e.g. start_date is "2022-5" and end_date is "2022-07", then older subset start date is "2022-02 and "2022-04)
#' @param saveTable logical value indicating whether to save the tables to the specified saveFolder
#' @param saveFolder file path containing folder the finished tables will be saved in
#' @param height desired height of table .png images
#' @param width desired width of table .png images; if the width is not large enough, text will wrap or table will be cut off
#' @param colWidth1 width of "Pangolin Lineage" column in px; if NULL, widths are automatically assigned to try to fit all text
#' @param colWidth2 width of "5" column in px; if NULL, widths are automatically assigned to try to fit all text
#' @param colWidth3 width of ">=6" column in px; if NULL, widths are automatically assigned to try to fit all text
#' @param colWidth4 width of "Sequences by Lineage" column in px; if NULL, widths are automatically assigned to try to fit all text
#' @param colWidth5 width of all other columns in px; if NULL, widths are automatically assigned to try to fit all text
#'
#' @export
#'
makeCoOccTablesCondensed <- function(DB, variantFile, startDate = NULL, endDate = NULL, archiveDBFile = NULL, saveTable = FALSE, saveFolder = NULL, height = NULL, width = NULL, colWidth1 = NULL, colWidth2 = NULL, colWidth3 = NULL, colWidth4 = NULL, colWidth5 = NULL) {
  `Pangolin Lineage` <- Subgroup <- NULL

  setStrings2AssayDB()

  res <- getCoOccCounts_2(DB, HR = FALSE, startDate, endDate, archiveDBFile)
  co_occ_all <- res$df
  co_occ_all <- makeDFcoOcc_2(co_occ_all)

  hr_res <- getCoOccCounts_2(DB, HR = TRUE, startDate, endDate, archiveDBFile)
  co_occ_hr <- hr_res$df
  co_occ_hr <- makeDFcoOcc_2(co_occ_hr)


  # Order rows based on Displayed Lineage Order column in lineage input file
  WHO_lin_des <- readr::read_csv(variantFile, show_col_types = FALSE)
  rowname_order <- WHO_lin_des$Displayed_Lineage_Order
  rowname_order <- rowname_order[!is.na(rowname_order)]
  rowname_order <- c(rowname_order, " ")

  co_occ_all <- co_occ_all %>%
    dplyr::mutate(`Pangolin Lineage` = forcats::fct_relevel(`Pangolin Lineage`, rowname_order)) %>%
    dplyr::arrange(`Pangolin Lineage`)

  co_occ_hr <- co_occ_hr %>%
    dplyr::mutate(`Pangolin Lineage` = forcats::fct_relevel(`Pangolin Lineage`, rowname_order)) %>%
    dplyr::arrange(`Pangolin Lineage`)


  col_names <- c("0", "1", "2", "3", "4", "5", {{ s$NEW_COL_NAME }}, "Sequences by Lineage")
  col_names_pct <- c("0_%", "1_%", "2_%", "3_%", "4_%", "5_%", {{ s$NEW_PERCENT_COL_NAME }}, "Total_%")



  # write to proportion change file
  prop_df <- data.frame(all_prop_change = res$prop_change, hr_prop_change = hr_res$prop_change)

  if (!dir.exists(saveFolder)) {
    dir.create(saveFolder, recursive = TRUE)
  }

  readr::write_csv(
    prop_df,
    paste0(saveFolder, "7_assay_prop_change_", startDate, "_to_", endDate, ".csv")
  )



  # all mutations --------------------------------------------------------------------

  # Change the percent columns to percent format
  gt_tbl <- gt::gt(co_occ_all, groupname_col = "rowgrp_col") %>%
    gt::fmt_percent(
      columns = gtExtras::all_of(col_names_pct),
      rows = gt::everything(),
      drop_trailing_zeros = TRUE
    )

  # Merge count and percent columns
  for (i in seq_along(col_names)) {
    gt_tbl <- gt_tbl %>% gt::cols_merge_n_pct(
      col_n = col_names[i],
      col_pct = col_names_pct[i],
      autohide = TRUE
    )
  }

  # Add title
  gt_tbl <- gt_tbl %>%
    gt::tab_header(
      title = "Frequencies of Sequences Containing Identifiable Co-occurring Mutations"
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bolder"),
      location = gt::cells_title()
    )

  # Add formatting
  gt_tbl <- format_co_gt_tbl_2(gt_tbl, colWidth1, colWidth2, colWidth3, colWidth4, colWidth5, HR = FALSE)

  print(gt_tbl)

  if (saveTable) {
    if (!dir.exists(saveFolder)) {
      dir.create(saveFolder, recursive = TRUE)
    }

    gt_tbl %>% gt::gtsave(paste0(saveFolder, "co_occ_all.png"), vwidth = width)
  }



  # 3 prime end mutations  -----------------------------------------------------------------

  # Change the percent columns to percent format
  gt_tbl_hr <- gt::gt(co_occ_hr, groupname_col = "rowgrp_col") %>%
    gt::fmt_percent(
      columns = col_names_pct,
      rows = gt::everything(),
      drop_trailing_zeros = TRUE
    )

  # Merge count and percent columns
  for (i in seq_along(col_names)) {
    gt_tbl_hr <- gt_tbl_hr %>% gt::cols_merge_n_pct(
      col_n = col_names[i],
      col_pct = col_names_pct[i],
      autohide = TRUE
    )
  }


  # Add title
  gt_tbl_hr <- gt_tbl_hr %>%
    gt::tab_header(
      title = "Frequencies of Sequences Containing Identifiable 3' End Co-occurring Mutations"
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bolder"),
      location = gt::cells_title()
    )

  # Add formatting
  gt_tbl_hr <- format_co_gt_tbl_2(gt_tbl_hr, colWidth1, colWidth2, colWidth3, colWidth4, colWidth5, HR = TRUE)

  print(gt_tbl_hr)

  if (saveTable) {
    if (!dir.exists(saveFolder)) {
      dir.create(saveFolder, recursive = TRUE)
    }

    gt_tbl_hr %>% gt::gtsave(paste0(saveFolder, "co_occ_hr.png"), vwidth = width)
  }
}
