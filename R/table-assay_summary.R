###### Make Assay Summary Tables for COVID Inclusivity Reports ######


# ---- Helper functions ----
color_good_cells <- function(cols, condition) {
  gt::cells_body(
    columns = !!rlang::sym(cols),
    rows = !!rlang::sym(cols) / !!rlang::sym(condition) < 0.01
  )
}


color_bad_cells <- function(cols, condition) {
  gt::cells_body(
    columns = !!rlang::sym(cols),
    rows = !!rlang::sym(cols) / !!rlang::sym(condition) > 0.05
  )
}


#' Format Assay Summary GT Table
#'
#' @param in_tbl GT table object
#' @param colWidth1 width for "Pangolin Lineage" and "WHO Label" columns
#' @param colWidth2 width for all other columns
#'
#' @return formatted GT table
#' @noRd
#'
format_as_gt_tbl <- function(in_tbl, colWidth1 = NULL, colWidth2 = NULL) {
  `Pangolin Lineage` <- `Sequences by Lineage` <- `2h` <- `7` <- where <- NULL

  out_tbl <- in_tbl %>%
    gt::tab_spanner(
      label = "Mutated Sequences by Assay",
      columns = c("2a", "2c", "2d", "2e", "2f", "2g", "2h")
    ) %>%
    gt::tab_spanner(
      label = "Co-occurring Mutated Sequences",
      columns = c("5", "6", "7")
    ) %>%
    gt::tab_footnote(
      footnote = "Mutated sequences co-occurring on multiple assays, by # assays",
      locations = gt::cells_column_spanners(
        spanners = "Co-occurring Mutated Sequences"
      )
    ) %>%
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

  # Add color formatting, optional additional source notes, and note explaining the formatting
  cols <- c("2a", "2c", "2d", "2e", "2f", "2g", "2h", "5", "6", "7")
  last_col <- c("Sequences by Lineage")

  out_tbl <- out_tbl %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#CDE9EB"), # light blue
      locations = lapply(cols, color_good_cells, condition = "Sequences by Lineage")
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#FFE699"), # light yellow
      locations = lapply(cols, color_bad_cells, condition = "Sequences by Lineage")
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#D5D8DC"), # grey
      locations = gt::cells_body(columns = `Sequences by Lineage`)
    )

  # add source notes
  out_tbl <- out_tbl %>%
    gt::tab_source_note( # use paste0 to print >= with no space between it and the number
      source_note = paste0("Color guide - Yellow: mutations occur at ", "\U2265", "5%; Blue: mutations occur at <1%")
    ) %>%
    gt::tab_source_note(
      source_note = "Note: Assay names are abbreviated as SARS-COV-2a (2a), SARS-COV-2d (2d), etc."
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
      columns = cols,
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
      columns = c(`2h`, `7`),
      sides = "right",
      color = "grey",
      weight = gt::px(3)
    )

  if (!is.null(colWidth1)) {
    out_tbl <- out_tbl %>% gt::cols_width(
      c(`Pangolin Lineage`) ~ gt::px(colWidth1) # 110
    )
  }

  if (!is.null(colWidth2)) {
    out_tbl <- out_tbl %>% gt::cols_width(
      c(`2a`, `2c`, `2d`, `2e`, `2f`, `2g`, `2h`, `5`, `6`, `7`, `Sequences by Lineage`) ~ gt::px(colWidth2) # 130
    )
  }


  out_tbl <- out_tbl %>% gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = list(
      gt::cells_row_groups(groups = "Summary: All Sequences by Assays"),
      gt::cells_body(columns = gt::everything(), rows = `Pangolin Lineage` == " ")
    ) # the All row
  )

  return(out_tbl)
}




# ---- Main functions ----


#' Get Assay Summary Counts
#'
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param HR if TRUE, returns 3 prime end, i.e. high risk, version
#' @param BASE_PAIR_LIMIT number of base pairs from the 3 prime end to be considered high risk
#'
#' @return single_assays, a dataframe with mutated sequence counts per assay
#' @noRd
#'
getAssaySummaryCounts <- function(DB, HR = FALSE, BASE_PAIR_LIMIT = 10) {
  `Assay` <- `All` <- `WHO_Designation` <- `Subgroup` <- `Lineage` <- `Displayed_Lineage` <- `Displayed_Group` <- `Counts` <- `1` <- `2` <-
    `3` <- `4` <- `5` <- `6` <- `7` <- NULL

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  idMap <- dplyr::tbl(con, "ID_Map") %>% dplyr::collect()
  varType <- dplyr::tbl(con, "Variant_Type") %>% dplyr::collect()
  primers <- dplyr::tbl(con, "Primer_Info") %>% dplyr::collect()

  DBI::dbDisconnect(con)

  primer_ids <- primers %>% dplyr::select(s$PRIMER_ID, s$ASSAY_ID)
  assayNames <- primers %>%
    dplyr::mutate(Assay = substr(get(s$ASSAY_NAME), (nchar(get(s$ASSAY_NAME)) - 1), nchar(get(s$ASSAY_NAME)))) %>%
    dplyr::select(s$ASSAY_ID, Assay) %>%
    dplyr::distinct()

  idMap <- idMap %>% dplyr::select(-c(Lineage))

  if (HR) { # mutations on the 3 prime end
    VarIds3P <- varType %>%
      dplyr::filter(get(s$POS_FROM_SIDE) <= BASE_PAIR_LIMIT) %>%
      dplyr::select(s$VAR_ID) %>%
      dplyr::distinct()
    idMap <- VarIds3P %>% dplyr::inner_join(idMap, by = "Var_ID")
  }

  total_counts <- idMap %>%
    dplyr::select(-c(s$VAR_ID, s$PRIMER_ID)) %>%
    dplyr::distinct() %>%
    dplyr::group_by("Assay_ID" = get(s$ASSAY_ID)) %>%
    dplyr::summarise(`All` = dplyr::n()) %>%
    dplyr::inner_join(assayNames, by = "Assay_ID") %>%
    dplyr::select(Assay, `All`) %>%
    tidyr::pivot_wider(names_from = Assay, values_from = All)

  lin_counts <- idMap %>%
    dplyr::select(-c(s$VAR_ID, s$PRIMER_ID)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group, "Assay_ID" = get(s$ASSAY_ID)) %>%
    dplyr::summarise(Counts = dplyr::n()) %>%
    tidyr::pivot_wider(names_from = s$ASSAY_ID, values_from = Counts) %>%
    dplyr::mutate_all(~ replace(., is.na(.), 0)) %>%
    dplyr::select(c(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group, `1`, `2`, `3`, `4`, `5`, `6`, `7`)) %>%
    dplyr::rename_at(dplyr::vars(`1`, `2`, `3`, `4`, `5`, `6`, `7`), ~ dplyr::pull(assayNames, Assay))

  single_assays <- lin_counts %>% rbind(total_counts)

  single_assays["WHO_Designation"][is.na(single_assays["WHO_Designation"])] <- "All"

  # No subcategory is needed for the "All" and "Other" groups; set to empty strings
  single_assays[is.na(single_assays)] <- ""
  single_assays["Subgroup"][single_assays["Subgroup"] == "Other"] <- ""


  single_assays["Displayed_Group"][single_assays["WHO_Designation"] == "All"] <- "Summary: All Sequences by Assays"

  return(single_assays)
}



#' Get Co-occurring Counts for 5 to 7 Assays
#'
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param HR if TRUE, returns 3 prime end, i.e. high risk, version
#'
#' @return df_5to7, a dataframe with mutated sequence counts co-occurring on multiple assays,
#' only including rows for 5-7 assays affected and the totals
#' @noRd
#'
getCoOccCounts5to7 <- function(DB, HR = FALSE) {
  WHO_Designation <- Subgroup <- Lineage <- Displayed_Lineage <- Displayed_Group <- HR_Assays_Affected <- Assays_Affected <- All <-
    count <- `:=` <- `5` <- `6` <- `7` <- `Total # Sequences` <- `.` <- NULL

  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  # Tables with Needed Information
  seqInfo <- dplyr::tbl(con, "Seq_Info") %>% dplyr::collect()
  assayInfo <- dplyr::tbl(con, "Assay_Info") %>%
    dplyr::select(s$ASSAY_ID, s$ASSAY_NAME) %>%
    dplyr::collect()

  DBI::dbDisconnect(con)

  total_seqs <- seqInfo %>%
    dplyr::summarise(total_seqs = dplyr::n()) %>%
    dplyr::collect() %>%
    dplyr::pull(total_seqs)

  seqInfo <- seqInfo %>% dplyr::select(-c(Lineage))


  lin_count <- seqInfo %>%
    dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group) %>%
    dplyr::summarise(`Total # Sequences` = dplyr::n()) %>%
    dplyr::collect() %>%
    dplyr::ungroup()
  lin_count <- rbind(lin_count, c(WHO_Designation = "All", Subgroup = "All", Displayed_Lineage = "All", Displayed_Group = "All", `Total # Sequences` = total_seqs))

  lin_df <- lin_count %>% dplyr::mutate_all(~ replace(., is.na(`Total # Sequences`), 0))



  if (HR) { # mutations on the 3 prime end
    # TODO: there is something wrong here--counts for 5+ assays are off (only HR version)
    total <- seqInfo %>%
      dplyr::group_by(HR_Assays_Affected) %>%
      dplyr::summarise(All = dplyr::n()) %>%
      dplyr::collect() %>%
      dplyr::rename(Assays_Affected = HR_Assays_Affected) %>%
      tidyr::pivot_wider(names_from = Assays_Affected, values_from = All)

    seqInfoC_HR <- seqInfo %>%
      dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group, HR_Assays_Affected) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::collect() %>%
      dplyr::ungroup()

    df <- seqInfoC_HR %>%
      dplyr::rename(Assays_Affected = HR_Assays_Affected) %>%
      tidyr::pivot_wider(names_from = Assays_Affected, values_from = count) %>%
      dplyr::mutate_all(~ replace(., is.na(.), 0)) %>%
      dplyr::relocate(sort(names(df)))

  } else { # all mutations
    total <- seqInfo %>%
      dplyr::group_by(Assays_Affected) %>%
      dplyr::summarise(All = dplyr::n()) %>%
      dplyr::collect() %>%
      tidyr::pivot_wider(names_from = Assays_Affected, values_from = All)

    seqInfo <- seqInfo %>%
      dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group, Assays_Affected) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::collect() %>%
      dplyr::ungroup()

    df <- seqInfo %>%
      tidyr::pivot_wider(names_from = Assays_Affected, values_from = count) %>%
      dplyr::mutate_all(~ replace(., is.na(.), 0))
  }

  # put columns in right order-----------------------------------------------
  # separate word columns from number columns
  non_numbered_names <- c("WHO_Designation", "Subgroup", "Displayed_Lineage", "Displayed_Group")
  non_numbered <- df %>% dplyr::select(dplyr::all_of(non_numbered_names))
  numbered <- df %>%
    dplyr::select(-non_numbered_names) %>%
    dplyr::select(sort(names(.)))
  # add back together
  df <- cbind(non_numbered, numbered)

  df$Displayed_Lineage[df$Displayed_Lineage == 0] <- NA

  # for this version of lin_df, remove bottom row so that we can left join with df
  lin_df_without_totals <- lin_df %>% dplyr::filter(dplyr::row_number() <= dplyr::n() - 1)

  df <- df %>% dplyr::full_join(lin_df_without_totals, by = c("WHO_Designation", "Subgroup", "Displayed_Lineage", "Displayed_Group"))

  df <- df %>% rbind(c("All", " ", " ", "Summary: All Sequences by Assays", unlist(total[1, ]), total_seqs))


  # Add columns 5-7 if needed
  for (col in 0:7) {
    if (!(col %in% colnames(df))) {
      df <- df %>% tibble::add_column(!!toString(col) := 0, .before = "Total # Sequences")
    }
  }


  df_5to7 <- df %>% dplyr::select(`5`, `6`, `7`, `Total # Sequences`)
  df_5to7[is.na(df_5to7)] <- ""

  return(df_5to7)
}



#' Make Combined Dataframe
#'
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param HR if TRUE, returns 3 prime end, i.e. high risk, version
#'
#' @return combined_df, a dataframe with counts and percentages both by assay and
#' by co-occurring assays
#' @noRd
#'
makeCombinedDataframe <- function(DB, HR = FALSE) {
  `2a` <- `7` <- `Total # Sequences` <- `Total_%` <- `2a_%` <- `2c` <- `2c_%` <-
    `2d` <- `2d_%` <- `2e` <- `2e_%` <- `2f` <- `2f_%` <- `2g` <- `2g_%` <- `2h` <-
    `2h_%` <- `5` <- `5_%` <- `6` <- `6_%` <- `7_%` <- `WHO_Designation` <-
    `Subgroup` <- `Lineage` <- `LinSubgroup` <- Displayed_Group <- Displayed_Lineage <- NULL

  if (HR) {
    single <- getAssaySummaryCounts(DB, HR = TRUE)
    co_oc <- getCoOccCounts5to7(DB, HR = TRUE)
  } else {
    single <- getAssaySummaryCounts(DB)
    co_oc <- getCoOccCounts5to7(DB) # issue here
  }

  num_cols <- ncol(co_oc)
  co_oc <- co_oc %>%
    dplyr::mutate(dplyr::across(1:num_cols, as.double))

  combined_df <- cbind(single, co_oc) %>% dplyr::mutate_all(~ replace(., is.na(.), 0))

  # Add percent columns
  combined_df <- combined_df %>%
    dplyr::mutate_at(dplyr::vars(`2a`:`7`), list("%" = ~ . / `Total # Sequences`))

  # Temporarily remove bottom row and round it to 1 decimal place (for percentages)
  bottom_totals <- utils::tail(combined_df, n = 1) %>% dplyr::mutate_if(is.numeric, ~ round(., 3))

  combined_df <- utils::head(combined_df, -1)

  # Round everything else to nearest integer (for percentages)
  combined_df <- combined_df %>% dplyr::mutate(dplyr::across(dplyr::everything(), round, 2))

  # Get total number of sequences
  overall_total <- bottom_totals %>%
    dplyr::select(`Total # Sequences`) %>%
    dplyr::pull()

  # Add bottom row back
  combined_df <- combined_df %>% rbind(bottom_totals)

  # Add Total # sequences value to bottom row
  combined_df$`Total # Sequences`[is.na(combined_df$`Total # Sequences`)] <- overall_total

  # Add the last percent column, which is the totals for each variant divided by the overall total
  # Round the last column to 1 decimal place (for percentages)
  combined_df <- combined_df %>%
    dplyr::mutate(`Total_%` = `Total # Sequences` / overall_total) %>%
    dplyr::mutate(dplyr::across(`Total_%`, round, 3))

  combined_df <- combined_df %>% dplyr::ungroup()

  # Order columns correctly (count, percentage, count, percentage, etc)
  combined_df <- combined_df %>% dplyr::select(
    `Displayed_Group`, `Displayed_Lineage`, `WHO_Designation`, `2a`, `2a_%`, `2c`, `2c_%`, `2d`, `2d_%`, `2e`, `2e_%`, `2f`, `2f_%`, `2g`,
    `2g_%`, `2h`, `2h_%`, `5`, `5_%`, `6`, `6_%`, `7`, `7_%`, `Total # Sequences`, `Total_%`
  )


  # Make displayed lin match who designation for other and all values
  combined_df <- combined_df %>%
    dplyr::mutate(Displayed_Lineage = ifelse(WHO_Designation == "Other",
      WHO_Designation,
      Displayed_Lineage
    ))
  combined_df <- combined_df %>%
    dplyr::mutate(Displayed_Lineage = ifelse(WHO_Designation == "All",
      WHO_Designation,
      Displayed_Lineage
    ))


  # for gt table formatting
  combined_df[combined_df == "Other"] <- ""
  combined_df[combined_df == "All"] <- " "


  # Rename some columns for gt table
  combined_df <- combined_df %>% dplyr::rename("rowgrp_col" = Displayed_Group, "Pangolin Lineage" = Displayed_Lineage, "WHO Label" = WHO_Designation, "Sequences by Lineage" = `Total # Sequences`)


  combined_df[is.na(combined_df)] <- 0

  return(combined_df)
}



#' Make Assay Summary Tables
#'
#' @description Make 2 tables, a version with all mutations from the subset database, and a version with high-risk mutations; displays mutated sequence counts categorized by assay and number of assays affected on top and by variant on the left
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param variantFile csv file containing "Displayed_Lineage_Order" column for ordering row names/row groups in the GT table
#' @param extraNotes vector of strings that are extra source notes to be displayed at the bottom of the table
#' @param saveTable logical value indicating whether to save the tables to the specified saveFolder
#' @param saveFolder file path containing folder the finished tables will be saved in
#' @param height desired height of table .png images
#' @param width desired width of table .png images; if the width is not large enough, text will wrap or table will be cut off
#' @param colWidth1 width of "Pangolin Lineage" column; if NULL, widths are automatically assigned to try to fit all text
#' @param colWidth2 width of all other columns; if NULL, widths are automatically assigned to try to fit all text
#'
#' @export
#'
makeAssaySummaryTables <- function(DB, variantFile, extraNotes = NULL, saveTable = FALSE, saveFolder = NULL, height = NULL, width = NULL, colWidth1 = NULL, colWidth2 = NULL) {
  setStrings2AssayDB()

  `Pangolin Lineage` <- Subgroup <- `Total # Sequences` <- `2h` <- `7` <- NULL

  combined <- makeCombinedDataframe(DB)
  combined_hr <- makeCombinedDataframe(DB, HR = TRUE)


  # No need for WHO Label column anymore
  combined <- combined %>% dplyr::select(-c(3))
  combined_hr <- combined_hr %>% dplyr::select(-c(3))


  # Order rows based on Displayed Lineage Order column in lineage input file
  WHO_lin_des <- readr::read_csv(variantFile, show_col_types = FALSE)
  rowname_order <- WHO_lin_des$Displayed_Lineage_Order
  rowname_order <- rowname_order[!is.na(rowname_order)]
  rowname_order <- c(rowname_order, " ")

  combined <- combined %>%
    dplyr::mutate(`Pangolin Lineage` = forcats::fct_relevel(`Pangolin Lineage`, rowname_order)) %>%
    dplyr::arrange(`Pangolin Lineage`)

  combined_hr <- combined_hr %>%
    dplyr::mutate(`Pangolin Lineage` = forcats::fct_relevel(`Pangolin Lineage`, rowname_order)) %>%
    dplyr::arrange(`Pangolin Lineage`)


  col_names <- c("2a", "2c", "2d", "2e", "2f", "2g", "2h", "5", "6", "7", "Sequences by Lineage")
  col_names_pct <- c("2a_%", "2c_%", "2d_%", "2e_%", "2f_%", "2g_%", "2h_%", "5_%", "6_%", "7_%", "Total_%")

  # make GT table for all mutations---------------------------------------------------------------------
  gt_tbl <- gt::gt(combined, groupname_col = "rowgrp_col") %>%
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
      title = "Frequencies of Sequences Containing Identifiable Mutations Across Multiple Assays"
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bolder"),
      location = gt::cells_title()
    )


  formatted_gt_tbl <- format_as_gt_tbl(gt_tbl, colWidth1, colWidth2)

  print(formatted_gt_tbl)

  if (saveTable) {
    if (!dir.exists(saveFolder)) {
      dir.create(saveFolder, recursive = TRUE)
    }

    formatted_gt_tbl %>% gt::gtsave(paste0(saveFolder, "assay_summary_all.png"), vwidth = width)
  }






  # make GT table for high risk mutations--------------------------------------------------------------


  gt_tbl_hr <- gt::gt(combined_hr, groupname_col = "rowgrp_col") %>%
    gt::fmt_percent(
      columns = gtExtras::all_of(col_names_pct),
      rows = gt::everything(),
      drop_trailing_zeros = TRUE
    )


  for (i in seq_along(col_names)) {
    gt_tbl_hr <- gt_tbl_hr %>% gt::cols_merge_n_pct(
      col_n = col_names[i],
      col_pct = col_names_pct[i],
      autohide = TRUE
    )
  }


  # Add title and hr footnote
  gt_tbl_hr <- gt_tbl_hr %>%
    gt::tab_header(
      title = "Frequencies of Sequences with High-Risk Mutations"
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bolder"),
      location = gt::cells_title()
    ) %>%
    gt::tab_footnote(
      footnote = "Non-ambiguous mutations under the primer binding regions that fall
    within 10bp of 3' end of the primer were considered in this analysis.",
      locations = gt::cells_title()
    )

  formatted_gt_tbl_hr <- format_as_gt_tbl(gt_tbl_hr, colWidth1, colWidth2)

  print(formatted_gt_tbl_hr)

  if (saveTable) {
    formatted_gt_tbl_hr %>% gt::gtsave(file = paste0(saveFolder, "assay_summary_hr.png"), vwidth = width)
  }
}
