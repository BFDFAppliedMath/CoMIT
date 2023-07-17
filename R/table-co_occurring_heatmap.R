###### Make Co-Occurring Mutation Combinations Heatmap for COVID Inclusivity Reports ######


# ---- Helper functions ----

# https://jokergoo.github.io/2020/07/06/block-annotation-over-several-slices/
#' Add Group Block Annotation
#'
#' @param group the column groups this annotation will cover
#' @param empty_anno name of the anno_empty object created from HeatmapAnnotation
#' @param gp vector of graphic parameters for the rectangle
#' @param label text label for this annotation
#' @param label_gp vector of graphic parameters for the text
#'
group_block_anno <- function(group, empty_anno, gp = grid::gpar(),
                             label = NULL, label_gp = grid::gpar()) {
  grid::seekViewport(GetoptLong::qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 <- grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
  grid::seekViewport(GetoptLong::qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 <- grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

  grid::seekViewport("global")
  grid::grid.rect(loc1$x, loc1$y,
    width = loc2$x - loc1$x, height = loc2$y - loc1$y,
    just = c("left", "bottom"), gp = gp
  )
  if (!is.null(label)) {
    grid::grid.text(label, x = (loc1$x + loc2$x) * 0.5, y = (loc1$y + loc2$y) * 0.5, gp = label_gp)
  }
}


#' Get Three-Assay Proportion
#'
#' @description calculate proportion of sequences with mutations affecting all three of a, d, and e assays for post-market surveillance of 3-assay test
#'
#' @param coOccDF3Assay a dataframe returned by \code{makeDFcoOccHeatmap()}
#'
#' @return the number of sequences with mutations affecting all three assays divided by the total number of sequences in the dataset
#'
getThreeAssayProportion <- function(coOccDF3Assay) {
  totalSeqs <- coOccDF3Assay[c("Freq"), ] %>% rowSums()

  cols <- which(coOccDF3Assay["Assays_Affected", ] == 3)
  # all_three_affected <- coOccDF3Assay[, cols]
  all_three_affected <- coOccDF3Assay %>% dplyr::select(dplyr::all_of(cols))

  # if(dim(all_three_affected)[2] != 0){ # if there are rows, i.e. there are combinations of mutations where all three assays are affected
  all_three_affected$totals <- rowSums(all_three_affected)
  rownames(all_three_affected) <- NULL

  mutated_seqs_all_three_affected <- all_three_affected$totals[5]
  # }

  prop <- mutated_seqs_all_three_affected / totalSeqs

  return(prop)
}


# ---- Main functions ----


#' Make Dataframe for Co-occurring Mutations Heatmap
#'
#' @description Create dataframe to be used by \code{makeCoOccHeatmaps()} in creating visulization of co-occurring mutations on different combinations of assays, with the sequences organized according to WHO designations/Pango lineages
#'
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param HR if TRUE, returns 3 prime end, i.e. high risk, version
#' @param BASE_PAIR_LIMIT number of base pairs from the 3 prime end to be considered high risk
#' @param threeAssay if TRUE, makes three-assay version of dataframe (a, d, and e assays) for 3-assay version of Covid test
#' @param startDate start date for subset DB, if NULL then the .1 percent analysis won't be run
#' @param endDate end date for subset DB, if NULL then the .1 percent analysis won't be run
#' @param archiveDBFile File path to archive db (used to retrieve older subset of data for comparison; older subset uses same length of time as input subset db, e.g. start_date is "2022-5" and end_date is "2022-07", then older subset start date is "2022-02 and "2022-04)
#'
#' @return comboFreq, a dataframe with counts and percentages of each co-occurring mutation combo of assays
#'
makeDFcoOccHeatmap <- function(DB, HR = FALSE, BASE_PAIR_LIMIT = 10, threeAssay = FALSE, startDate = NULL, endDate = NULL, archiveDBFile = NULL) {
  Assays_Affected <- . <- Freq <- value <- Collection_Date <- Accession_ID <- Var_ID <- Assay_ID <- NULL

  if (is.null(startDate) && is.null(endDate)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), DB)
    # Tables with Needed Information
    idMap <- dplyr::tbl(con, "ID_Map")
    varType <- dplyr::tbl(con, "Variant_Type")
    primers <- dplyr::tbl(con, "Primer_Info")

    # Get Tables
    idMap <- idMap %>% dplyr::collect()
    varType <- varType %>%
      dplyr::select(s$VAR_ID, s$POS_FROM_SIDE) %>%
      dplyr::collect()
    primers <- primers %>%
      dplyr::select(s$PRIMER_ID, s$ASSAY_NAME, s$ASSAY_ID) %>%
      dplyr::collect()

    DBI::dbDisconnect(con)
  } else { # looking at past data
    monthDiff <- lubridate::interval(lubridate::ym(startDate), lubridate::ym(endDate)) %/% months(1)
    last_end_date <- lubridate::ym(startDate) - months(1)
    last_start_date <- lubridate::ym(startDate) - months(1 + monthDiff)
    last_end_date <- toString(format(last_end_date, "%Y-%m"))
    last_start_date <- toString(format(last_start_date, "%Y-%m"))

    con <- DBI::dbConnect(RSQLite::SQLite(), archiveDBFile)
    pastAccs <- dplyr::tbl(con, "Seq_Info") %>%
      dplyr::filter(Collection_Date <= last_end_date & Collection_Date >= last_start_date) %>%
      dplyr::collect() %>%
      dplyr::select(Accession_ID) %>%
      dplyr::pull()
    idMap <- dplyr::tbl(con, "ID_Map") %>%
      dplyr::filter(Accession_ID %in% pastAccs) %>%
      dplyr::collect()
    IDs <- idMap %>%
      dplyr::select(Var_ID) %>%
      dplyr::distinct() %>%
      dplyr::pull()
    varType <- dplyr::tbl(con, "Variant_Type") %>%
      dplyr::filter(Var_ID %in% IDs) %>%
      dplyr::select(s$VAR_ID, s$POS_FROM_SIDE) %>%
      dplyr::collect()
    primers <- dplyr::tbl(con, "Primer_Info") %>%
      dplyr::select(s$PRIMER_ID, s$ASSAY_NAME, s$ASSAY_ID) %>%
      dplyr::collect()

    DBI::dbDisconnect(con)
  }


  within10bp <- varType %>%
    dplyr::filter(get(s$POS_FROM_SIDE) <= BASE_PAIR_LIMIT) %>%
    dplyr::select(s$VAR_ID)

  if (isTRUE(HR)) {
    allCombos <- within10bp %>% dplyr::inner_join(idMap, by = "Var_ID")
  } else {
    allCombos <- idMap
  }

  if (isTRUE(threeAssay)) { # assay IDs for a, d, and e assays are 1, 3, and 4
    allCombos <- allCombos %>% dplyr::filter(Assay_ID == 1 | Assay_ID == 3 | Assay_ID == 4)
  }

  # Every row represents a sequence
  allCombos <- allCombos %>%
    dplyr::select(s$ACC_ID, s$ASSAY_ID) %>%
    dplyr::distinct() %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(
      names_from = s$ASSAY_ID, names_sort = TRUE,
      values_from = value, values_fill = 0
    ) %>%
    dplyr::select(-s$ACC_ID)

  # Aggregate by combination of assays affected and get frequencies of each combo
  comboFreq <- stats::aggregate(list(Freq = rep(1, nrow(allCombos))), allCombos, length)

  # Have to remove freq column to get # assays affected column
  withAssaysAffected <- comboFreq %>%
    dplyr::select(-Freq) %>%
    dplyr::mutate(Assays_Affected = rowSums(.))

  # Add assays affected back to combo freq
  comboFreq <- base::cbind(comboFreq, Assays_Affected = withAssaysAffected$Assays_Affected)

  comboFreq <- comboFreq %>% dplyr::arrange(Assays_Affected, dplyr::desc(Freq))


  if (isTRUE(threeAssay)) {
    names(comboFreq) <- c("2a", "2d", "2e", "Freq", "Assays_Affected")
  } else {
    names(comboFreq) <- c("2a", "2c", "2d", "2e", "2f", "2g", "2h", "Freq", "Assays_Affected")
  }


  # Move Assays_Affected to the front
  comboFreq <- comboFreq %>%
    dplyr::select(Assays_Affected, dplyr::everything())

  # Transpose
  comboFreq <- as.data.frame(t(comboFreq[, ]))

  return(comboFreq)
}


createMatrix <- function(inDF, totalSeqs, HR = FALSE) { # (dfHR, totalSeqs, HR=TRUE)
  # Remove Assays_Affected and Frequency rows from dataframe:

  assays_affected <- inDF["Assays_Affected", ]

  # assays_affected is used to split the columns into groups
  assays_affected <- as.numeric(assays_affected[1, ])

  # col_groups is used to create labels for column groups (# assays affected)
  col_groups <- unique(assays_affected)

  remove_row <- c("Assays_Affected")
  inDF <- inDF[!(row.names(inDF) %in% remove_row), ]

  counts <- utils::tail(inDF, n = 1)
  counts <- as.numeric(counts[1, ])

  # Add commas for heatmap display
  prettyCounts <- format(counts, big.mark = ",", trim = TRUE)
  prettyTotal <- format(totalSeqs, big.mark = ",", trim = TRUE)

  percents <- round(((counts / totalSeqs) * 100), 1)

  inDF <- utils::head(inDF, -1) # Remove Freq row

  outMatrix <- data.matrix(inDF, rownames.force = TRUE) # Convert to matrix for heatmap creation

  res <- list(
    outMatrix = outMatrix,
    col_groups = col_groups,
    assays_affected = assays_affected,
    counts = counts,
    prettyCounts = prettyCounts,
    percents = percents,
    prettyTotal = prettyTotal
  )

  return(res)
}



formatHeatmap <- function(inMatrix, col_groups, assays_affected, counts, prettyCounts, percents, prettyTotal, counts_font_size, HR = FALSE) {
  # Set colors for heatmap
  colors <- c(colorspace::sequential_hcl(5, palette = "Purple-Blue"))
  heatmapColors <- structure(c(colors[5], colors[2]), names = c("0", "1"))

  # For making "# Assays Affected" a group label instead of part of the title:
  top_labels <- ComplexHeatmap::HeatmapAnnotation(
    label = ComplexHeatmap::anno_empty(border = FALSE),
    foo = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = colors[3]),
      labels = col_groups
    )
  )

  count_bargraph <- ComplexHeatmap::HeatmapAnnotation(" " = ComplexHeatmap::anno_barplot(counts,
    add_numbers = TRUE,
    border = TRUE,
    height = grid::unit(40, "mm"),
    numbers_rot = 0
  ))

  # To use percentages instead of counts for bargraph:
  # percent_bargraph <- ComplexHeatmap::HeatmapAnnotation("% of Total Sequences" = ComplexHeatmap::anno_barplot(percents,
  #                                                          add_numbers=TRUE,
  #                                                          border=TRUE,
  #                                                          height=grid::unit(40, "mm"),
  #                                                          numbers_rot=0),
  #                                       annotation_name_side="left")

  percent_bargraph <- ComplexHeatmap::HeatmapAnnotation(
    "% of Total Sequences" = ComplexHeatmap::anno_barplot(percents,
      border = TRUE,
      height = grid::unit(40, "mm")
    ),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = 13)
  ) # fontsize for % of total seqs section label

  if (HR) {
    heatmap <- ComplexHeatmap::Heatmap(inMatrix,
      name = "hm2",
      col = heatmapColors,
      rect_gp = grid::gpar(col = "white", lwd = 2),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title = "Assay",
      row_title_side = "left",
      row_names_side = "left",
      column_title = "Co-occurring 3' End Mutation Combinations",
      column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
      column_title_side = "top",
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      column_labels = prettyCounts,
      column_names_side = "bottom",
      column_names_centered = TRUE,
      column_names_rot = 0,
      column_names_gp = grid::gpar(fontsize = counts_font_size), # fontsize for number labels of bargraph
      heatmap_legend_param = list(
        title = " ",
        labels = c("Assay Not Affected", "Assay Affected"),
        labels_gp = grid::gpar(fontsize = 12),
        border = "black"
      ),
      column_split = assays_affected,
      border = TRUE,
      top_annotation = top_labels,
      # bottom_annotation=count_bargraph
      bottom_annotation = percent_bargraph
    )
  } else {
    heatmap <- ComplexHeatmap::Heatmap(inMatrix,
      name = "hm1",
      col = heatmapColors,
      rect_gp = grid::gpar(col = "white", lwd = 2),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title = "Assay",
      row_title_side = "left",
      row_names_side = "left",
      column_title = "Co-occurring Mutation Combinations",
      column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
      column_title_side = "top",
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      column_labels = prettyCounts,
      column_names_side = "bottom",
      column_names_centered = TRUE,
      column_names_rot = 0,
      column_names_gp = grid::gpar(fontsize = counts_font_size), # fontsize for number labels of bargraph
      heatmap_legend_param = list(
        title = " ",
        labels = c("Assay Not Affected", "Assay Affected"),
        labels_gp = grid::gpar(fontsize = 12),
        border = "black"
      ),
      column_split = assays_affected,
      border = TRUE,
      top_annotation = top_labels,
      # bottom_annotation=count_bargraph
      bottom_annotation = percent_bargraph
    )
  }

  return(heatmap)
}


saveHeatmap <- function(saveTable, saveFolder, inHeatmap, HR = FALSE, prettyTotal, col_groups, w, h) {
  if (HR) {
    saveName <- "heatmap_hr.png"
    heatmapID <- "hm2"
  } else {
    saveName <- "heatmap_all.png"
    heatmapID <- "hm1"
  }

  if (saveTable) {
    if (!dir.exists(saveFolder)) {
      dir.create(saveFolder, recursive = TRUE)
    }

    grDevices::png(file = paste0(saveFolder, saveName), width = w, height = h, units = "px", res = 100)

    ComplexHeatmap::draw(inHeatmap, padding = grid::unit(c(10, 2, 2, 8), "mm")) # bottom, left, top, right paddings
  } else {
    ComplexHeatmap::draw(inHeatmap, padding = grid::unit(c(10, 2, 2, 2), "mm"))
  }

  # Total sequences text box at bottom of heatmap
  ComplexHeatmap::decorate_heatmap_body(heatmapID, {
    grid::grid.text(paste0("Total sequences: ", prettyTotal),
      grid::unit(22, "mm"), # 345 22
      grid::unit(-50, "mm"), # -20 -50
      gp = grid::gpar(fontsize = 11)
    )
  })

  # # assays affected header
  group_block_anno(1:length(col_groups), "label",
    gp = grid::gpar(col = "white"),
    label_gp = grid::gpar(fontsize = 14),
    label = "# Assays Affected"
  )

  if (saveTable) {
    grDevices::dev.off()
  }
}





#' Make Heatmaps for Co-occurring mutations
#'
#' @description Create heatmaps using the ComplexHeatmap package (one version with all mutated sequences and one version with sequences with high-risk mutations); heatmaps are colored according to whether an assay is affected by mutated sequences
#'
#' @param DB File path to subset database file. In the case of making summary tables, we use the
#' 3M version from CoMIT, returned from \code{pull_subset_DB()}
#' @param startDate start date for subset DB, if NULL then the .1 percent analysis won't be run
#' @param endDate end date for subset DB, if NULL then the .1 percent analysis won't be run
#' @param archiveDBFile File path to archive db (used to retrieve older subset of data for comparison; older subset uses same length of time as input subset db, e.g. start_date is "2022-5" and end_date is "2022-07", then older subset start date is "2022-02 and "2022-04)
#' @param saveTable logical value indicating whether to save the tables to the specified saveFolder;
#' # NOTE: It is easiest to always set saveTable to TRUE. If it is set to FALSE, there is a chance the
#' legend will be cut off, and the user will have to adjust the proportions in the Plots pane for each
#' table before group_block_anno is run. There is a pause in execution added to the code right before
#' the subtitle is added in order to give the user time to do this.
#' @param saveFolder file path containing folder the finished tables will be saved in
#' @param all_width width of cells for heatmap showing combinations for all mutations
#' @param all_height height of cells for heatmap showing combinations for all mutations
#' @param hr_width width of cells for heatmap showing combinations for high-risk mutations
#' @param hr_height height of cells for heatmap showing combinations for high-risk mutations
#' @param counts_font_size font size for counts at the bottom of the bar chart
#' @import grDevices
#' @export
#'
makeCoOccHeatmaps <- function(DB, startDate = NULL, endDate = NULL, archiveDBFile = NULL, saveTable = FALSE, saveFolder = NULL, all_width = NULL, all_height = NULL, hr_width = NULL, hr_height = NULL, counts_font_size = NULL) {
  n <- NULL

  setStrings2AssayDB()

  # Get total num sequences:
  con <- DBI::dbConnect(RSQLite::SQLite(), DB)
  totalSeqs <- dplyr::tbl(con, "Seq_Info") %>%
    dplyr::count() %>%
    dplyr::collect() %>%
    dplyr::pull(n)
  DBI::dbDisconnect(con)

  # make dataframes to be used in heatmaps:
  dfAll <- makeDFcoOccHeatmap(DB)
  dfAll3Assay <- makeDFcoOccHeatmap(DB, threeAssay = TRUE)
  dfAll3AssayPast <- makeDFcoOccHeatmap(DB, threeAssay = TRUE, startDate = startDate, endDate = endDate, archiveDBFile = archiveDBFile)

  dfHR <- makeDFcoOccHeatmap(DB, HR = TRUE)
  dfHR3Assay <- makeDFcoOccHeatmap(DB, HR = TRUE, threeAssay = TRUE)
  dfHR3AssayPast <- makeDFcoOccHeatmap(DB, HR = TRUE, threeAssay = TRUE, startDate = startDate, endDate = endDate, archiveDBFile = archiveDBFile)


  if (!dir.exists(saveFolder)) {
    dir.create(saveFolder, recursive = TRUE)
  }

  # -----------------------------------------------------------------------------
  # 3-assay .1% test (for conclusion slide in SPC presentations) ---------------

  # write 3-assay dataframes to csv file:
  readr::write_csv(
    dfAll3Assay,
    paste0(saveFolder, "3_assay_heatmaps_", startDate, "_to_", endDate, ".csv")
  )

  readr::write_csv(dfAll3AssayPast,
    paste0(saveFolder, "3_assay_heatmaps_", startDate, "_to_", endDate, ".csv"),
    append = TRUE
  )

  readr::write_csv(dfHR3Assay,
    paste0(saveFolder, "3_assay_heatmaps_", startDate, "_to_", endDate, ".csv"),
    append = TRUE
  )

  readr::write_csv(dfHR3AssayPast,
    paste0(saveFolder, "3_assay_heatmaps_", startDate, "_to_", endDate, ".csv"),
    append = TRUE
  )

  # print(dfAll3Assay)
  # print(dfAll3AssayPast)
  # print(dfHR3Assay)
  # print(dfHR3AssayPast)


  # calculate change in proportion of mutations affecting all three assays and write to file:
  allCurProp <- getThreeAssayProportion(dfAll3Assay)
  allPastProp <- getThreeAssayProportion(dfAll3AssayPast)
  HRcurProp <- getThreeAssayProportion(dfHR3Assay)
  HRpastProp <- getThreeAssayProportion(dfHR3AssayPast)

  allPropChange <- allCurProp - allPastProp
  HRpropChange <- HRcurProp - HRpastProp

  print(paste0("Change in proportion of seqs with mutations affecting a,d,e assays, all mutations: ", allPropChange))
  print(paste0("Change in proportion of seqs with mutations affecting a,d,e assays, HR mutations: ", HRpropChange))

  prop_df <- data.frame(all_prop_change = allPropChange, hr_prop_change = HRpropChange)
  readr::write_csv(
    prop_df,
    paste0(saveFolder, "3_assay_prop_change_", startDate, "_to_", endDate, ".csv")
  )

  # ------------------------------------------------------------------------------


  # Make heatmaps --------------------------------------------------------------

  # All mutations ------------------------------

  # Remove Assays_Affected and Frequency rows from dataframe:

  assays_affected <- dfAll["Assays_Affected", ]

  # assays_affected is used to split the columns into groups
  assays_affected <- as.numeric(assays_affected[1, ])

  # col_groups is used to create labels for column groups (# assays affected)
  col_groups <- unique(assays_affected)

  remove_row <- c("Assays_Affected")
  dfAll <- dfAll[!(row.names(dfAll) %in% remove_row), ]

  counts <- utils::tail(dfAll, n = 1)
  counts <- as.numeric(counts[1, ])

  counts3Assay <- utils::tail(dfAll3Assay, n = 1)
  counts3Assay <- as.numeric(counts3Assay[1, ])

  # Add commas for heatmap display
  prettyCounts <- format(counts, big.mark = ",", trim = TRUE)
  prettyTotal <- format(totalSeqs, big.mark = ",", trim = TRUE)

  percents <- round(((counts / totalSeqs) * 100), 1)

  dfAll <- utils::head(dfAll, -1) # Remove Freq row

  matrixAll <- data.matrix(dfAll, rownames.force = TRUE) # Convert to matrix for heatmap creation



  heatmapAll <- formatHeatmap(matrixAll, col_groups, assays_affected, counts, prettyCounts, percents, prettyTotal, counts_font_size, HR = FALSE)


  if (saveTable) {
    grDevices::png(file = paste0(saveFolder, "heatmap_all.png"), width = all_width, height = all_height, units = "px", res = 100)

    ComplexHeatmap::draw(heatmapAll, padding = grid::unit(c(20, 2, 2, 8), "mm")) # bottom, left, top, right paddings
  } else {
    ComplexHeatmap::draw(heatmapAll, padding = grid::unit(c(10, 2, 2, 2), "mm"))
  }

  # Total sequences text box at bottom of heatmap
  ComplexHeatmap::decorate_heatmap_body("hm1", {
    grid::grid.text(
      label = paste0("Total sequences: ", prettyTotal),
      x = grid::unit(10, "cm"), # 345 22
      y = grid::unit(0, "cm"), # -50 -20
      gp = grid::gpar(fontsize = 11)
    )
  })

  print(heatmapAll)

  # if (!saveTable) { # Give the user some time to adjust proportions of table in Plots window (group_block_anno calculates position of subtitle based on those proportions)
  #   Sys.sleep(20)
  # }

  # assays affected header
  group_block_anno(1:length(col_groups), "label",
    gp = grid::gpar(col = "white"),
    label_gp = grid::gpar(fontsize = 14),
    label = "# Assays Affected"
  )

  if (saveTable) {
    grDevices::dev.off()
  }




  # 3' end mutations --------------------------

  # Remove Assays_Affected and Frequency rows from dataframe:

  assays_affected <- dfHR["Assays_Affected", ]

  # assays_affected is used to split the columns into groups
  assays_affected <- as.numeric(assays_affected[1, ])

  # col_groups is used to create labels for column groups (# assays affected)
  col_groups <- unique(assays_affected)

  remove_row <- c("Assays_Affected")
  dfHR <- dfHR[!(row.names(dfHR) %in% remove_row), ]

  counts <- utils::tail(dfHR, n = 1)
  counts <- as.numeric(counts[1, ])

  # Add commas for heatmap display
  prettyCounts <- format(counts, big.mark = ",", trim = TRUE)
  prettyTotal <- format(totalSeqs, big.mark = ",", trim = TRUE)

  percents <- round(((counts / totalSeqs) * 100), 1)

  dfHR <- utils::head(dfHR, -1) # Remove Freq row

  matrixHR <- data.matrix(dfHR, rownames.force = TRUE) # Convert to matrix for heatmap creation



  heatmapHR <- formatHeatmap(matrixHR, col_groups, assays_affected, counts, prettyCounts, percents, prettyTotal, counts_font_size, HR = TRUE)



  if (saveTable) {
    grDevices::png(file = paste0(saveFolder, "heatmap_hr.png"), width = hr_width, height = hr_height, units = "px", res = 100)

    ComplexHeatmap::draw(heatmapHR, padding = grid::unit(c(10, 2, 2, 8), "mm")) # bottom, left, top, right paddings)
  } else {
    ComplexHeatmap::draw(heatmapHR, padding = grid::unit(c(10, 2, 2, 2), "mm"))
  }


  # Total sequences textbox at bottom of heatmap
  ComplexHeatmap::decorate_heatmap_body("hm2", {
    grid::grid.text(
      label = paste0("Total sequences: ", prettyTotal),
      x = grid::unit(22, "mm"), # 345
      y = grid::unit(-20, "mm"), # -50
      gp = grid::gpar(fontsize = 11)
    )
  })

  print(heatmapHR)
  print("total sequences:")
  print(prettyTotal)

  # if (!saveTable) { # Give the user some time to adjust proportions of table in Plots window (group_block_anno calculates position of subtitle based on those proportions)
  #   Sys.sleep(20)
  # }

  # assays affected header
  group_block_anno(1:length(col_groups), "label",
    gp = grid::gpar(col = "white"),
    label_gp = grid::gpar(fontsize = 14),
    label = "# Assays Affected"
  )

  if (saveTable) {
    grDevices::dev.off()
  }

}
