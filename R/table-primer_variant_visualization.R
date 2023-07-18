###### Make Primer Variant Visualization Table for COVID Inclusivity Reports ######

#' Reorder vect - change the position of what in vect to where
#'
#' @param vect vector containing elements to reorder
#' @param what element of vector to move
#' @param where new position of what
#'
#' @return vect
#' @noRd
#'
reorderVect <- function(vect, what, where) {
  idx <- which(vect == what)
  vect <- append(vect[-idx], vect[idx], where - 1)
  return(vect)
}



#' Make Primer Variant Mutation Dataframe
#'
#' @param start_date Start date for subset db
#' @param end_date End date for subset db
#' @param varIds Vector of variant IDs to be used in the visualization, provided if subsetDBFile is NULL (not using Sig_Var_IDs from the subsetDBFile)
#' @param archiveDBFile File path to archive db (used to retrieve older subset of data for comparison; older subset uses same length of time as input subset db, e.g. start_date is "2022-5" and end_date is "2022-07", then older subset start date is "2022-02 and "2022-04)
#' @param subsetDBFile File path to subset db, provided if using variant IDs from Sig_Var_IDs as those in the visualization (varIds parameter should be NULL if this is provided)
#' @param printableDoubles Boolean, default is FALSE, if TRUE will show double mutations for a single variant ID in a single column of the visualization
#'
#' @return unique_var_table, a table with mutations and trends for each significant variant, comparing this subset to an older subset
#' @noRd
#'
makePrimerVarMutationDF <- function(start_date, end_date,
                                    varIds = NULL,
                                    archiveDBFile,
                                    subsetDBFile = NULL,
                                    printableDoubles = FALSE) {
  Trend <- `Mutation Count` <- Pos_from_3P <- Last_Prop <- Curr_Prop <- Curr_Q_Freq <-
    Last_Q_Freq <- Primer_ID <- Primer <- Primer_Dir <- Primer_Pos <- Direction <-
    Reaction <- Seq_Change <- Mismatch_Type <- Var_ID <- n <- Accession_ID <-
    Collection_Date <- Assay_Name <- Indel_Length <- NULL

  if (!is.null(subsetDBFile)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), subsetDBFile)
    varIds <- dplyr::tbl(con, "Sig_Var_IDs") %>%
      dplyr::collect() %>%
      dplyr::pull()
    DBI::dbDisconnect(con)
  } else {
    if (is.null(varIds)) {
      warning("Neither a subset DB or Var List was provided. Please check your inputs.")
      return()
    }
  }

  # previous time period does not overlap current
  monthDiff <- lubridate::interval(lubridate::ym(start_date), lubridate::ym(end_date)) %/% months(1)
  last_end_date <- lubridate::ym(start_date) - months(1)
  last_start_date <- lubridate::ym(start_date) - months(1 + monthDiff)
  last_end_date <- toString(format(last_end_date, "%Y-%m"))
  last_start_date <- toString(format(last_start_date, "%Y-%m"))

  con <- DBI::dbConnect(RSQLite::SQLite(), archiveDBFile)
  primers <- dplyr::tbl(con, "Primer_Info") %>% dplyr::collect()
  idMap <- dplyr::tbl(con, "ID_Map") %>% dplyr::collect()
  varType <- dplyr::tbl(con, "Variant_Type") %>% dplyr::collect()
  seqInfo <- dplyr::tbl(con, "Seq_Info") %>% dplyr::collect()
  varInfo <- dplyr::tbl(con, "Variant_Info") %>% dplyr::collect()
  DBI::dbDisconnect(con)

  # get primer counts by assay
  assays <- primers %>%
    dplyr::group_by(Assay_Name) %>%
    dplyr::tally()

  # get all accessions between start and end date
  currentAccs <- seqInfo %>%
    dplyr::filter(Collection_Date <= end_date & Collection_Date >= start_date) %>%
    dplyr::select(Accession_ID)

  # get number of current accessions
  totalSeqs <- currentAccs %>%
    dplyr::count() %>%
    dplyr::pull(n)
  currTotSeqs <- totalSeqs

  # get var ID for each accession ID, filtering for only var IDs in varIDs vector;
  # then get number of accessions associated with each var ID (counts will total more than number of accessions because an accession can have multiple var IDs)
  var_counts <- currentAccs %>%
    dplyr::inner_join(idMap, by = "Accession_ID") %>%
    dplyr::filter(Var_ID %in% varIds) %>%
    dplyr::group_by(Var_ID) %>%
    dplyr::summarise(Curr_Q_Freq = dplyr::n_distinct(Accession_ID))

  #-------------------------------------------
  # Previous time period
  pastAccs <- seqInfo %>%
    dplyr::filter(Collection_Date <= last_end_date & Collection_Date >= last_start_date) %>%
    dplyr::select(Accession_ID)
  pastTotalSeqs <- pastAccs %>%
    dplyr::count() %>%
    dplyr::pull(n)

  # get var ID for each accession ID, filtering for only var IDs in varIDs vector; then get number of accessions associated with each var ID (counts will total more than number of accessions because an accession can have multiple var IDs)
  var_counts_Q_last <- pastAccs %>%
    dplyr::inner_join(idMap, by = "Accession_ID") %>%
    dplyr::filter(Var_ID %in% varIds) %>%
    dplyr::group_by(Var_ID) %>%
    dplyr::summarise(Last_Q_Freq = dplyr::n_distinct(Accession_ID))
  #-------------------------------------------

  # Join both time periods
  full_freq_table <- var_counts %>% dplyr::left_join(var_counts_Q_last, by = "Var_ID")

  # Convert NA frequency values to 0
  full_freq_table <- full_freq_table %>% dplyr::mutate_at(c("Curr_Q_Freq", "Last_Q_Freq"), ~ tidyr::replace_na(., 0))

  # Get mutation information associated with each var ID
  mutData <- varType %>% dplyr::filter(Var_ID %in% varIds)

  # If indel, make Mismatch_Type column value match Seq_Change column value
  mutData <- mutData %>%
    dplyr::mutate(Mismatch_Type = ifelse(Mismatch_Type == "", Seq_Change, Mismatch_Type)) %>%
    dplyr::distinct()

  # Add Primer column which combines direction and reaction into the standard 2 letters, select columns
  primerDir <- primers %>%
    dplyr::mutate(Primer_Pos = substr(Reaction, 1, 1)) %>%
    dplyr::mutate(Primer_Dir = substr(Direction, 1, 1)) %>%
    tidyr::unite("Primer", Primer_Pos:Primer_Dir, remove = TRUE, sep = "") %>%
    dplyr::select(Assay_Name, Primer, Primer_ID)

  # Combine full_freq_table (last and current accession frequencies for each variant ID) with
  # mutData table (mutation info for each variant ID) and with primerDir (assay name and 2-letter
  # direction/reaction for each primer ID); also add Last_Prop and Curr_Prop, which are
  # accession frequencies / total num accessions for each variant ID
  unique_var_table <- full_freq_table %>%
    dplyr::left_join(mutData, by = "Var_ID") %>%
    dplyr::left_join(primerDir, by = "Primer_ID") %>%
    dplyr::mutate(Last_Prop = (Last_Q_Freq / pastTotalSeqs), Curr_Prop = (Curr_Q_Freq / currTotSeqs))

  # Get mutation count for each variant ID (a mutation count > 1 means there was more than 1 mismatch/indel associated with the same variant ID ie on the same primer)
  mutCounts <- unique_var_table %>%
    dplyr::distinct() %>%
    dplyr::group_by(Var_ID) %>%
    dplyr::summarise(`Mutation Count` = dplyr::n())

  # Add Trend column, which compares Last_Prop with Curr_Prop--"p" if at least .1% increase,
  # "q" if at least .1% decrease, and "D" if no change or not big enough change; will be put in Wingdings 3 in Excel
  unique_var_table <- unique_var_table %>%
    dplyr::left_join(mutCounts, by = "Var_ID") %>%
    dplyr::mutate(Trend = ifelse((Curr_Prop - Last_Prop) > .001, "p",
      ifelse((Curr_Prop - Last_Prop) < -.001, "q", "D")
    )) %>%
    dplyr::select(Assay_Name, Var_ID, Primer, Pos_from_3P, Mismatch_Type, Indel_Length, `Mutation Count`, Last_Q_Freq, Last_Prop, Curr_Q_Freq, Curr_Prop, Trend)

  # If printableDoubles is True, combine rows for mutations associated with same variant IDs, separating values in Pos_from_3P and Mismatch_Type with a comma
  if (printableDoubles) {
    if (sum(unique_var_table$`Mutation Count`) > nrow(unique_var_table)) { # Check if there are multiple mutations for any variant IDs
      multMutsDF <- unique_var_table %>%
        dplyr::filter(`Mutation Count` > 1) %>%
        dplyr::group_by(Var_ID) %>%
        tidyr::nest(Pos_from_3P = Pos_from_3P, Mismatch_Type = Mismatch_Type) %>%
        dplyr::rowwise() %>%
        dplyr::mutate_if(is.list, ~ paste(unlist(.), collapse = ","))

      unique_var_table <- unique_var_table %>%
        dplyr::filter(`Mutation Count` <= 1) %>%
        rbind(multMutsDF) # add multiple mutations df to single mutations df
    }
  }

  # Order dataframe by assay name first, and descending current frequency second
  unique_var_table <- unique_var_table %>% dplyr::arrange(Assay_Name, dplyr::desc(Curr_Q_Freq))

  return(unique_var_table)
}




#' Make primer variant frequency and lineage dataframes
#'
#' @param dbFile File path for subset db
#'
#' @return freqDF, which has frequencies for each significant var ID, and linDF, which has lineage information for each significant var ID (a single var ID can be marked as part of multiple lineages in the database)
#' @noRd
#'
makePrimerVarFreqLinDFs <- function(dbFile) {
  Percent_Total <- Var_Count <- Var_by_Lin <- Var_ID <- Lineage_Count <- Subgroup <-
    WHO_Designation <- Primer_ID <- is_Masked <- Assay <- Primer <- Primer_Dir <-
    Primer_Pos <- Direction <- Reaction <- Displayed_Lineage <- Displayed_Group <-
    Assay_Name <- Pos_from_3P <- NULL

  con <- DBI::dbConnect(RSQLite::SQLite(), dbFile)
  # Tables with Needed Information
  idMap <- dplyr::tbl(con, "ID_Map")
  seqInfo <- dplyr::tbl(con, "Seq_Info")
  varInfo <- dplyr::tbl(con, "Variant_Info")
  sigVars <- dplyr::tbl(con, "Sig_Var_Ids") %>% dplyr::collect()
  primers <- dplyr::tbl(con, "Primer_Info") %>% dplyr::collect()
  varType <- dplyr::tbl(con, "Variant_Type") %>% dplyr::collect()

  # Get Tables
  seqInfo <- seqInfo %>%
    dplyr::select(s$ACC_ID, s$GROUP, s$SUBGROUP, Displayed_Lineage, Displayed_Group, s$LOC_ACC) %>%
    dplyr::collect()
  idMap <- idMap %>% dplyr::collect()
  varCounts <- varInfo %>%
    dplyr::select(s$VAR_ID, s$SEQ_VAR_COUNT) %>%
    dplyr::collect()
  var2Primer <- varInfo %>%
    dplyr::select(s$VAR_ID, s$PRIMER_ID) %>%
    dplyr::collect()

  DBI::dbDisconnect(con)




  # Get keep and remove lists---------------------------------------------------

  sig_var_IDs <- sigVars %>%
    dplyr::select(Var_ID) %>%
    unique() %>%
    dplyr::pull()

  # Get all significant var ID's where E IR is affected (for keep list)
  e_ir_primer_IDs <- primers %>%
    dplyr::filter(Assay_Name == "SARS-CoV-2e" &
      Direction == "Reverse" &
      Reaction == "Inner") %>%
    dplyr::select(Primer_ID) %>%
    dplyr::pull()

  e_ir_sig_var_IDs <- var2Primer %>%
    dplyr::filter(Primer_ID %in% e_ir_primer_IDs &
      Var_ID %in% sig_var_IDs) %>%
    dplyr::select(Var_ID) %>%
    dplyr::pull() %>%
    unique()

  print("for keep list")
  print(e_ir_sig_var_IDs)

  # Get all significant var ID's where F IR is affected (for remove list)
  f_ir_primer_IDs <- primers %>%
    dplyr::filter(Assay_Name == "SARS-CoV-2f" &
      Direction == "Reverse" &
      Reaction == "Inner") %>%
    dplyr::select(Primer_ID) %>%
    dplyr::pull()

  f_ir_sig_var_IDs <- var2Primer %>%
    dplyr::filter(Primer_ID %in% f_ir_primer_IDs &
      Var_ID %in% sig_var_IDs) %>%
    dplyr::select(Var_ID) %>%
    dplyr::pull() %>%
    unique()

  print("for remove list")
  print(f_ir_sig_var_IDs)

  # There are some locations on the left side of E and the right side of F that
  # don't overlap, we want to remove those from lists---------------------------
  # Primer ID for E IR: 16
  # Primer ID for F IR: 21

  # Get significant var ID's from right side of the F IR that doesn't overlap with E
  e_ir_var_IDs <- varType %>%
    dplyr::filter(Primer_ID == 16) %>%
    dplyr::select(Var_ID) %>%
    dplyr::pull()
  no_overlap_f_ir_var_IDs <- varType %>%
    dplyr::filter(Pos_from_3P < 3 &
      Var_ID %in% e_ir_var_IDs &
      Var_ID %in% sig_var_IDs) %>%
    dplyr::select(Var_ID) %>%
    dplyr::pull()

  # Remove those ID's from keep list
  keep_list <- e_ir_sig_var_IDs[!e_ir_sig_var_IDs %in% no_overlap_f_ir_var_IDs]

  print("final keep list")
  print(keep_list)

  # Get significant var ID's from left side of E IR that doesn't overlap with F
  f_ir_var_IDs <- varType %>%
    dplyr::filter(Primer_ID == 21) %>%
    dplyr::select(Var_ID) %>%
    dplyr::pull()
  no_overlap_e_ir_var_IDs <- varType %>%
    dplyr::filter(Pos_from_3P > 16 &
      Var_ID %in% f_ir_var_IDs &
      Var_ID %in% sig_var_IDs) %>%
    dplyr::select(Var_ID) %>%
    dplyr::pull()

  # Remove those ID's from remove list
  remove_list <- f_ir_sig_var_IDs[!f_ir_sig_var_IDs %in% no_overlap_e_ir_var_IDs]

  print("final remove list")
  print(remove_list)


  # Prep for same assay filtering-------------------------------------------
  # Map each primer ID to its direction and reaction (e.g. OF)
  pNames <- primers %>%
    dplyr::mutate(Primer_Pos = substr(Reaction, 1, 1)) %>%
    dplyr::mutate(Primer_Dir = substr(Direction, 1, 1)) %>%
    tidyr::unite("Primer", Primer_Pos:Primer_Dir, remove = FALSE, sep = "") %>%
    dplyr::select(s$PRIMER_ID, Primer)

  # Map each assay ID to its assay letter name
  aNames <- primers %>%
    dplyr::mutate(Assay = toupper(substr(get(s$ASSAY_NAME), (nchar(get(s$ASSAY_NAME))), nchar(get(s$ASSAY_NAME))))) %>%
    dplyr::select(s$ASSAY_ID, Assay) %>%
    dplyr::distinct()

  primers <- primers %>%
    dplyr::inner_join(pNames, by = "Primer_ID") %>%
    dplyr::inner_join(aNames, by = "Assay_ID") %>%
    dplyr::select(Assay, s$ASSAY_ID, Primer, s$PRIMER_ID, s$IS_COSPOT_STRING, is_Masked)

  # OF, IF, and OR primers are the same sequences in C and D assays and in E and F assays
  # However, they have different primer ID's and the nucleotide positions are slightly different
  # To account for this overlap, we keep only the D and E data where the primers overlap in this visualization

  # first, get ID's for primers that have same sequences
  c_same <- primers %>%
    dplyr::filter(Assay == "C") %>%
    dplyr::filter(Primer != "IR") %>%
    dplyr::pull(get(s$PRIMER_ID))
  d_same <- primers %>%
    dplyr::filter(Assay == "D") %>%
    dplyr::filter(Primer != "IR") %>%
    dplyr::pull(get(s$PRIMER_ID))
  e_same <- primers %>%
    dplyr::filter(Assay == "E") %>%
    dplyr::filter(Primer != "IR") %>%
    dplyr::pull(get(s$PRIMER_ID))
  f_same <- primers %>%
    dplyr::filter(Assay == "F") %>%
    dplyr::filter(Primer != "IR") %>%
    dplyr::pull(get(s$PRIMER_ID))

  equal_primers <- c(c_same, d_same, e_same, f_same)

  # Get all significant Var ID's
  all_ids <- sigVars %>% dplyr::pull(get(s$VAR_ID))

  # Get Var ID's that are significant and that match with primer ID's in the equal primers list
  equal_Var_IDs <- sigVars %>%
    dplyr::inner_join(var2Primer, by = "Var_ID") %>%
    dplyr::filter(get(s$PRIMER_ID) %in% equal_primers) %>%
    dplyr::pull(get(s$VAR_ID))

  # All other Var ID's
  Non_equal_Var_IDs <- base::setdiff(all_ids, equal_Var_IDs)

  seqInfoCount <- seqInfo %>% dplyr::summarise(n = dplyr::n())
  totalSeqs <- seqInfoCount$n[1]


  LinSeqCounts <- seqInfo %>%
    dplyr::group_by(WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group) %>%
    dplyr::summarise(Lineage_Count = dplyr::n()) %>%
    dplyr::mutate(Per_Var_by_Lin = round(Lineage_Count / totalSeqs * 100, 2)) %>%
    dplyr::mutate(Var_ID = "All Seqs") %>%
    dplyr::mutate(Var_Count = totalSeqs) %>%
    dplyr::mutate(Assay = "None")

  colnames(LinSeqCounts) <- c("WHO_Designation", "Subgroup", "Displayed_Lineage", "Displayed_Group", "Var_by_Lin", "Per_Var_by_Lin", "Var_ID", "Var_Count", "Assay")

  if (totalSeqs != sum(LinSeqCounts$Var_by_Lin)) {
    print("Error: All sequences not in the analysis")
  }

  filtered <- sigVars %>% dplyr::inner_join(idMap, by = "Var_ID")


  LinCount <- filtered %>%
    dplyr::group_by(Var_ID, WHO_Designation, Subgroup, Displayed_Lineage, Displayed_Group) %>%
    dplyr::summarise(Var_by_Lin = dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(Var_by_Lin))

  # Var_by_Lin is new name for lineage count in LinSeqCounts table, which is count for each WHO_Designation-Subgroup combo
  # Var_by_Lin is Var_ID-WHO_Designation-Subgroup combo count in LinCount table
  # Per_Var_by_Lin is Var_by_Lin count / Var_Count (which is count for just the variant ID)
  LinCount <- LinCount %>%
    dplyr::left_join(varCounts, by = "Var_ID") %>%
    dplyr::mutate(Per_Var_by_Lin = round(Var_by_Lin / Var_Count * 100, 2)) %>%
    dplyr::arrange(Var_Count)

  LinCount$Var_ID <- as.character(LinCount[[s$VAR_ID]])

  var2Primer$Var_ID <- as.character(var2Primer[[s$VAR_ID]])

  # Map from primer ID to its assay letter
  primerId2Assay <- primers %>% dplyr::select(s$PRIMER_ID, Assay)

  # Mark shared sequences in C/D and E/F assays (and the assignment was D or E) with both letters; mark C and F-assigned assays to be removed
  lin_count_by_Var_Same <- LinCount %>%
    dplyr::inner_join(var2Primer, by = "Var_ID") %>%
    dplyr::inner_join(primerId2Assay, by = "Primer_ID") %>%
    dplyr::mutate(Assay = dplyr::case_when(
      Primer_ID %in% d_same ~ "C/D",
      Primer_ID %in% e_same ~ "E/F",
      Primer_ID %in% c_same | Primer_ID %in% f_same ~ "Rm",
      TRUE ~ Assay
    )) %>%
    dplyr::select(-Primer_ID) %>%
    dplyr::ungroup()

  # Mark sequences in keep list (for E and F IR primer overlap) as E/F; mark sequences in remove list to be removed
  lin_graph <- lin_count_by_Var_Same %>% dplyr::mutate(Assay = dplyr::case_when(
    Var_ID %in% keep_list ~ "E/F",
    Var_ID %in% remove_list ~ "Rm",
    TRUE ~ Assay
  ))

  lin_graph <- lin_graph %>% dplyr::filter(Assay != "Rm")


  # Add LinSeqCounts information, which contains total var count and counts for each WHO-Designation-Subgroup mapping, without selecting Var ID's
  lin_graph <- lin_graph %>%
    dplyr::full_join(LinSeqCounts,
      by = c(
        "Var_ID",
        "WHO_Designation",
        "Subgroup",
        "Displayed_Lineage",
        "Displayed_Group",
        "Var_by_Lin",
        "Per_Var_by_Lin",
        "Var_Count",
        "Assay"
      )
    ) %>%
    dplyr::arrange(Assay, dplyr::desc(s$SEQ_VAR_COUNT))


  var2Primer[[s$VAR_ID]] <- as.numeric(var2Primer[[s$VAR_ID]])

  # Get counts for Var ID's that are significant and occur in lin_graph
  sig_var_counts <- sigVars %>%
    dplyr::inner_join(varCounts, by = "Var_ID") %>%
    dplyr::filter(Var_ID %in% unique(lin_graph$Var_ID))

  # Map from Assay letter to Var ID
  var2Assay <- lin_graph %>%
    dplyr::select(Assay, s$VAR_ID) %>%
    dplyr::distinct()

  sig_var_counts$Var_ID <- as.character(sig_var_counts$Var_ID)

  # Add assay letter to sig var counts
  sig_var_counts <- sig_var_counts %>%
    dplyr::inner_join(var2Assay, by = "Var_ID") %>%
    rbind(c("All Seqs", totalSeqs, "None"))

  sig_var_counts$Var_Count <- as.numeric(sig_var_counts$Var_Count)

  # Add percentage of total sequences for every significant variant row
  sig_var_counts <- sig_var_counts %>% dplyr::mutate(Percent_Total = round(Var_Count / totalSeqs * 100, 2))

  return(list(freqDF = sig_var_counts, linDF = lin_graph))
}




#' Make primer variant visualization
#'
#' @description Create a figure with three parts: a bar graph with variant (set of mutations at a primer binding region) counts organized by assay on top, a table with specific variant information in the middle, and a bar graph with lineage makeup per variant on the bottom
#' @param save_table Whether or not to save the table
#' @param save_folder Folder name for saving table
#' @param start_date Start date for subset db
#' @param end_date End date for subset db
#' @param archive_file File path to archive db (used to retrieve older subset of data for comparison; older subset uses same length of time as input subset db, e.g. start_date is "2022-5" and end_date is "2022-07", then older subset start date is "2022-02 and "2022-04)
#' @param subset_file File path to subset db
#' @param who_var_file File path to WHO variant file (for getting Displayed Lineage order for legend in bottom stacked plot)
#' @param color_palette Letter value A - H corresponding to a viridis color palette; palettes are [A: magma, B: inferno, C: plasma, D: viridis, E: cividis, F: rocket, G: mako, H: turbo]; if color value doesn't match, it defaults to "viridis" palette
#' @param section_heights Optional vector containing the heights of the three sections of the visualization (top graph, middle table, and bottom graph). If no value is supplied, default is c(2,1.25,2).
#' @param figure_res Resolution of saved png figure
#' @param figure_width Width of saved png figure
#' @param figure_height Height of saved png figure
#'
#' @export
#'
makePrimerVarTable <- function(save_table, save_folder, start_date, end_date, archive_file, subset_file, who_var_file, color_palette, section_heights = NULL, figure_res = NULL, figure_width = NULL, figure_height = NULL) {
  setStrings2AssayDB()

  Percent_Total <- Per_Var_by_Lin <- Subgroup <- n <- Trend <- Curr_Prop <- Last_Prop <-
    Last_Q_Freq <- `Mutation Count` <- Var_ID <- Assay <- Assay_Name <- Curr_Q_Freq <-
    Mismatch_Type <- Pos_from_3P <- Primer <- lastSubDate <- Indel_Length <- Displayed_Lineage <-
    Displayed_Lineage_Order <- NULL

  varPrimerTable <- makePrimerVarMutationDF(start_date, end_date,
    archiveDBFile = archive_file,
    subsetDBFile = subset_file,
    printableDoubles = TRUE
  )

  varPrimerTableSave <- varPrimerTable %>%
    dplyr::group_by(Primer, Pos_from_3P, Mismatch_Type, Curr_Q_Freq) %>%
    dplyr::mutate(Assay = substr(Assay_Name, (nchar(Assay_Name) - 1), nchar(Assay_Name))) %>%
    dplyr::select(-Assay_Name) %>%
    tidyr::nest(Assay = Assay, Var_ID = Var_ID) %>%
    dplyr::rowwise() %>%
    dplyr::mutate_if(is.list, ~ paste(unlist(.), collapse = ", ")) %>%
    dplyr::select(Assay, Var_ID, Primer, Pos_from_3P, Mismatch_Type, Indel_Length, `Mutation Count`, Last_Q_Freq, Last_Prop, Curr_Q_Freq, Curr_Prop, Trend)

  if (save_table) {
    if (!dir.exists(save_folder)) {
      dir.create(save_folder, recursive = TRUE)
    }

    readr::write_csv(
      varPrimerTableSave,
      paste0(save_folder, "unique_vars_", start_date, "_to_", end_date, "_all.csv")
    )
  }

  print(varPrimerTableSave %>% dplyr::filter(Primer == "IR") %>% dplyr::filter(Assay == "2e" | Assay == "2f") %>% dplyr::arrange(as.integer(Var_ID)))

  graphDFs <- makePrimerVarFreqLinDFs(subset_file)

  # BUILD GRAPHS-----------------------------------------------------------------------

  lin_graph <- graphDFs$linDF
  sig_var_counts <- graphDFs$freqDF

  # Build ordering levels
  vars <- unique(lin_graph$Var_ID)
  var_levs <- c(vars[1:length(vars) - 1], "All Seqs")

  BAR_WIDTH <- 0.75

  assayNames <- sort(unique(sig_var_counts$Assay))
  assayNames[which(assayNames == "None")] <- "All"

  x_line_pos <- c()
  for (assay in assayNames) {
    count <- sig_var_counts %>%
      dplyr::filter(Assay == assay) %>%
      dplyr::count() %>%
      dplyr::pull(n)
    if (length(x_line_pos) == 0) {
      last <- count + .5
      x_line_pos <- c(last)
    } else {
      pos <- last + count
      last <- pos
      x_line_pos <- c(x_line_pos, pos)
    }
  }
  assayPlaces <- c()
  last <- .5
  for (pos in x_line_pos) {
    new <- ((pos - last) / 2) + last
    assayPlaces <- c(assayPlaces, new)
    last <- pos
  }

  assayPlaces[length(assayPlaces)] <- assayPlaces[length(assayPlaces)] + .5


  # Pull lineage order from who_var_file
  legend_lineage_order <- readr::read_csv(who_var_file, show_col_types = FALSE) %>%
    dplyr::select(Displayed_Lineage_Order) %>%
    dplyr::pull()
  legend_lineage_order <- legend_lineage_order[!is.na(legend_lineage_order)]

  # Make Displayed Lineage a factor so that order of legend values in bottom figure matches legend_lineage_order
  lin_graph$Displayed_Lineage <- factor(lin_graph$Displayed_Lineage, levels=legend_lineage_order)

  bottom_stacked <- ggplot2::ggplot(lin_graph, ggplot2::aes(fill = Displayed_Lineage, y = Per_Var_by_Lin, x = factor(Var_ID, levels = var_levs), width = BAR_WIDTH)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::geom_vline(
      xintercept = x_line_pos,
      colour = "#1520A6",
      alpha = 0.5
    ) +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    viridis::scale_fill_viridis(name="Lineage", discrete=TRUE, option=color_palette)


  top_counts <- ggplot2::ggplot(sig_var_counts, ggplot2::aes(y = Percent_Total, x = factor(Var_ID, levels = var_levs), width = BAR_WIDTH)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_bw() +
    ggplot2::ylab("") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line.y.left = ggplot2::element_line(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_vline(
      xintercept = x_line_pos,
      colour = "#1520A6",
      alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = c(100),
      colour = "#1520A6",
      alpha = 0.5
    ) +
    ggplot2::annotate("text",
      x = assayPlaces,
      y = 107,
      label = assayNames,
      colour = "grey10",
      size = 4
    )



  # Middle table----------------------------------
  table_data <- varPrimerTable %>%
    dplyr::select(Var_ID, Curr_Prop, Trend, Primer, Pos_from_3P, Mismatch_Type, Indel_Length) %>%
    dplyr::filter(Var_ID %in% sig_var_counts$Var_ID) %>%
    dplyr::mutate(Curr_Prop = round(Curr_Prop * 100, 1)) %>%
    dplyr::mutate(Trend = dplyr::case_when(
      Trend == "p" ~ "\U25B2",
      Trend == "D" ~ "\U21C6",
      Trend == "q" ~ "\U25BC"
    ))

  # Indels
  for (i in seq_along(table_data$Mismatch_Type)) {
    misType <- table_data$Mismatch_Type[i]
    if (startsWith(misType, "Deletion")) {
      table_data$Mismatch_Type[i] <- paste0("DL", table_data$Indel_Length[i])
    } else if (startsWith(misType, "Insertion")) {
      table_data$Mismatch_Type[i] <- paste0("IN", table_data$Indel_Length[i])
    }
  }

  # We don't need Indel_Length anymore
  table_data <- table_data %>% dplyr::select(-Indel_Length)

  table_data <- rbind(table_data, c("All Seqs", 100, NA, NA, NA, NA))

  print("All primer variants")
  print(table_data)

  names <- table_data$Var_ID
  table_data.T <- as.data.frame(as.matrix(t(table_data[, -1])))
  colnames(table_data.T) <- names

  ordered_table <- table_data.T %>%
    dplyr::select(dplyr::all_of(var_levs)) %>%
    dplyr::rename("All" = "All Seqs")

  p2 <- gridExtra::tableGrob(ordered_table, rows = NULL, cols = NULL)
  # Set widths/heights to 'fill whatever space I have'
  p2$widths <- grid::unit(rep(1, ncol(p2)), "null")
  p2$heights <- grid::unit(rep(1, nrow(p2)), "null")

  p3 <- ggplot2::ggplot() +
    ggplot2::annotation_custom(p2)

  if(!is.null(section_heights)) {
    CombinedPlot <- egg::ggarrange(top_counts, p3, bottom_stacked, heights = section_heights, draw = FALSE)
  } else {
    CombinedPlot <- egg::ggarrange(top_counts, p3, bottom_stacked, heights = c(2, 1.25, 2), draw = FALSE)
  }

  print(CombinedPlot)

  if (save_table) {
    dev.copy(png, res = figure_res, width = figure_width, height = figure_height, paste0(save_folder, "primer_variant_vis.png"))
    dev.off()
  }
}
