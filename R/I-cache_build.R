# strings list
s <- new.env()
# string list for pull_subset_DB
subsetStrings <- new.env()
# Program cache
cache <- new.env()

clearCache <- function() {
  rm(list = ls(envir = cache), envir = cache)
}


clearStrings <- function() {
  rm(list = ls(envir = s), envir = s)
}

setCacheConsts <- function() {
  cache$alignmentBreaktime <- .5
}

setStrings2AssayDB <- function() {
  # TABLE NAMES--------------------------------------
  s$AMBIGUOUS_SEQS_TABLE <- "Ambiguous_Seqs"
  s$PRIMER_TABLE <- "Primer_Info"
  s$ASSAY_TABLE <- "Assay_Info"
  s$COSPOT_KEY <- "Cospot_Key"
  s$HIGH_N_TABLE <- "High_N_Seqs"
  s$PRIMER_LOCS_TABLE <- "Primer_Locations"
  s$SEQ_INFO_TABLE <- "Seq_Info"
  s$SEQUENCES_TABLE <- "Sequences"
  s$VAR_INFO_TABLE <- "Variant_Info"
  s$ID_MAP_TABLE <- "ID_Map"
  s$VAR_TYPE_TABLE <- "Variant_Type"
  s$AMPLICON_TABLE <- "Var_Amplicon_Info"
  s$SEQ_META_TABLE <- "Seq_Meta"
  s$SIG_VAR_TABLE <- "Sig_Var_Ids"

  # Column Names-------------------------------------
  s$ACC_ID <- "Accession_ID"
  s$PRIMER_ID <- "Primer_ID"
  s$ASSAY_ID <- "Assay_ID"
  s$VAR_ID <- "Var_ID"

  s$VAR_SEQ_COL <- "Var_Seq"
  s$INDELS_PRESENT_VAR <- "Has_Indel"

  # Cospot Map
  s$SPOT1 <- "Spot"
  s$SPOT2 <- "Cospot"

  # Primers
  s$SEARCHABLE_SEQ <- "Searchable_Seq"
  s$REACTION_COL <- "Reaction"
  s$IS_COSPOT_STRING <- "is_Cospot"
  s$DIRECTION_COL <- "Direction"
  s$P_SEQ_COL <- "Primer_Sequence"
  s$S_SEQ_COL <- "Searchable_Seq"

  # Primers Ref only
  s$PRIMER_NAME <- "Primer_Name"
  s$ASSAY_DEV_NAME <- "Development_Name"
  s$ASSAY_NAME <- "Assay_Name"

  # Location
  s$LOC_START <- "Start"
  s$LOC_END <- "End"

  # Var Info
  s$MISMATCH_NUM <- "Num_Mismatches"
  s$VAR_COUNT <- "Var_Count"

  # Variant Type
  s$VAR_TYPE <- "Var_Type"
  s$SEQ_CHANGE <- "Seq_Change"
  s$POS_FROM_SIDE <- "Pos_from_3P"
  s$MISMATCH_TYPE <- "Mismatch_Type"
  s$INDEL_LEN_COL <- "Indel_Length"

  # SEQ INFO
  s$LOC_ACC <- "Location"
  s$FASTA_HEADER <- "FASTA_Header"
  s$PER_N <- "Percent_N"
  s$HAS_AMB_COL <- "Has_Ambiguities"
  s$SEQ_VAR_COUNT <- "Var_Count"
  s$ASSAYS_AFFECTED <- "Assays_Affected"
  s$BP10_AA <- "HR_Assays_Affected"
  s$WITH_AMBS_AA <- "WC_Assays_Affected"
  s$HIGH_RISK_POS <- 10
  s$COLLECTION_DATE <- "Collection_Date"
  s$COLLECTION_DAY <- "Collection_Day"
  s$PULL_DATE <- "Pull_Date"
  s$LINEAGE <- "Lineage"
  s$CLADE <- "Clade"
  s$GROUP <- "WHO_Designation"
  s$SUBGROUP <- "Subgroup"
  s$DISPLAY_LIN <- "Displayed_Lineage"
  s$DISPLAY_GROUP <- "Displayed_Group"

  # ASSAY INFO AND BAD AMPLICONS
  s$INNER_AMP_LEN <- "Inner_Amplicon_Length"
  s$OUTER_AMP_LEN <- "Outer_Amplicon_Length"
  s$INNER_GC <- "GC_Content"

  # BAD AMPLICONS
  s$OUT_OF_ORDER <- "Out_of_Order"
  s$DUPLICATED <- "Duplicated"

  # AMBIGOUS SEQS
  s$AMB_CHARS_COL <- "Ambiguity_Chars"
  s$SEQ_PER_N <- "Seq_Percent_N"

  # METADATA
  s$CONTINENT <- "Continent"
  s$COUNTRY <- "Country"
  s$LOC_DETAILS <- "Location_Details"

  # SIG VAR IDs (subset DB)
  s$SIG_VAR_CUTOFF <- 0.001

  # Primer index df in cache
  s$ACC_PRIMER_S <- "Acc_Primer_Start"
  s$ACC_PRIMER_E <- "Acc_Primer_End"

  # Adaptions
  s$P_VAR_ID <- s$VAR_ID
  s$P_VAR_INFO_TABLE <- s$VAR_INFO_TABLE
  s$LOC_ID <- s$ASSAY_ID

  s$NEW_COL_NAME <- paste0("\U2265", "6")
  s$NEW_PERCENT_COL_NAME <- paste0({{ s$NEW_COL_NAME }}, "_", "%")
}

setGISAIDstrings <- function() {
  s$METADATA_ACC_ID <- "Accession ID"
  s$GISAID_CLADE <- "Clade"
  s$PANGO_LINEAGE <- "Lineage"

  s$NON_HUMAN_CHECK <- c(
    "bat", "pangolin", "env", "lion", "dog", "monkey", "cat",
    "gorilla", "hamster", "mouse", "mink", "tiger", "leopard", "canine", "deer"
  )

  # GISAID Index Parsing Constants
  s$LOCATION_INDEX <- 2
  s$ACCESSION_INDEX <- 5 # Not currently in use. Regular expression below is used
  s$SUBMISSION_DATE_INDEX <- 6 # Not currently in use. Regular expression below is used

  # REGULAR EXPRESSIONS
  s$ACCESSION_REG_EXP <- "[|](EPI[[^|]]+)"
  s$COL_DATE_REG_EXP <- "[[^|]]+$"
  s$STRING_SPLIT_EXPRESSION <- "[/|]+"
}

#' Set strings in "subset" environment for pull_subset_DB function
#'
setSubsetStrings <- function() {
  # Add Assay and Primers
  subsetStrings$build_primers <- glue::glue("INSERT INTO {s$ASSAY_TABLE}
SELECT *
  FROM archive.{s$ASSAY_TABLE};

INSERT INTO {s$PRIMER_TABLE}
SELECT {s$PRIMER_ID}, {s$ASSAY_NAME}, {s$ASSAY_DEV_NAME}, {s$ASSAY_ID},
{s$PRIMER_NAME}, {s$DIRECTION_COL}, {s$REACTION_COL}, {s$P_SEQ_COL}, {s$IS_COSPOT_STRING},
{s$S_SEQ_COL}, is_Masked
FROM archive.{s$PRIMER_TABLE};")

  # Seq_Info Table Insert
  subsetStrings$build_seq_info <- glue::glue("INSERT INTO {s$SEQ_INFO_TABLE}
SELECT {s$ACC_ID}, {s$LOC_ACC}, {s$COLLECTION_DATE}, {s$COLLECTION_DAY},
{s$LINEAGE}, {s$CLADE}, {s$FASTA_HEADER}, {s$PULL_DATE},
{s$PER_N}, {s$HAS_AMB_COL}, {s$INDELS_PRESENT_VAR},
{s$SEQ_VAR_COUNT}, {s$ASSAYS_AFFECTED}, {s$BP10_AA}, {s$WITH_AMBS_AA}
  FROM archive.{s$SEQ_INFO_TABLE}")

  subsetStrings$addCol <- glue::glue("ALTER TABLE {s$SEQ_INFO_TABLE} ADD COLUMN {s$GROUP};")

  subsetStrings$addCol2 <- glue::glue("ALTER TABLE {s$SEQ_INFO_TABLE} ADD COLUMN {s$SUBGROUP};")

  subsetStrings$addCol3 <- glue::glue("ALTER TABLE {s$SEQ_INFO_TABLE} ADD COLUMN {s$DISPLAY_LIN};")

  subsetStrings$addCol4 <- glue::glue("ALTER TABLE {s$ID_MAP_TABLE} ADD COLUMN {s$DISPLAY_LIN};")

  subsetStrings$addCol5 <- glue::glue("ALTER TABLE {s$SEQ_INFO_TABLE} ADD COLUMN {s$DISPLAY_GROUP};")

  subsetStrings$addCol6 <- glue::glue("ALTER TABLE {s$ID_MAP_TABLE} ADD COLUMN {s$DISPLAY_GROUP};")

  subsetStrings$addIDmap <- glue::glue("PRAGMA foreign_keys=off;

INSERT INTO {s$ID_MAP_TABLE} ({s$ACC_ID}, {s$GROUP}, {s$SUBGROUP}, {s$LINEAGE}, {s$VAR_ID}, {s$ASSAY_ID}, {s$PRIMER_ID})
SELECT SL.{s$ACC_ID}, {s$GROUP}, {s$SUBGROUP}, {s$LINEAGE}, {s$VAR_ID}, {s$ASSAY_ID}, {s$PRIMER_ID}
FROM (SELECT {s$ACC_ID}, {s$COLLECTION_DATE}, {s$GROUP}, {s$SUBGROUP}, {s$LINEAGE} FROM {s$SEQ_INFO_TABLE}) AS SL
	INNER JOIN archive.{s$ID_MAP_TABLE}  AS VL ON SL.{s$ACC_ID} = VL.{s$ACC_ID};

PRAGMA foreign_keys=on;")


  subsetStrings$addDisplayIDmap <- glue::glue(
    "UPDATE {s$ID_MAP_TABLE}
  SET {s$DISPLAY_LIN}=(SELECT {s$DISPLAY_LIN} FROM {s$SEQ_INFO_TABLE} WHERE {s$ID_MAP_TABLE}.{s$ACC_ID}={s$SEQ_INFO_TABLE}.{s$ACC_ID}),
  {s$DISPLAY_GROUP}=(SELECT {s$DISPLAY_GROUP} FROM {s$SEQ_INFO_TABLE} WHERE {s$ID_MAP_TABLE}.{s$ACC_ID}={s$SEQ_INFO_TABLE}.{s$ACC_ID});"
  )

  # The new Sig_Var_Ids table contains 1 column with all the Var_ID's that occur at a user-defined frequency.
  # It has Id's with var counts (from var info table) that are >= total number of accessions * user-defined sig var cutoff
  subsetStrings$insertOtherVarInfo <- glue::glue("INSERT INTO {s$VAR_INFO_TABLE}
SELECT VI.{s$VAR_ID}, {s$VAR_SEQ_COL}, {s$VAR_COUNT}, {s$MISMATCH_NUM},
{s$INDELS_PRESENT_VAR}, {s$PRIMER_ID}, {s$ASSAY_ID}
FROM archive.{s$VAR_INFO_TABLE} AS VI JOIN  (
	SELECT {s$VAR_ID}, COUNT(*) AS {s$VAR_COUNT}
	FROM {s$ID_MAP_TABLE}
	GROUP BY {s$VAR_ID}) AS FV ON VI.{s$VAR_ID} = FV.{s$VAR_ID};

INSERT INTO {s$VAR_TYPE_TABLE}
SELECT DISTINCT *
FROM archive.{s$VAR_TYPE_TABLE}
WHERE {s$VAR_ID} IN (
	SELECT {s$VAR_ID} FROM {s$VAR_INFO_TABLE});


INSERT INTO {s$SIG_VAR_TABLE} ({s$VAR_ID})
SELECT {s$VAR_ID}
FROM {s$VAR_INFO_TABLE}
WHERE {s$VAR_COUNT} >= (
	SELECT COUNT({s$ACC_ID}) * {s$SIG_VAR_CUTOFF} AS Total_Seqs
	FROM {s$SEQ_INFO_TABLE});

INSERT INTO {s$AMBIGUOUS_SEQS_TABLE}
SELECT DISTINCT *
FROM archive.{s$AMBIGUOUS_SEQS_TABLE}
WHERE {s$ACC_ID} IN (SELECT {s$ACC_ID} FROM {s$SEQ_INFO_TABLE});

INSERT INTO {s$SEQ_META_TABLE}
SELECT DISTINCT *
FROM archive.{s$SEQ_META_TABLE}
WHERE {s$ACC_ID} IN (SELECT {s$ACC_ID} FROM {s$SEQ_INFO_TABLE});")
}
