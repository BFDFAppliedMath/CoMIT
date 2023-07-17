-- !preview conn=DBI::dbConnect(RSQLite::SQLite())

PRAGMA foreign_keys=off;

DROP TABLE IF EXISTS Seq_Info;

DROP TABLE IF EXISTS Assay_Info;

DROP TABLE IF EXISTS Primer_Info;

DROP TABLE IF EXISTS Variant_Info;

DROP TABLE IF EXISTS ID_Map;

DROP TABLE IF EXISTS Ambiguous_Seqs;

DROP TABLE IF EXISTS Variant_Type;

DROP TABLE IF EXISTS Seq_Meta;

DROP TABLE IF EXISTS Sig_Var_Ids;

CREATE TABLE Seq_Info (
  Accession_ID varchar,
  Location varchar,
  Collection_Date varchar,
  Collection_Day varchar,
  Lineage varchar,
  Clade varchar,
  FASTA_Header varchar,
  Pull_Date varchar,
  Percent_N numeric DEFAULT 0,
  Has_Ambiguities integer DEFAULT 0,
  Has_Indel integer DEFAULT 0,
  Var_Count integer DEFAULT 0,
  Assays_Affected integer DEFAULT 0,
  HR_Assays_Affected integer DEFAULT 0,
  WC_Assays_Affected integer DEFAULT 0,
  PRIMARY KEY(Accession_ID)
);


CREATE TABLE Assay_Info (
  Assay_ID integer PRIMARY KEY AUTOINCREMENT,
  Assay_Name varchar,
  Development_Name varchar,
  is_Masked integer NOT NULL DEFAULT 0,
  Inner_Amplicon_Length integer,
  Outer_Amplicon_Length integer,
  GC_Content real
);

CREATE TABLE Primer_Info (
  Primer_ID integer PRIMARY KEY AUTOINCREMENT,
  Assay_Name varchar,
  Development_Name varchar,
  Assay_ID integer,
  Primer_Name varchar,
  Direction varchar,
  Reaction varchar,
  Primer_Sequence varchar,
  is_Cospot integer NOT NULL DEFAULT 0,
  Searchable_Seq varchar,
  is_Masked integer NOT NULL DEFAULT 0,
  FOREIGN KEY(Assay_ID) REFERENCES Assay_Info(Assay_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE
);

CREATE TABLE Variant_Info (
  Var_ID integer PRIMARY KEY AUTOINCREMENT,
  Var_Seq text,
  Var_Count integer NOT NULL DEFAULT 0,
  Num_Mismatches integer NOT NULL DEFAULT 0,
  Has_Indel integer NOT NULL DEFAULT 0,
  Primer_ID integer,
  Assay_ID integer,
  FOREIGN KEY(Primer_ID) REFERENCES Primer_Info(Primer_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(Assay_ID) REFERENCES Assay_Info(Assay_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE
);

CREATE TABLE ID_Map (
  Accession_ID varchar,
  WHO_Designation varchar,
  Subgroup varchar,
  Lineage varchar,
  Var_ID integer,
  Assay_ID integer,
  Primer_ID integer,
  FOREIGN KEY(Accession_ID) REFERENCES Seq_Info(Accession_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(Primer_ID) REFERENCES Primer_Info(Primer_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(Var_ID) REFERENCES Variant_Info(Var_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(Assay_ID) REFERENCES Assay_Info(Assay_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE
);

CREATE TABLE Ambiguous_Seqs (
  Accession_ID varchar,
  Var_Seq text,
  Ambiguity_Chars varchar,
  Seq_Percent_N numeric NOT NULL DEFAULT 0,
  Primer_ID integer,
  FOREIGN KEY(Primer_ID) REFERENCES Primer_Info(Primer_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(Accession_ID) REFERENCES Seq_Info(Accession_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE
);

CREATE TABLE Variant_Type (
  Var_ID integer,
  Var_Type varchar,
  Seq_Change varchar,
  Mismatch_Type varchar,
  Pos_from_3P integer,
  Primer_ID integer,
  Indel_Length integer,
  FOREIGN KEY(Var_ID) REFERENCES Variant_Info(Var_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(Primer_ID) REFERENCES Primer_Info(Primer_ID)
  ON DELETE CASCADE
  ON UPDATE CASCADE
);

CREATE TABLE Seq_Meta (
  Accession_ID varchar,
  Continent text,
  Country text,
  Location_Details text
);

CREATE TABLE Sig_Var_Ids (
	Var_ID integer
	,FOREIGN KEY(Var_ID) REFERENCES Variant_Info(Var_ID)
		ON DELETE CASCADE
		ON UPDATE CASCADE
);

PRAGMA foreign_keys=on;
