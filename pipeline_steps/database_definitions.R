studies_db <- list(
  study_sources = list(
    name = "study_sources",
    query = "CREATE TABLE study_sources (
      id INTEGER PRIMARY KEY,
      source TEXT NOT NULL,
      name TEXT NOT NULL,
      url TEXT NOT NULL,
      doi TEXT NOT NULL
    )"
  ),
  ld_blocks = list(
    name = "ld_blocks",
    query = "CREATE TABLE ld_blocks (
      id INTEGER PRIMARY KEY,
      chr INTEGER NOT NULL,
      start INTEGER NOT NULL,
      stop INTEGER NOT NULL,
      ancestry TEXT NOT NULL,
      ld_block TEXT NOT NULL
    )"
  ),
  studies = list(
    name = "studies",
    query = "CREATE TABLE studies (
      id INTEGER PRIMARY KEY,
      data_type TEXT NOT NULL,
      data_format TEXT NOT NULL,
      study_name TEXT NOT NULL,
      trait TEXT NOT NULL,
      ancestry TEXT NOT NULL,
      sample_size INTEGER NOT NULL,
      category TEXT,
      study_location TEXT,
      extracted_location TEXT,
      probe TEXT,
      tissue TEXT,
      source_id INTEGER ,
      variant_type TEXT,
      p_value_threshold DOUBLE,
      gene TEXT,
      FOREIGN KEY (source_id) REFERENCES study_sources(id)
    )"
  ),
  snp_annotations = list(
    name = "snp_annotations",
    query = "CREATE TABLE snp_annotations (
      id INTEGER PRIMARY KEY,
      snp TEXT,
      chr INTEGER,
      bp INTEGER,
      ea TEXT,
      oa TEXT,
      gene TEXT,
      feature_type TEXT,
      consequence TEXT,
      cdna_position TEXT,
      cds_position TEXT,
      protein_position TEXT,
      amino_acids TEXT,
      codons TEXT,
      rsid TEXT,
      impact TEXT,
      symbol TEXT,
      biotype TEXT,
      strand TEXT,
      canonical TEXT,
      all_af REAL CHECK (all_af BETWEEN 0 AND 1),
      eur_af REAL CHECK (eur_af BETWEEN 0 AND 1),
      amr_af REAL CHECK (amr_af BETWEEN 0 AND 1),
      eas_af REAL CHECK (eas_af BETWEEN 0 AND 1),
      sas_af REAL CHECK (sas_af BETWEEN 0 AND 1),
      afr_af REAL CHECK (afr_af BETWEEN 0 AND 1)
    )"
  ),
  study_extractions = list(
    name = "study_extractions",
    query = "CREATE TABLE study_extractions (
      id INTEGER PRIMARY KEY,
      study_id INTEGER,
      snp_id INTEGER,
      snp TEXT NOT NULL,
      ld_block_id INTEGER,
      unique_study_id TEXT NOT NULL,
      study TEXT NOT NULL,
      file TEXT NOT NULL,
      chr INTEGER NOT NULL,
      bp INTEGER NOT NULL,
      min_p DOUBLE CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block TEXT,
      known_gene TEXT,
      FOREIGN KEY (study_id) REFERENCES studies(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  ),
  colocalisations = list(
    name = "colocalisations",
    query = "CREATE TABLE colocalisations (
      study_extraction_id INTEGER,
      snp_id INTEGER,
      ld_block_id INTEGER,
      coloc_group_id INTEGER,
      iteration INTEGER,
      unique_study_id TEXT,
      posterior_prob REAL CHECK (posterior_prob BETWEEN 0 AND 1),
      regional_prob REAL CHECK (regional_prob BETWEEN 0 AND 1),
      posterior_explained_by_snp REAL CHECK (posterior_explained_by_snp BETWEEN 0 AND 1),
      candidate_snp TEXT,
      study_id INTEGER,
      chr INTEGER NOT NULL,
      bp INTEGER NOT NULL,
      min_p DOUBLE CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block TEXT,
      known_gene TEXT,
      PRIMARY KEY (study_extraction_id, snp_id),
      FOREIGN KEY (study_extraction_id) REFERENCES study_extractions(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  ),
  rare_results = list(
    name = "rare_results",
    query = "CREATE TABLE rare_results (
      rare_result_group_id INTEGER,
      study_extraction_id INTEGER,
      snp_id INTEGER,
      ld_block_id INTEGER,
      unique_study_id TEXT,
      candidate_snp TEXT,
      study_id INTEGER,
      chr INTEGER,
      bp INTEGER,
      min_p DOUBLE CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      known_gene TEXT,
      FOREIGN KEY (study_extraction_id) REFERENCES study_extractions(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  ),
  results_metadata = list(
    name = "results_metadata",
    query = "CREATE TABLE results_metadata (
      ld_block_id INTEGER,
      ld_block TEXT,
      number_extracted INTEGER,
      number_standardised INTEGER,
      mean_snps_removed_by_reference_panel REAL,
      number_imputed INTEGER,
      significant_snps_imputed INTEGER,
      significant_imputed_snps_filtered INTEGER,
      number_finemapped INTEGER,
      finemapped_per_imputed REAL,
      num_finemap_failed INTEGER,
      standardised_time_taken REAL,
      imputed_time_taken REAL,
      finemapped_time_taken REAL,
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )"
  )
)

ld_table <- list(
  name = "ld",
  query = "CREATE TABLE ld (
    lead_snp_id INTEGER,
    variant_snp_id INTEGER,
    ld_block_id INTEGER,
    r REAL CHECK (r BETWEEN -1 AND 1)
  )"
)

#TODO: add in se > 0 check later
associations_table <- list(
  name = "associations",
  query = "CREATE TABLE associations (
    snp_id INTEGER,
    study_id INTEGER,
    beta REAL,
    se REAL,
    p DOUBLE CHECK (p BETWEEN 0 AND 1),
    eaf REAL CHECK (eaf BETWEEN 0 AND 1),
    imputed BOOLEAN
  )"
)

gwas_upload_db <- list(
  gwas_upload = list(
    name = "gwas_upload",
    query = "CREATE SEQUENCE id_sequence START 1; CREATE TABLE gwas_upload (
      id INTEGER PRIMARY KEY DEFAULT nextval('id_sequence'),
      guid TEXT UNIQUE NOT NULL,
      email TEXT NOT NULL,
      name TEXT NOT NULL,
      sample_size INTEGER NOT NULL,
      ancestry TEXT NOT NULL,
      category TEXT NOT NULL,
      is_published BOOLEAN NOT NULL,
      doi TEXT,
      should_be_added BOOLEAN,
      status TEXT NOT NULL,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )"
  ),
  study_extractions = list(
    name = "study_extractions",
    query = "CREATE TABLE study_extractions (
      id INTEGER PRIMARY KEY,
      study_id INTEGER,
      snp_id INTEGER,
      snp TEXT NOT NULL,
      ld_block_id INTEGER,
      unique_study_id TEXT NOT NULL,
      study TEXT NOT NULL,
      file TEXT NOT NULL,
      chr INTEGER NOT NULL,
      bp INTEGER NOT NULL,
      min_p DOUBLE CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block TEXT,
      known_gene TEXT)"
  ),
  colocalisations = list(
    name = "colocalisations",
    query = "CREATE TABLE colocalisations (
      study_extraction_id INTEGER,
      snp_id INTEGER,
      ld_block_id INTEGER,
      coloc_group_id INTEGER,
      iteration INTEGER,
      unique_study_id TEXT,
      posterior_prob REAL CHECK (posterior_prob BETWEEN 0 AND 1),
      regional_prob REAL CHECK (regional_prob BETWEEN 0 AND 1),
      posterior_explained_by_snp REAL CHECK (posterior_explained_by_snp BETWEEN 0 AND 1),
      candidate_snp TEXT,
      study_id INTEGER,
      chr INTEGER NOT NULL,
      bp INTEGER NOT NULL,
      min_p DOUBLE CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block TEXT,
      known_gene TEXT)"
  )
)