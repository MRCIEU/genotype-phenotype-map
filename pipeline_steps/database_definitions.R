studies_db <- list(
  study_sources = list(
    name = "study_sources",
    persist_id_from = "source",
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
    persist_id_from = "ld_block",
    query = "CREATE TABLE ld_blocks (
      id INTEGER PRIMARY KEY,
      chr INTEGER NOT NULL,
      start INTEGER NOT NULL,
      stop INTEGER NOT NULL,
      ancestry TEXT NOT NULL,
      ld_block TEXT NOT NULL
    )"
  ),
  gene_annotations = list(
    name = "gene_annotations",
    persist_id_from = "gene",
    query = "CREATE TABLE gene_annotations (
      id INTEGER PRIMARY KEY,
      ensembl_id TEXT NOT NULL,
      gene TEXT NOT NULL,
      description TEXT,
      gene_biotype TEXT,
      chr INTEGER NOT NULL,
      start INTEGER NOT NULL,
      stop INTEGER NOT NULL,
      strand INTEGER,
      source TEXT
    )"
  ),
  traits = list(
    name = "traits",
    persist_id_from = "trait",
    query = glue::glue("CREATE TABLE traits (
      id INTEGER PRIMARY KEY,
      data_type TEXT NOT NULL,
      trait TEXT NOT NULL,
      trait_name TEXT NOT NULL,
      trait_category TEXT
    )")
  ),
  studies = list(
    name = "studies",
    persist_id_from = "study_name",
    query = glue::glue("CREATE TABLE studies (
      id INTEGER PRIMARY KEY,
      data_type TEXT NOT NULL,
      study_name TEXT NOT NULL,
      trait_id INTEGER NOT NULL,
      ancestry TEXT NOT NULL,
      sample_size INTEGER NOT NULL,
      category TEXT NOT NULL,
      probe TEXT,
      tissue TEXT,
      cell_type TEXT,
      source_id INTEGER,
      variant_type TEXT,
      p_value_threshold DOUBLE,
      gene_id INTEGER,
      gene TEXT,
      ensg TEXT,
      heritability REAL,
      heritability_se REAL,
      FOREIGN KEY (trait_id) REFERENCES traits(id),
      FOREIGN KEY (source_id) REFERENCES study_sources(id),
      FOREIGN KEY (gene_id) REFERENCES gene_annotations(id)
    )")
  ),
  snp_annotations = list(
    name = "snp_annotations",
    persist_id_from = "snp",
    query = "CREATE TABLE snp_annotations (
      id INTEGER PRIMARY KEY,
      snp TEXT NOT NULL,
      display_snp TEXT,
      chr INTEGER NOT NULL,
      bp INTEGER NOT NULL,
      ea TEXT NOT NULL,
      oa TEXT NOT NULL,
      ref_allele TEXT,
      flipped BOOLEAN,
      gene_id INTEGER,
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
      afr_af REAL CHECK (afr_af BETWEEN 0 AND 1),
      FOREIGN KEY (gene_id) REFERENCES gene_annotations(id)
    )"
  ),
  study_extractions = list(
    name = "study_extractions",
    persist_id_from = "unique_study_id",
    query = "CREATE TABLE study_extractions (
      id INTEGER PRIMARY KEY,
      study_id INTEGER NOT NULL,
      snp_id INTEGER NOT NULL,
      display_snp TEXT NOT NULL NOT NULL,
      rsid TEXT NOT NULL NOT NULL,
      ld_block_id INTEGER NOT NULL,
      unique_study_id TEXT NOT NULL,
      study TEXT NOT NULL,
      file TEXT NOT NULL,
      svg_file TEXT,
      file_with_lbfs TEXT,
      chr INTEGER NOT NULL,
      bp INTEGER NOT NULL,
      min_p DOUBLE CHECK (min_p BETWEEN 0 AND 1),
      cis_trans TEXT,
      ld_block TEXT,
      gene_id INTEGER,
      situated_gene_id INTEGER,
      FOREIGN KEY (study_id) REFERENCES studies(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id),
      FOREIGN KEY (gene_id) REFERENCES gene_annotations(id),
      FOREIGN KEY (situated_gene_id) REFERENCES gene_annotations(id)
    )",
    indexes = "CREATE INDEX idx_study_extractions_study_id ON study_extractions(study_id);"
  ),
  coloc_groups = list(
    name = "coloc_groups",
    persist_id_from = NA,
    query = glue::glue("CREATE TABLE coloc_groups (
      coloc_group_id INTEGER NOT NULL,
      study_id INTEGER NOT NULL,
      study_extraction_id INTEGER NOT NULL,
      snp_id INTEGER NOT NULL,
      ld_block_id INTEGER NOT NULL,
      h4_connectedness REAL,
      h3_connectedness REAL,
      PRIMARY KEY (coloc_group_id, study_extraction_id),
      FOREIGN KEY (study_id) REFERENCES studies(id),
      FOREIGN KEY (study_extraction_id) REFERENCES study_extractions(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id)
    )")
  ),
  rare_results = list(
    name = "rare_results",
    persist_id_from = NA,
    query = "CREATE TABLE rare_results (
      rare_result_group_id INTEGER,
      study_id INTEGER NOT NULL,
      study_extraction_id INTEGER NOT NULL,
      snp_id INTEGER NOT NULL,
      gene_id INTEGER,
      situated_gene_id INTEGER,
      ld_block_id INTEGER NOT NULL,
      FOREIGN KEY (study_id) REFERENCES studies(id),
      FOREIGN KEY (study_extraction_id) REFERENCES study_extractions(id),
      FOREIGN KEY (snp_id) REFERENCES snp_annotations(id),
      FOREIGN KEY (ld_block_id) REFERENCES ld_blocks(id),
      FOREIGN KEY (gene_id) REFERENCES gene_annotations(id),
      FOREIGN KEY (situated_gene_id) REFERENCES gene_annotations(id)
    )"
  ),
  snp_pleiotropy = list(
    name = "snp_pleiotropy",
    persist_id_from = NA,
    query = "CREATE TABLE snp_pleiotropy (
      snp_id INTEGER PRIMARY KEY NOT NULL,
      distinct_trait_categories INTEGER NOT NULL,
      distinct_protein_coding_genes INTEGER NOT NULL
    )"
  ),
  gene_pleiotropy = list(
    name = "gene_pleiotropy",
    persist_id_from = NA,
    query = "CREATE TABLE gene_pleiotropy (
      gene_id INTEGER PRIMARY KEY NOT NULL,
      distinct_trait_categories INTEGER NOT NULL,
      distinct_protein_coding_genes INTEGER NOT NULL
    )"
  )
)

additional_studies_tables <- list(
  coloc_groups_wide = list(
    name = "coloc_groups_wide",
    query = "CREATE TABLE coloc_groups_wide AS 
      SELECT coloc_groups.*, 
        snp_annotations.chr, snp_annotations.bp, study_extractions.min_p, study_extractions.cis_trans,
        study_extractions.ld_block, snp_annotations.display_snp, snp_annotations.rsid, gene_annotations.gene, gene_annotations.id as gene_id,
        traits.id as trait_id, traits.trait_name, traits.trait_category, studies.data_type, studies.tissue, studies.cell_type,
        study_sources.id as source_id, study_sources.name as source_name, study_sources.url as source_url
      FROM coloc_groups 
      JOIN studies ON coloc_groups.study_id = studies.id
      JOIN snp_annotations on coloc_groups.snp_id = snp_annotations.id 
      JOIN study_extractions ON coloc_groups.study_extraction_id = study_extractions.id
      LEFT JOIN gene_annotations on studies.gene_id = gene_annotations.id
      JOIN traits ON studies.trait_id = traits.id
      JOIN study_sources ON studies.source_id = study_sources.id",
    indexes = "CREATE INDEX idx_coloc_groups_wide_study_id ON coloc_groups_wide(study_id);
      CREATE INDEX idx_coloc_groups_wide_study_extraction_id ON coloc_groups_wide(study_extraction_id);"
  ),
  rare_results_wide = list(
    name = "rare_results_wide",
    query = "CREATE TABLE rare_results_wide AS 
      SELECT rare_results.*,
        study_extractions.chr, study_extractions.bp, study_extractions.min_p, study_extractions.cis_trans, snp_annotations.display_snp, snp_annotations.rsid,
        gene_annotation.gene AS gene, situated_gene_annotation.gene AS situated_gene, traits.id as trait_id, traits.trait_name, traits.trait_category, studies.data_type, studies.tissue, studies.cell_type,
        ld_blocks.ld_block, study_sources.id as source_id, study_sources.name as source_name, study_sources.url as source_url
      FROM rare_results
      JOIN studies ON rare_results.study_id = studies.id
      JOIN snp_annotations ON rare_results.snp_id = snp_annotations.id
      JOIN study_extractions ON rare_results.study_extraction_id = study_extractions.id
      LEFT JOIN gene_annotations AS gene_annotation ON rare_results.gene_id = gene_annotation.id
      LEFT JOIN gene_annotations AS situated_gene_annotation ON rare_results.situated_gene_id = situated_gene_annotation.id
      JOIN traits ON studies.trait_id = traits.id
      JOIN study_sources ON studies.source_id = study_sources.id
      JOIN ld_blocks ON rare_results.ld_block_id = ld_blocks.id",
    indexes = "CREATE INDEX idx_rare_results_wide_study_id ON rare_results_wide(study_id);
      CREATE INDEX idx_rare_results_wide_study_extraction_id ON rare_results_wide(study_extraction_id);"
  )
)

ld_table <- list(
  name = "ld",
  query = "CREATE TABLE ld (
    lead_snp_id INTEGER,
    variant_snp_id INTEGER,
    ld_block_id INTEGER,
    r REAL CHECK (r BETWEEN -1 AND 1),
    PRIMARY KEY (lead_snp_id, variant_snp_id)
  )"
)

coloc_pairs_full_table <- list(
  name = "coloc_pairs",
  query = "CREATE TABLE coloc_pairs (
    study_extraction_a_id INTEGER NOT NULL,
    study_extraction_b_id INTEGER NOT NULL,
    ld_block_id INTEGER NOT NULL,
    h0 REAL CHECK (h0 BETWEEN 0 AND 1) NOT NULL,
    h1 REAL CHECK (h1 BETWEEN 0 AND 1) NOT NULL,
    h2 REAL CHECK (h2 BETWEEN 0 AND 1) NOT NULL,
    h3 REAL CHECK (h3 BETWEEN 0 AND 1) NOT NULL,
    h4 REAL CHECK (h4 BETWEEN 0 AND 1) NOT NULL,
    false_positive BOOLEAN NOT NULL,
    false_negative BOOLEAN NOT NULL,
    PRIMARY KEY (study_extraction_a_id, study_extraction_b_id)
  )"
)

coloc_pairs_significant_table <- list(
  name = "coloc_pairs",
  query = "CREATE TABLE coloc_pairs (
    snp_id INTEGER NOT NULL,
    study_extraction_a_id INTEGER NOT NULL,
    study_extraction_b_id INTEGER NOT NULL,
    ld_block_id INTEGER NOT NULL,
    h3 REAL CHECK (h3 BETWEEN 0 AND 1) NOT NULL,
    h4 REAL CHECK (h4 BETWEEN 0 AND 1) NOT NULL,
    false_positive BOOLEAN NOT NULL,
    false_negative BOOLEAN NOT NULL,
    PRIMARY KEY (study_extraction_a_id, study_extraction_b_id)
  )",
  indexes = "CREATE INDEX idx_coloc_pairs_study_extraction_a_id ON coloc_pairs (study_extraction_a_id);
    CREATE INDEX idx_coloc_pairs_study_extraction_b_id ON coloc_pairs (study_extraction_b_id);
    CREATE INDEX idx_coloc_pairs_ld_block_id ON coloc_pairs (ld_block_id);
    CREATE INDEX idx_coloc_pairs_snp_id ON coloc_pairs (snp_id);"
)

coloc_pairs_significant_db <- list(
  coloc_pairs_metadata = list(
    name = "coloc_pairs_metadata",
    query = "CREATE TABLE coloc_pairs_metadata (
      start_snp_id INTEGER,
      stop_snp_id INTEGER,
      coloc_pairs_table_name TEXT NOT NULL
    )"
  ),
  coloc_pairs = list(
  name = "coloc_pairs",
  query = "CREATE TABLE table_name (
    snp_id INTEGER NOT NULL,
    study_extraction_a_id INTEGER NOT NULL,
    study_extraction_b_id INTEGER NOT NULL,
    ld_block_id INTEGER NOT NULL,
    h3 REAL CHECK (h3 BETWEEN 0 AND 1) NOT NULL,
    h4 REAL CHECK (h4 BETWEEN 0 AND 1) NOT NULL,
    false_positive BOOLEAN NOT NULL,
    false_negative BOOLEAN NOT NULL,
    PRIMARY KEY (study_extraction_a_id, study_extraction_b_id)
  )",
  indexes = "CREATE INDEX idx_table_name_study_extraction_a_id ON table_name (study_extraction_a_id);
    CREATE INDEX idx_table_name_study_extraction_b_id ON table_name (study_extraction_b_id);
    CREATE INDEX idx_table_name_ld_block_id ON table_name (ld_block_id);
    CREATE INDEX idx_table_name_snp_id ON table_name (snp_id);"
  )
)

associations_db <- list(
  associations_metadata = list(
    name = "associations_metadata",
    query = "CREATE TABLE associations_metadata (
      start_snp_id INTEGER,
      stop_snp_id INTEGER,
      associations_table_name TEXT NOT NULL
    )"
  ),
  associations = list(
    name = "associations",
    query = "CREATE TABLE table_name (
      snp_id INTEGER,
      study_id INTEGER,
      beta REAL NOT NULL,
      se REAL NOT NULL CHECK (se > 0),
      p DOUBLE CHECK (p BETWEEN 0 AND 1) NOT NULL,
      eaf REAL CHECK (eaf BETWEEN 0 AND 1) NOT NULL,
      imputed BOOLEAN NOT NULL,
      PRIMARY KEY (snp_id, study_id)
    )"
  )
)

gwas_upload_db <- list(
  gwas_upload = list(
    name = "gwas_upload",
    query = "CREATE SEQUENCE gwas_upload_id_sequence START 1; CREATE TABLE gwas_upload (
      id INTEGER PRIMARY KEY DEFAULT nextval('gwas_upload_id_sequence'),
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
      failure_reason TEXT,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )"
  ),
  study_extractions = list(
    name = "study_extractions",
    query = "CREATE SEQUENCE study_extractions_id_sequence START 1; CREATE TABLE study_extractions (
      id INTEGER PRIMARY KEY DEFAULT nextval('study_extractions_id_sequence'),
      gwas_upload_id INTEGER,
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
      ld_block TEXT
    )"
  ),
  coloc_pairs = list(
    name = "coloc_pairs",
    query = "CREATE TABLE coloc_pairs (
      gwas_upload_id INTEGER NOT NULL,
      existing_study_extraction_a BOOLEAN NOT NULL,
      study_extraction_a_id INTEGER NOT NULL,
      existing_study_extraction_b BOOLEAN NOT NULL,
      study_extraction_b_id INTEGER NOT NULL,
      ld_block_id INTEGER NOT NULL,
      h3 REAL CHECK (h3 BETWEEN 0 AND 1),
      h4 REAL CHECK (h4 BETWEEN 0 AND 1)
    )"
  ),
  coloc_groups = list(
    name = "coloc_groups",
    query = glue::glue("CREATE TABLE coloc_groups (
      gwas_upload_id INTEGER NOT NULL,
      coloc_group_id INTEGER NOT NULL,
      existing_study_extraction BOOLEAN NOT NULL,
      study_extraction_id INTEGER NOT NULL,
      snp_id INTEGER NOT NULL,
      ld_block_id INTEGER NOT NULL
    )")
  )
)