source('constants.R')
source('svg_helpers.R')

parser <- argparser::arg_parser('Create static web files')
parser <- argparser::add_argument(parser, '--static_web_dir', help = 'Directory for static web files (robots.txt, sitemap.xml)', type = 'character')
parser <- argparser::add_argument(parser, '--studies_db_file', help = 'Studies DB file', type = 'character')
parser <- argparser::add_argument(parser, '--static_web_files_ready_file', help = 'Static web files ready file', type = 'character')

args <- argparser::parse_args(parser)

main <- function() {
  dir.create(args$static_web_dir, recursive = TRUE, showWarnings = FALSE)

  create_opengwas_map(args$studies_db_file, args$static_web_dir)
  create_seo_files(args$studies_db_file, args$static_web_dir)
  prepare_svg_files_for_use(args$studies_db_file)

  rmarkdown::render('pipeline_summary.Rmd', output_file = file.path(args$static_web_dir, 'pipeline_summary.html'))
  vroom::vroom_write(data.frame(), args$static_web_files_ready_file)
}

create_opengwas_map <- function(studies_db_file, static_web_dir) {
  studies_conn <- duckdb::dbConnect(duckdb::duckdb(), studies_db_file)

  phenotype_studies <- DBI::dbGetQuery(studies_conn, "SELECT trait_id, trait_name, study_name FROM studies
    JOIN study_sources ON studies.source_id = study_sources.id
    JOIN traits ON studies.trait_id = traits.id
    WHERE source IN ('ebi_catalog', 'ukb')
    AND studies.data_type = 'phenotype' AND studies.variant_type = 'common'"
  )
  jsonlite::write_json(phenotype_studies$study_name, file.path(static_web_dir, 'opengwas_ids.json'), auto_unbox = TRUE, pretty = TRUE)
  phenotype_id_map <- setNames(as.list(phenotype_studies$trait_name), phenotype_studies$trait_id)
  jsonlite::write_json(phenotype_id_map, file.path(static_web_dir, 'phenotype_id_map.json'), auto_unbox = TRUE, pretty = TRUE)
  duckdb::dbDisconnect(studies_conn, shutdown = TRUE)
}

create_seo_files <- function(studies_db_file, static_web_dir) {
  studies_conn <- duckdb::dbConnect(duckdb::duckdb(), studies_db_file)

  studies <- DBI::dbGetQuery(studies_conn, "SELECT * FROM studies
    WHERE data_type = 'phenotype' AND variant_type = 'common'"
  )
  
  genes <- DBI::dbGetQuery(studies_conn, "
    SELECT gene, COUNT(*) as count 
    FROM coloc_groups_wide 
    WHERE gene IS NOT NULL 
    GROUP BY gene 
    ORDER BY count DESC 
    LIMIT 10000
  ")

  base_url <- gpm_website_data$url
  
  sitemap_file <- file.path(static_web_dir, "sitemap.xml")
  
  sitemap_header <- '<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">'
  sitemap_footer <- '</urlset>'
  
  static_pages <- c("", "/about", "/contact", "/faq", "/data")

  # Static pages
  static_urls <- glue::glue_data(
    list(page = static_pages),
    "
  <url>
    <loc>{base_url}{page}</loc>
    <changefreq>monthly</changefreq>
    <priority>0.8</priority>
  </url>"
  ) |> paste(collapse = "")
  
  # Study pages
  study_urls <- glue::glue_data(
    studies,
    "
  <url>
    <loc>{base_url}/traits.html?id={trait_id}</loc>
    <changefreq>monthly</changefreq>
    <priority>0.6</priority>
  </url>"
  ) |> paste(collapse = "")
  
  # Gene pages
  gene_urls <- glue::glue_data(
    genes,
    "
  <url>
    <loc>{base_url}/gene.html?id={gene}</loc>
    <changefreq>monthly</changefreq>
    <priority>0.6</priority>
  </url>"
  ) |> paste(collapse = "")
  
  writeLines(glue::glue("{sitemap_header}{static_urls}{study_urls}{gene_urls}\n{sitemap_footer}"), sitemap_file)

  robots_file <- file.path(static_web_dir, "robots.txt")
  robots_content <- glue::glue("
    User-agent: *
    Disallow: /api/
    Crawl-delay: 10
    Disallow: /static/
    User-agent: GPTBot
    Disallow: /
    User-agent: CCBot
    Disallow: /
    Allow: /

    Sitemap: {base_url}/static/sitemap.xml
  ")
  writeLines(robots_content, robots_file)

  duckdb::dbDisconnect(studies_conn, shutdown = TRUE)
}

main()