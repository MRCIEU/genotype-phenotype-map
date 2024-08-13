source('constants.R')
parser <- argparser::arg_parser('Standardise alleles after extraction')
parser <- argparser::add_argument(parser, '--gwas_file', help = 'GWAS File', type = 'character')
args <- argparser::parse_args(parser)

gwas <- vroom::vroom(gwas_file, show_col_types = F)
gwas <- standardise_alleles(gwas)

vroom::vroom_write(gwas, gwas_file)