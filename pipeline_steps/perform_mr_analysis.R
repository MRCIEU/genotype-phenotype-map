source('constants.R')

#Interesting example: ukb-b-7953_10_12300790

#could be really easy if we store the chosen SNPs from hyprcoloc to run the MR on

parser <- argparser::arg_parser('Perform MR on coloc results')
parser <- argparser::add_argument(parser, '--studies_processed', help = 'Data of all studies processed', type = 'character')
parser <- argparser::add_argument(parser, '--all_study_regions', help = 'All study regions', type = 'character')
parser <- argparser::add_argument(parser, '--coloc_result_file', help = 'Coloc result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--mr_results', help = 'Coloc result file to save', type = 'character')
args <- argparser::parse_args(parser)

