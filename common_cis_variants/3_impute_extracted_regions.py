import click
import pandas as pd
import raiss.ld_matrix
import scipy
from raiss import ImputationLauncher

click.command()
click.option('--ld_matrix', help='LD matrix file')
click.option('--chr', help='Chromosome')
click.option('--extracted_region_file', help='Extracted region file')
click.option('--output_file', help='output file of imputed summary statistics')


def main(ld_matrix, chr, extracted_region, output_file):

    mat_ld = raiss.ld_matrix.load_plink_ld(ld_matrix)

    #stolen from pipes.save_chromosome_imputation
    print("Imputation of {0} gwas for chromosome {1}".format(gwas, chrom))
    # Imputer settings
    imputer = ImputationLauncher( window_size=int(window_size), buf=int(buffer_size), lamb= float(l2_regularization), pinv_rtol = float(eigen_threshold), ld_type = ld_type)
    # Reading of inputs
    z_file = "{0}/z_{1}_{2}.txt".format(zscore_folder, gwas, chrom)
    zscore = pd.read_csv(z_file, index_col=0, sep="\t")
    ref_panel_file = ref_folder + "/"+ ref_panel_prefix + chrom + ref_panel_suffix
    ref_panel = pd.read_csv(ref_panel_file, sep="\t", names=['chr', "nothing", 'pos', 'Ref_all', 'alt_all'], index_col = 1)

    # imputation
    imputed_zscore = imputer.chromosome_imputation(chrom, zscore, ref_panel, ld_folder)
    print("Imputation DONE")

    # Formatting and filtering
    # and Saving results
    z_fo = "{0}/z_{1}_{2}{3}.txt".format(output_folder, gwas, chrom, tag)
    filter_output(imputed_zscore, z_fo, float(R2_threshold), minimum_ld=float(minimum_ld))
    print("Save imputation done at {0}".format(z_fo))

    scipy_sparse_file = 1
    #actually, dont use this, create it yourself, but use it as a way to know how to format it
    raiss.LD.load_sparse_matrix()
    raiss.LD.generate_genome_matrices()
    raiss.pipes.save_chromosome_imputation()


if __name__ == '__main__':
    main()