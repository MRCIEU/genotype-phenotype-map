from __future__ import annotations
from typing import Any

import click
import glob
import os
import numpy as np
import pandas as pd
import scipy.linalg

DATA_DIR = os.getenv("DATA_DIR")
imputed_r2_threshold = 0.9
ld_score_threshold = 5

@click.command(name="Impute studies by region")
@click.option('--ld_region_prefix', help='LD region used for imputation', required=True)
@click.option('--ld_block_dir', help='List of GWAS Summary Statistics file', required=True)
def main(ld_region_prefix, ld_block_dir):
    ld_matrix = pd.read_csv(ld_region_prefix + '.ld', header=None, delimiter=' ')
    ld_matrix = np.array(ld_matrix)

    ld_region_from_reference_panel = pd.read_csv(ld_region_prefix + '.tsv', delimiter='\t')
    extracted_studies = pd.read_csv(ld_block_dir + 'extracted_studies.tsv', delimiter='\t')

    for extracted_study in extracted_studies.iterrows():
        gwas_file = extracted_study['file']
        imputed_file = gwas_file.replace('original', 'imputed')
        if os.path.isfile(imputed_file):
            print('Imputed file exists, skipping.')
            continue

        gwas = pd.read_csv(gwas_file, delimiter='\t')
        gwas.drop_duplicates(subset=['RSID'], inplace=True)
        gwas['Z'] = gwas.apply(lambda row: row.BETA / row.SE, axis=1)

        rsids_in_gwas = np.isin(ld_region_from_reference_panel.RSID, gwas.RSID)
        rsids_in_ld_block = np.isin(gwas.RSID, ld_region_from_reference_panel.RSID)
        gwas_filter = np.where(rsids_in_ld_block)[0]
        known = np.where(rsids_in_gwas)[0]
        unknown = np.where(rsids_in_gwas == False)[0]

        sig_t = ld_matrix[known, :][:, known]
        sig_i_t = ld_matrix[unknown, :][:, known]
        z = np.array(gwas.iloc[:, -1])
        z = z[gwas_filter]

        imputation_results = SummaryStatisticsImputation.raiss_model(
            z, sig_t, sig_i_t, lamb=0.01, rtol=0.01
        )

        rsids_to_add = (imputation_results["imputation_r2"] >= imputed_r2_threshold) * (
                imputation_results["ld_score"] >= ld_score_threshold
        )
        if sum(rsids_to_add) >= 1:
            ld_region_from_reference_panel = ld_region_from_reference_panel.iloc[unknown]
            ld_region_from_reference_panel['Z'] = imputation_results['mu']
            imputed_gwas_data_to_add = ld_region_from_reference_panel[rsids_to_add]
            gwas = pd.concat([gwas,imputed_gwas_data_to_add], axis=0, ignore_index=True)

        print(f'Imputed {sum(rsids_to_add)} SNPs')
        gwas.to_csv(imputed_file, sep='\t')
        os.symlink(imputed_file, ld_block_dir + 'imputed/', target_is_directory=True)
        #TODO: update extracted_studies.tsv with imputed file, and how many rows were imputed / added
        #TODO: OR create a new file, which is what snakemake could look at?


class SummaryStatisticsImputation:
    """Implementation of RAISS summary statstics imputation model."""

    @staticmethod
    def raiss_model(
            z_scores_known: np.ndarray,
            ld_matrix_known: np.ndarray,
            ld_matrix_known_missing: np.ndarray,
            lamb: float = 0.01,
            rtol: float = 0.01,
    ) -> dict[str, Any]:
        """Compute the imputation of the z-score using the RAISS model.

        Args:
            z_scores_known (np.ndarray): the vector of known Z scores
            ld_matrix_known (np.ndarray) : the matrix of known LD correlations
            ld_matrix_known_missing (np.ndarray): LD matrix of known SNPs with other unknown SNPs in large matrix (similar to ld[unknowns, :][:,known])
            lamb (float): size of the small value added to the diagonal of the covariance matrix before inversion. Defaults to 0.01.
            rtol (float): threshold to filter eigenvectos by its eigenvalue. It makes an inversion biased but much more numerically robust. Default to 0.01.

        Returns:
            dict[str, Any]:
                - var (np.ndarray): variance of the imputed SNPs
                - mu (np.ndarray): the estimation of the zscore of the imputed SNPs
                - ld_score (np.ndarray): the linkage disequilibrium score of the imputed SNPs
                - condition_number (np.ndarray): the condition number of the correlation matrix
                - correct_inversion (np.ndarray): a boolean array indicating if the inversion was successful
                - imputation_r2 (np.ndarray): the R2 of the imputation
        """
        sig_t_inv = SummaryStatisticsImputation._invert_sig_t(
            ld_matrix_known, lamb, rtol
        )
        if sig_t_inv is None:
            return {
                "var": None,
                "mu": None,
                "ld_score": None,
                "condition_number": None,
                "correct_inversion": None,
                "imputation_r2": None,
            }
        else:
            condition_number = np.array(
                [np.linalg.cond(ld_matrix_known)] * ld_matrix_known_missing.shape[0]
            )
            correct_inversion = np.array(
                [
                    SummaryStatisticsImputation._check_inversion(
                        ld_matrix_known, sig_t_inv
                    )
                ]
                * ld_matrix_known_missing.shape[0]
            )

            var, ld_score = SummaryStatisticsImputation._compute_var(
                ld_matrix_known_missing, sig_t_inv, lamb
            )

            mu = SummaryStatisticsImputation._compute_mu(
                ld_matrix_known_missing, sig_t_inv, z_scores_known
            )
            var_norm = SummaryStatisticsImputation._var_in_boundaries(var, lamb)

            R2 = (1 + lamb) - var_norm

            mu = mu / np.sqrt(R2)
            return {
                "var": var,
                "mu": mu,
                "ld_score": ld_score,
                "condition_number": condition_number,
                "correct_inversion": correct_inversion,
                "imputation_r2": 1 - var,
            }

    @staticmethod
    def _compute_mu(
            sig_i_t: np.ndarray, sig_t_inv: np.ndarray, zt: np.ndarray
    ) -> np.ndarray:
        """Compute the estimation of z-score from neighborring snp.

        Args:
            sig_i_t (np.ndarray) : correlation matrix with line corresponding to unknown Snp (snp to impute) and column to known SNPs
            sig_t_inv (np.ndarray): inverse of the correlation matrix of known matrix
            zt (np.ndarray): Zscores of known snp
        Returns:
            np.ndarray: a vector of length i containing the estimate of zscore

        """
        return np.dot(sig_i_t, np.dot(sig_t_inv, zt))

    @staticmethod
    def _compute_var(
            sig_i_t: np.ndarray, sig_t_inv: np.ndarray, lamb: float
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compute the expected variance of the imputed SNPs.

        Args:
            sig_i_t (np.ndarray) : correlation matrix with line corresponding to unknown Snp (snp to impute) and column to known SNPs
            sig_t_inv (np.ndarray): inverse of the correlation matrix of known matrix
            lamb (float): regularization term added to matrix

        Returns:
            tuple[np.ndarray, np.ndarray]: a tuple containing the variance and the ld score
        """
        var = (1 + lamb) - np.einsum(
            "ij,jk,ki->i", sig_i_t, sig_t_inv, sig_i_t.transpose()
        )
        ld_score = (sig_i_t ** 2).sum(1)

        return var, ld_score

    @staticmethod
    def _check_inversion(sig_t: np.ndarray, sig_t_inv: np.ndarray) -> bool:
        """Check if the inversion is correct.

        Args:
            sig_t (np.ndarray): the correlation matrix
            sig_t_inv (np.ndarray): the inverse of the correlation matrix
        Returns:
            bool: True if the inversion is correct, False otherwise
        """
        return np.allclose(sig_t, np.dot(sig_t, np.dot(sig_t_inv, sig_t)))

    @staticmethod
    def _var_in_boundaries(var: np.ndarray, lamb: float) -> np.ndarray:
        """Forces the variance to be in the 0 to 1+lambda boundary. Theoritically we shouldn't have to do that.

        Args:
            var (np.ndarray): the variance of the imputed SNPs
            lamb (float): regularization term added to the diagonal of the sig_t matrix

        Returns:
            np.ndarray: the variance of the imputed SNPs
        """
        id_neg = np.where(var < 0)
        var[id_neg] = 0
        id_inf = np.where(var > (0.99999 + lamb))
        var[id_inf] = 1

        return var

    @staticmethod
    def _invert_sig_t(sig_t: np.ndarray, lamb: float, rtol: float) -> np.ndarray:
        """Invert the correlation matrix. If the provided regularization values are not enough to stabilize the inversion process for the given matrix, the function calls itself recursively, increasing lamb and rtol by 10%.

        Args:
            sig_t (np.ndarray): the correlation matrix
            lamb (float): regularization term added to the diagonal of the sig_t matrix
            rtol (float): threshold to filter eigenvector with a eigenvalue under rtol make inversion biased but much more numerically robust

        Returns:
            np.ndarray: the inverse of the correlation matrix
        """
        try:
            np.fill_diagonal(sig_t, (1 + lamb))
            print(sig_t.dtype)
            sig_t_inv = scipy.linalg.pinv(sig_t, rtol=rtol, atol=0)
            return sig_t_inv
        except np.linalg.LinAlgError:
            return SummaryStatisticsImputation._invert_sig_t(
                sig_t, lamb * 1.1, rtol * 1.1
            )


if __name__ == "__main__":
    main()
