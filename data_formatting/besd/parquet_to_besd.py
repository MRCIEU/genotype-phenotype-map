import click
import pandas as pd
import numpy as np
import subprocess
import tempfile

@click.command(name='Convert Parquet file into besd format')
@click.option('--input_parquet_file', help='Parquet File to read', required=True)
@click.option('--output_besd_prefix', help='BESD File to save', required=True)
def main(input_parquet_file, output_besd_prefix):
    csv_file = tempfile.TemporaryFile()

    convert_parquet_to_csv(input_parquet_file, csv_file)
    run_smr_to_besd(csv_file, output_besd_prefix)

# Function to convert Parquet to CSV (intermediate format for SMR)
def convert_parquet_to_csv(parquet_file, csv_file):
    try:
        df = pd.read_parquet(parquet_file)
        print(len(df.index))
        print(df.head(10))

        
        one_gene = np.where(df['gene_id'] == 'ENSG00000227232.7')
        print(len(df.iloc[one_gene]))
        exit()
        # Assuming the necessary columns are present, rename columns as needed by SMR.
        # Update 'SNP', 'effect_size', etc., based on the actual column names in your data.
        df.rename(columns={
            'snp_column': 'SNP',
            'effect_size_column': 'Effect_Size',
            'p_value_column': 'P_Value',
            'other_needed_column': 'Other'
        }, inplace=True)

        # Save as CSV
        df.to_csv(csv_file, sep='\t', index=False)
        print(f"Parquet data converted to CSV at: {csv_file}")
    except Exception as e:
        print(f"Error converting Parquet to CSV: {e}")

# Function to run SMR tool to generate BESD file
def run_smr_to_besd(csv_file, output_prefix):
    try:
        # Command to invoke the SMR tool for creating a BESD file
        smr_command = [
            "smr",  # Adjust the path if smr is not in PATH
            "--make-besd",
            "--gwas-summary", csv_file,
            "--out", output_prefix
        ]
        
        # Execute the command
        subprocess.run(smr_command, check=True)
        print(f"BESD file generated with prefix: {output_prefix}")
    except subprocess.CalledProcessError as e:
        print(f"SMR tool failed: {e}")
    except FileNotFoundError:
        print("SMR tool not found. Ensure it is installed and in your PATH.")


if __name__ == "__main__":
    main()

