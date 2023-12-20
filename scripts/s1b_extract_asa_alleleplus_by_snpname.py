import os
import sys
import pandas as pd
from scripts.constants.raw_asa_dtype import raw_asa_dtype
from scripts.constants import config

def extractAsaAlleleplusBySnpname(idir, odir):
    ## Specify the path of mastervariant 
    mastervariant_path = config.mastervariant_file

    ## Import and read mastervariant as string
    try:
        mastervariant_df = pd.read_csv(mastervariant_path, sep='\t', header=0, skip_blank_lines=True).astype('string')
    except FileNotFoundError as e:
        print("File does not exist:", e, "\nPlease provide the correct mastervariant path in scripts/constants/config.py")
        sys.exit()

    ## Data cleaning for masterfile
    # Drop rows with all empty cells
    # Strip leading and trailing whitespaces and drop duplicates based on all columns
    mastervariant_df = mastervariant_df.dropna(
        axis=0,
        how='all',
        subset=None
    ).applymap(lambda x: x.strip()).drop_duplicates()
    
    ## Get a list of all the files in the directory
    asa_samples_list = os.listdir(idir)

    ## Loop every file in the directory into a dataframe
    for asa_sample in asa_samples_list:
        if asa_sample.endswith('.txt'):
            raw_asa_df = pd.read_csv(idir + asa_sample, sep='\t', header=0, dtype=raw_asa_dtype, skip_blank_lines=True)

            ## Merge dataframes by performing an inner merge based on SNP Name
            asa_mastervariant_merged_by_snpname_df = pd.merge(mastervariant_df, raw_asa_df[[config.asa_raw_snpname, config.asa_raw_sample_name, config.asa_raw_allele1plus, config.asa_raw_allele2plus]], how='inner', left_on= config.snpname, right_on=config.asa_raw_snpname).astype('string')

            ## Drop column ["Sample Name"] extracted from ASA raw file
            asa_mastervariant_merged_by_snpname_df = asa_mastervariant_merged_by_snpname_df.drop(columns=config.asa_raw_snpname)

            ## Rename columns
            rename_dict = {config.asa_raw_sample_name: config.asa_sample_name,
                        config.asa_raw_allele1plus: config.allele_1_plus,
                        config.asa_raw_allele2plus: config.allele_2_plus}
            
            asa_mastervariant_merged_by_snpname_df = asa_mastervariant_merged_by_snpname_df.rename(columns=rename_dict)

            ## Combine allele plus 1 & 2 to form ["asa_genotype"] 
            asa_mastervariant_merged_by_snpname_df[config.asa_genotype] = asa_mastervariant_merged_by_snpname_df[config.allele_1_plus] + asa_mastervariant_merged_by_snpname_df[config.allele_2_plus].astype('string')

            ## Drop 'fill missing values with '.'
            asa_mastervariant_merged_by_snpname_df = asa_mastervariant_merged_by_snpname_df.fillna('.')

            ## Output directory
            base_name = asa_sample.split('.')[0]
            output_file_name = f"{base_name}_mastervariant.txt"
            output_file_path = os.path.join(odir, output_file_name)
            asa_mastervariant_merged_by_snpname_df.to_csv(output_file_path, sep='\t', index=False)