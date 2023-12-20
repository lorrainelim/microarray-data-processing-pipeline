import os
import pandas as pd
from scripts.constants import config

## Define a function to Compare ["ALLELE_1_PLUS"] with ["REF"] and ["ALT"]
def compare_allele1_refalt(row):
    if row[config.allele_1_plus] == row[config.ref]:
        return '0'
    elif row[config.allele_1_plus] in row[config.alt].split(','):
        return str(row[config.alt].split(',').index(row[config.allele_1_plus]) + 1)
    elif row[config.allele_1_plus] == 'I':
        refalt_diff = len(row[config.ref]) - len(row[config.alt].split(',')[0])
        if refalt_diff > 0:
            return '0'
        else:
            return '1'
    elif row[config.allele_1_plus] == 'D':
        refalt_diff = len(row[config.ref]) - len(row[config.alt].split(',')[0])
        if refalt_diff > 0:
            return '1'
        else:
            return '0'
    else:
        return None

## Define a function to Compare ["ALLELE_2_PLUS"] with ["REF"] and ["ALT"]
def compare_allele2_refalt(row):
    if row[config.allele_2_plus] == row[config.ref]:
        return '0'
    elif row[config.allele_2_plus] in row[config.alt].split(','):
        return str(row[config.alt].split(',').index(row[config.allele_2_plus]) + 1)
    elif row[config.allele_2_plus] == 'I':
        refalt_diff = len(row[config.ref]) - len(row[config.alt].split(',')[0])
        if refalt_diff > 0:
            return '0'
        else:
            return '1'
    elif row[config.allele_2_plus] =='D':
        refalt_diff = len(row[config.ref]) - len(row[config.alt].split(',')[0])
        if refalt_diff > 0:
            return '1'
        else:
            return '0'
    else:
        return None

## Create ["ASA_GT"] column with sorted genotype in 0/1 format
def generate_genotype(row):
    allele_1_gt = row[config.allele_1_gt]
    allele_2_gt = row[config.allele_2_gt]

    if allele_1_gt is not None and allele_2_gt is not None:
        combined_genotype = sorted([int(allele_1_gt), int(allele_2_gt)])
        return f"{combined_genotype[0]}/{combined_genotype[1]}"
    else:
        # Handle the case when either allele is None
        combined_genotype = None  # or any other appropriate value
        return '.'
    

def generateVcfFormat(idir, odir):
    ## Get a list of all the files in the directory
    asa_mastervariant_list = os.listdir(idir)

    ## Loop every file in the directory into a dataframe
    for asa_mastervariant in asa_mastervariant_list:
        if asa_mastervariant.endswith('.txt'):
            asa_mastervariant_df = pd.read_csv(idir + asa_mastervariant, sep='\t', header=0, skip_blank_lines=True).astype('string')

            ## Data cleaning for mastervariant
            # Drop rows with all empty cells
            asa_mastervariant_df = asa_mastervariant_df.dropna(
                axis=0,
                how='all',
                subset=None
            ).applymap(lambda x: x.strip()).drop_duplicates()

            ## Create ["Allele1_GT"] column
            asa_mastervariant_df[config.allele_1_gt] = asa_mastervariant_df.apply(lambda row: compare_allele1_refalt(row), axis=1)

            ## Create ["Allele2_GT"] column
            asa_mastervariant_df[config.allele_2_gt] = asa_mastervariant_df.apply(lambda row: compare_allele2_refalt(row), axis=1)

            ## Create ["ASA_GT"] column
            asa_mastervariant_df[config.asa_gt] = asa_mastervariant_df.apply(generate_genotype, axis=1)
            
            ## fill missing values with '.'
            asa_mastervariant_df = asa_mastervariant_df.fillna('.')

            # ## Drop columns: Allele1_GT and Allele2_GT (if not required)
            # asa_mastervariant_df = asa_mastervariant_df.drop([config.allele_1_gt, config.allele_2_gt], axis=1)
            
            ## Output
            base_name = asa_mastervariant.split('.')[0]
            output_file_name = f"{base_name}_gt.txt"
            output_file_path = os.path.join(odir, output_file_name)
            asa_mastervariant_df.to_csv(output_file_path, sep='\t', index=False)