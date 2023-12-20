import os
import pandas as pd
import numpy as np
from scripts.constants import config

## Objective 1: For traits that have 0 TOTAL MARKER, it should only be output once.
## a) When ["MARKER_TOTAL"] is 0, replace ["GENE"] and ["ASA_GENOTYPE"] with '.'
## b) Deduplicate so for one trait, it will only output once

## Objective 2: For traits that have ["MARKER_TOTAL"] is > 0, only output variant that has MARKER detected. Don't need to output variants that have 0 MARKER detected

def filter_marker_row(df):
    duplicate_df = df.copy()

    # Replace ["GENE"] and ["ASA_GENOTYPE"] with '.' for rows where ["MARKER_TOTAL"] is 0
    duplicate_df.loc[duplicate_df[config.marker_total] == '0', [config.gene, config.asa_genotype]] = '.'

    # Remove rows where ["MARKER_TOTAL"] is > 0 and ["MARKER_NUM"] is 0
    # tilde(~) sign works as a NOT(!) operator
    return duplicate_df[~((duplicate_df[config.marker_total] > '0') & (duplicate_df[config.marker_num] == '0'))]

# Aggregate ["ASA_GENOTYPE"] values for each unique combination of ["SECTION_N_TRAIT"], ["MARKER_TOTAL"],  ["GAUGE"], ["GENOTYPIC_EFFECT"] and ["GENE"]
def aggregate_asa_genotype(df):
    return df.groupby([config.section_n_trait, config.marker_total, config.gauge, config.genotypic_effect, config.gene])[config.asa_genotype].apply(','.join).reset_index()

def extractFinalResults(idir, odir):
    ## Get a list of all the files in the directory
    asa_mastervariant_gt_grr_traitscr_marker_gauge_list = os.listdir(idir)
    
    ## Loop every file in the directory into a dataframe
    for asa_mastervariant_gt_grr_traitscr_marker_gauge in asa_mastervariant_gt_grr_traitscr_marker_gauge_list:
        if asa_mastervariant_gt_grr_traitscr_marker_gauge.endswith('.txt'):
            asa_mastervariant_gt_grr_traitscr_marker_gauge_df = pd.read_csv(idir + asa_mastervariant_gt_grr_traitscr_marker_gauge, sep='\t', header=0, skip_blank_lines=True).astype('string')

            ## Strip leading and trailing whitespaces
            asa_mastervariant_gt_grr_traitscr_marker_gauge_df = asa_mastervariant_gt_grr_traitscr_marker_gauge_df.applymap(lambda x: x.strip()).drop_duplicates()

            ## Extract columns asa_mastervariant_gt_grr_traitscr_marker_gauge_df
            extracted_origene_results_df = asa_mastervariant_gt_grr_traitscr_marker_gauge_df[[config.section_n_trait, config.marker_num, config.marker_total, config.gauge, config.genotypic_effect, config.gene, config.asa_genotype]]

            ## fill missing values with '.'
            extracted_origene_results_df = extracted_origene_results_df.fillna('.')

            ## Call the functions
            filtered_df = filter_marker_row(extracted_origene_results_df).drop_duplicates()

            final_df = aggregate_asa_genotype(filtered_df)

            ## Strip leading and trailing whitespaces
            final_df = final_df.applymap(lambda x: x.strip()).drop_duplicates().sort_values(by=[config.section_n_trait])

            ## Rearrange extracted columns
            final_df = final_df[[config.section_n_trait, config.marker_total, config.gauge, config.genotypic_effect, config.gene, config.asa_genotype]]

            ## DATA POST-PROCESSING

            # Replace blank values as NaN, then fill NaN with '.'
            final_df = final_df.replace(r'^\s*$', np.nan, regex=True).fillna('.')

            # Strip leading and tailing whitespaces
            final_df = final_df.applymap(lambda x: x.strip()).drop_duplicates()

            # Drop rows if whole row is '.'
            final_df = final_df.drop(final_df.index[(final_df == '.').all(axis=1)])

            ## Output
            base_name = asa_mastervariant_gt_grr_traitscr_marker_gauge[:asa_mastervariant_gt_grr_traitscr_marker_gauge.index('_mastervariant')]
            output_file_name = f"{base_name}_finalresults.txt"
            output_file_path = os.path.join(odir, output_file_name)
            final_df.to_csv(output_file_path, sep='\t', index=False)