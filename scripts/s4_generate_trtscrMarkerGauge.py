import os
import pandas as pd
import numpy as np
from scripts.constants import config

### OBJECTIVE: Generate Trait Score, Marker, Total Marker, and assigning Gauge and Genotypic Effect

## Calculate number of markers in asa_mastervariant_df
def count_markers(row):
    effect_allele = row[config.risk_code]
    asa_gt = row[config.asa_gt]
    if effect_allele == '1':
        if asa_gt == '0/0':
            num_of_marker = '0'
            return num_of_marker 
        elif asa_gt.startswith('0/'):
            num_of_marker = '1'
            return num_of_marker
        elif len(asa_gt) >= 3 and asa_gt[0] != '0' and asa_gt[2] != '0':
            num_of_marker = '2'
            return num_of_marker
        elif asa_gt == '.':
            num_of_marker = '0'
            return num_of_marker
        else:
            return None  
    elif effect_allele == '0':
        if asa_gt == '0/0':
            num_of_marker = '2'
            return num_of_marker
        elif asa_gt.startswith('0/'):
            num_of_marker = '1'
            return num_of_marker
        elif len(asa_gt) >= 3 and asa_gt[0] != '0' and asa_gt[2] != '0':
            num_of_marker = '0'
            return num_of_marker
        elif asa_gt == '.':
            num_of_marker = '0'
            return num_of_marker
        else: 
            return None
    else:
        return None

## Gauge Ranges are below:
# Gauge 1 == {0.9 - 1.09}
# Gauge 2 == {1.10 - 1.29}
# Gauge 3 == {> 1.30}
# Gauge 4 == {< 0.89}

## Assign gauge based on ["TRAIT_SCORE"]
def assign_gauge(row):
    trait_score = row[config.trait_score]
    if 0.9 <= trait_score <= 1.09:
        gauge = '1'
    elif 1.10 <= trait_score <= 1.29:
        gauge = '2'
    elif trait_score > 1.30:
        gauge = '3'
    elif trait_score < 0.89:
        gauge = '4'
    else:
        gauge = None
    return gauge

## Extract gauge measurement based on ["GAUGE"]
def extract_genotypic_effect(row):
    gauge = row[config.gauge]
    genotypic_effect = ''
    if gauge == '1':
        genotypic_effect = row[config.ge1]
    elif gauge == '2':
        genotypic_effect = row[config.ge2]
    elif gauge == '3':
        genotypic_effect = row[config.ge3]
    elif gauge == '4':
        genotypic_effect = row[config.ge4]
    else: 
        gauge = None
    return genotypic_effect

def generateTrtscrMarkerGauge(traitid_file, idir, odir):
    ## Load trait id table as dataframe
    trait_id_df = pd.read_csv(traitid_file, sep='\t', header=0, skip_blank_lines=True).astype('string')

    ## Strip leading and trailing whitespaces
    trait_id_df.applymap(lambda x: x.strip()).drop_duplicates()

    ## Get a list of all the files in the directory
    asa_mastervariant_list = os.listdir(idir)

    ## Loop every file in the directory into a dataframe
    for asa_mastervariant in asa_mastervariant_list:
        if asa_mastervariant.endswith('.txt'):
            ### DF1 (asa mastervariant columns) = asa_mastervariant_df
            asa_mastervariant_df = pd.read_csv(idir + asa_mastervariant, sep='\t', header=0, skip_blank_lines=True).astype('string')

            ## Strip leading and trailing whitespaces
            asa_mastervariant_df.applymap(lambda x: x.strip()).drop_duplicates()

            ## Create new column for number of markers in asa_mastervariant_df
            asa_mastervariant_df[config.marker_num] = asa_mastervariant_df.apply(lambda row: count_markers(row), axis=1)

            ### DF2 (asa mastervariant columns with same SNPNAME with trait ID file) = asa_mastervariant_in_traitID_list_df 
            ## DF2: Merged asa_mastervariant + trait id table

            asa_mastervariant_in_traitID_list_df = pd.merge(trait_id_df[[config.snpname, config.trait_id, config.section, config.section_n_trait, config.assoc_type, config.ge1, config.ge2, config.ge3, config.ge4]], asa_mastervariant_df, on=[config.snpname, config.section_n_trait, config.assoc_type], how='inner')

            # Convert data type ["GRR"] and ["MARKER_NUM"] to float
            asa_mastervariant_in_traitID_list_df[config.grr] = asa_mastervariant_in_traitID_list_df[config.grr].astype(float)
            asa_mastervariant_in_traitID_list_df[config.marker_num] = asa_mastervariant_in_traitID_list_df[config.marker_num].astype(int)
            
            ### GROUPBY (columns: trait_id, trait_score, marker_total, gauge)
            ## Groupby trait_id from asa_mastervariant_in_traitID_list_df
            # calculate the product of the ["GRR"] column in trait_score_df
            # calculate the sum of the ["MARKER_NUM"] column in trait_score_df
            
            traitScore_markerTotal_gauge_groupby = asa_mastervariant_in_traitID_list_df.groupby(config.trait_id).agg({config.grr: 'prod', config.marker_num: 'sum'}).reset_index()

            # Rename ["GRR"] to ["TRAIT_SCORE"] and round up to 4 decimal points
            traitScore_markerTotal_gauge_groupby = traitScore_markerTotal_gauge_groupby.rename(columns={config.grr: config.trait_score}).round(4)
            # Rename ["MARKER_NUM"] to ["MARKER_TOTAL"]
            traitScore_markerTotal_gauge_groupby = traitScore_markerTotal_gauge_groupby.rename(columns={config.marker_num: config.marker_total})

            ## Convert data type of ["TRAIT_SCORE"] to float
            traitScore_markerTotal_gauge_groupby[config.trait_score] = traitScore_markerTotal_gauge_groupby[config.trait_score].astype(float)
            
            ## Call assign_gauge function based on ["TRAIT_SCORE"] and return results in new column ["GAUGE"]
            traitScore_markerTotal_gauge_groupby[config.gauge] = traitScore_markerTotal_gauge_groupby.apply(assign_gauge, axis=1)

            ## -- GROUPBY ended --

            asa_mastervariant_sec_merge_df = pd.merge(asa_mastervariant_in_traitID_list_df, traitScore_markerTotal_gauge_groupby, on=[config.trait_id], how='inner')

            ## Call extract_genotypic_effect function and return results in new column ["GENOTYPIC_EFFECT"]
            asa_mastervariant_sec_merge_df[config.genotypic_effect] = asa_mastervariant_sec_merge_df.apply(extract_genotypic_effect, axis=1)

            ## DATA POST-PROCESSING

            # Replace blank values as NaN, then fill NaN with '.'
            asa_mastervariant_sec_merge_df = asa_mastervariant_sec_merge_df.replace(r'^\s*$', np.nan, regex=True).fillna('.').astype('string')

            # Strip leading and tailing whitespaces
            asa_mastervariant_sec_merge_df = asa_mastervariant_sec_merge_df.applymap(lambda x: x.strip()).drop_duplicates()

            # Drop rows if whole row is '.'
            asa_mastervariant_sec_merge_df = asa_mastervariant_sec_merge_df.drop(asa_mastervariant_sec_merge_df.index[(asa_mastervariant_sec_merge_df == '.').all(axis=1)])

            ## Rearrange columns
            asa_mastervariant_sec_merge_df = asa_mastervariant_sec_merge_df[[config.snpname, config.trait_id, config.section, config.section_n_trait, config.assoc_type, config.ge1, config.ge2, config.ge3, config.ge4, config.gene, config.no_chr, config.chrom, config.pos37, config.pos38, config.rsid, config.ref, config.alt, config.risk_allele, config.risk_code, config.allele_type, config.odd_ratio, config.afr, config.amr, config.eas, config.eur, config.sas, config.asa_sample_name, config.allele_1_plus, config.allele_2_plus, config.asa_genotype, config.allele_1_gt, config.allele_2_gt, config.asa_gt, config.grr, config.trait_score, config.marker_num, config.marker_total, config.gauge, config.genotypic_effect]]

            base_name = asa_mastervariant.split('.')[0]
            output_file_name = f"{base_name}_trtscrMarkerGauge.txt"
            output_file_path = os.path.join(odir, output_file_name)
            asa_mastervariant_sec_merge_df.to_csv(output_file_path, sep='\t', index=False)