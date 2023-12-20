import os
import pandas as pd
from scripts.constants import config

## A function to get/identify sample name from filename
def get_sample_name(filename):
    basename = os.path.basename(filename)
    sample_name = basename.split('_')[1]
    return sample_name

## A function to generate relative risk based on ["RISK_CODE"] and ["ASA_GT"] 
def generate_relative_risk(df, population):
    effect_allele = df[config.risk_code]
    asa_gt = df[config.asa_gt]
    r_sq = df[config.odd_ratio] ** 2
    p_sq_r_sq = (df[config.odd_ratio] ** 2) * (df[population] ** 2)
    pqr2 = df[population] * (1 - df[population])*df[config.odd_ratio] * 2
    p_sq = df[population] ** 2
    q_sq = (1 - df[population]) ** 2
    avg_pop_risk_1_effect = p_sq_r_sq + pqr2 + q_sq
    avg_pop_risk_0_effect = (q_sq * r_sq) + pqr2 + p_sq

    if pd.Series(effect_allele).eq('1').any():
    # if effect_allele == '1':
        if len(asa_gt) == 3:
            if asa_gt == '0/0':
                relative_risk_gt00_nc = round(1 / avg_pop_risk_1_effect, 4)
                return relative_risk_gt00_nc
            elif asa_gt.startswith('0/'):
                relative_risk_gt01 = round(df[config.odd_ratio] / avg_pop_risk_1_effect, 4)
                return relative_risk_gt01
            elif asa_gt[0] != '0' and asa_gt[2] != '0':
                relative_risk_gt11_c = round(r_sq / avg_pop_risk_1_effect, 4)
                return relative_risk_gt11_c
            else:
                return None
        elif asa_gt == '.':
                return '1'
        else:
            return None  
        
    elif pd.Series(effect_allele).eq('0').any():
    # elif effect_allele == '0':
        if len(asa_gt) == 3:        
            if asa_gt == '0/0':
                relative_risk_gt00_c = round(r_sq / avg_pop_risk_0_effect, 4)
                return relative_risk_gt00_c
            elif asa_gt.startswith('0/'):
                relative_risk_gt01 = round(df[config.odd_ratio] / avg_pop_risk_0_effect, 4)
                return relative_risk_gt01
            elif asa_gt[0] != '0' and asa_gt[2] != '0':
                relative_risk_gt11_nc = round(1 / avg_pop_risk_0_effect, 4)
                return relative_risk_gt11_nc
            else:
                return None
        elif asa_gt == '.':
                return '1'        
        else:
            return None
    else:
        return None

def generateGrr(popfile, idir, odir):
    ## Load the population dataframe if popfile argument is parsed
    if popfile is not None: # File exists
        popfile_df = pd.read_csv(popfile, sep='\t', header=0)
        
        # Clean popfile_df["sample_name"]
        popfile_df[config.pop_sample_name] = popfile_df[config.pop_sample_name].str.strip().str.extract(r'(FM\d+)')

        popfile_dict = dict(zip(popfile_df[config.pop_sample_name], popfile_df[config.pop_population]))
    else:
        popfile_dict = {}

    ## Get a list of all the files from the input directory
    asa_mastervariant_gt_list = os.listdir(idir)

    ## Loop every file in the directory into a dataframe
    for filename in asa_mastervariant_gt_list:
        if filename.endswith('.txt'):
            asa_mastervariant_gt_df = pd.read_csv(idir + filename, sep='\t', header=0, skip_blank_lines=True).astype('string')
            sample_name = get_sample_name(filename)
            if sample_name in popfile_dict.keys():
                population = popfile_dict[sample_name]
            else: print("Error: Sample Name from filename does not match with any of the Sample Names in sample population file")

            ## Strip leading and trailing whitespaces and drop duplicates based on all columns
            asa_mastervariant_gt_df.applymap(lambda x: x.strip()).drop_duplicates()

            ## Convert data type string to float
            asa_mastervariant_gt_df[config.odd_ratio] = asa_mastervariant_gt_df[config.odd_ratio].astype(float)
            asa_mastervariant_gt_df[population] = asa_mastervariant_gt_df[population].astype(float)

            ## Create ["GRR"] column for the results from generate_relative_risk function
            asa_mastervariant_gt_df[config.grr] = asa_mastervariant_gt_df.apply(generate_relative_risk, population=population, axis=1)
            
            ## fill missing values with '.'
            asa_mastervariant_gt_df = asa_mastervariant_gt_df.fillna('.')

            ## Output
            base_name = filename.split('.')[0]
            output_file_name = f"{base_name}_grr.txt"
            output_file_path = os.path.join(odir, output_file_name)
            asa_mastervariant_gt_df.to_csv(output_file_path, sep='\t', index=False)