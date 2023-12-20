import os
import pandas as pd
from scripts.constants.raw_asa_dtype import raw_asa_dtype
from scripts.constants import config

def splitAsaBySampleName(asa_raw_file, popfile, odir):
    ## ASA run ID
    asa_run_id = os.path.basename(asa_raw_file).split('_')[1]

    ## Load the raw asa data file
    asa_raw_df = pd.read_csv(asa_raw_file, sep='\t', header=9, dtype=raw_asa_dtype, skip_blank_lines=True)

    ## Clean asa_raw_df["Sample Name"]
    asa_raw_df[config.asa_raw_sample_name] = asa_raw_df[config.asa_raw_sample_name].str.strip().str.extract(r'(FM\d+)')
    
    ## Split ASA raw file based on the sample name listed in popfile
    popfile_df = pd.read_csv(popfile, sep='\t', header=0)

    # Clean popfile_df["sample_name"]
    popfile_df[config.pop_sample_name] = popfile_df[config.pop_sample_name].str.strip().str.extract(r'(FM\d+)')

    popfile_dict = dict(zip(popfile_df[config.pop_sample_name], popfile_df[config.pop_population]))
    for sample_name in popfile_dict.keys(): 
        mask = asa_raw_df[config.asa_raw_sample_name] == sample_name
        if mask.any():
            # Select the rows that match the sample name and output them
            matching_rows = asa_raw_df[mask]
            ## Group the data by sample name
            grouped_by_samplename = matching_rows.groupby(config.asa_raw_sample_name)
            for sample_name, sample_data in grouped_by_samplename:
                ## Output directory
                output_file_name = f"{asa_run_id}_{sample_name}.txt"
                output_file_path = os.path.join(odir, output_file_name)
                sample_data.to_csv(output_file_path, sep='\t', index=False)
        else:
            None
