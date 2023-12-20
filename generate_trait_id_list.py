import pandas as pd
import os

## Specify the path of masterfile 
file_dir = os.path.dirname(__file__)
mastervariant_path = os.path.join(file_dir,"masterfile_sample.txt")

## Import and read masterfile as string
mastervariant_df = pd.read_csv(mastervariant_path, sep='\t', header=0, skip_blank_lines=True).astype('string')

## Data cleaning for masterfile
# Drop rows with all empty cells
mastervariant_df.dropna(
    axis=0,
    how='all',
    subset=None
).applymap(lambda x: x.strip())

# Drop duplicates based on all columns
mastervariant_df.drop_duplicates()

# Create a new dataframe (A_df) by extracting 3 columns from another dataframe (B_df)
trait_id_df = mastervariant_df[['ID', 'TRAIT_L', 'TRAIT_S']].copy().drop_duplicates()

# Output
trait_id_df.to_csv('trait_id_list.txt', sep='\t',index=False)

