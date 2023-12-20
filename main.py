import os
import glob
import logging
import argparse
from scripts.constants import config
from scripts.s1a_split_asa_by_sample_name import splitAsaBySampleName
from scripts.s1b_extract_asa_alleleplus_by_snpname import extractAsaAlleleplusBySnpname
from scripts.s2_generate_vcf_format import generateVcfFormat
from scripts.s3_generate_grr import generateGrr
from scripts.s4_generate_trtscrMarkerGauge import generateTrtscrMarkerGauge
from scripts.s5_extract_finalresults import extractFinalResults

### Can uncomment the following module and codes near to the end of the script to remove intermediate directories when necessary
# import shutil

## Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()

## Create an argument parser
parser = argparse.ArgumentParser()

## Define arguments
parser.add_argument('-asa_run_id', type=str, required=True, help='ASA run id, eg: 20XX391600XX')
parser.add_argument('-traitid_file', type=str, required=True, help='trait id file based on test')

## Parse the arguments
args = parser.parse_args()

## Extract ASA run id
asa_run_id = args.asa_run_id

## Declare input directory
main_input_dir = f"files/{asa_run_id}/input/"

## Extract ASA raw file
# Get the list of .txt files in the input directory
txt_files = glob.glob(main_input_dir + '*.txt')

# Filter the txt_files that contain the run id in the filename
asa_run_id_txt_files = [file for file in txt_files if asa_run_id in file]

# Assuming only one .txt file matches the criteria, extract ASA raw file (.txt)
if len(asa_run_id_txt_files) == 1:
    asa_raw_file = asa_run_id_txt_files[0]
    print('ASA raw file:', asa_raw_file)
else: print("Error: Please ensure only one text file (ASA raw file, .txt) in the input folder.")

## Extract sample population file (.tsv)
# Get the list of .tsv files in the input directory
tsv_files = glob.glob(main_input_dir + '*.tsv')
pop_tsv_files =  [file for file in tsv_files if asa_run_id in file and 'population' in file]

# Assuming only 1 .tsv file matches the criteria
if len(pop_tsv_files) == 1:
    popfile = pop_tsv_files[0]
    print("Sample population file:", popfile)
else: print("Error: Please ensure only one tab-separated values file (sample population file, .tsv) in the input folder.")

## Declare -traitid_file argument & mastervariant file
print("Trait id file:", args.traitid_file)
print("Mastervariant file:", config.mastervariant_file)

## Create the main output directory name
asa_main_output_dir = f"files/{asa_run_id}/output/"

## Create main output directory if it doesn't exist
if not os.path.exists(asa_main_output_dir):
    os.mkdir(asa_main_output_dir)

## Define input and output subdirectories 
script_1a_output_dir = os.path.join(asa_main_output_dir, "s1a_asa_sample/")

script_1b_input_dir = script_1a_output_dir
script_1b_output_dir = os.path.join(asa_main_output_dir, "s1b_asa_mastervariant/")

script_2_input_dir = script_1b_output_dir
script_2_output_dir = os.path.join(asa_main_output_dir, "s2_asa_mastervariant_gt/")

script_3_input_dir = script_2_output_dir
script_3_output_dir = os.path.join(asa_main_output_dir, "s3_asa_mastervariant_grr/")

script_4_input_dir = script_3_output_dir
script_4_output_dir = os.path.join(asa_main_output_dir, "s4_asa_mastervariant_TrtscrMkrGaugeMeas/")

script_5_input_dir = script_4_output_dir
script_5_output_dir = os.path.join(asa_main_output_dir, "s5_asa_mastervariant_finalresults/")

## Create sub output directories if they don't exist
if not os.path.exists(script_1a_output_dir):
    os.mkdir(script_1a_output_dir)
if not os.path.exists(script_1b_output_dir):
    os.mkdir(script_1b_output_dir)
if not os.path.exists(script_2_output_dir):
    os.mkdir(script_2_output_dir)
if not os.path.exists(script_3_output_dir):
    os.mkdir(script_3_output_dir)
if not os.path.exists(script_4_output_dir):
    os.mkdir(script_4_output_dir)
if not os.path.exists(script_5_output_dir):
    os.mkdir(script_5_output_dir)

logger.info('Starting to run Enterprise System Scripts')

## Run subscripts
logger.info('Split Asa raw file by sample name - Running')
splitAsaBySampleName(asa_raw_file, popfile, script_1a_output_dir)
logger.info('Split Asa raw file by sample name - Success!')

logger.info('Merge Allele Plus 1 & 2 to mastervariant by Snpname - Running')
extractAsaAlleleplusBySnpname(script_1b_input_dir, script_1b_output_dir)
logger.info('Merge Allele Plus 1 & 2 to mastervariant by Snpname - Success!')

logger.info('Generating VCF format - Running')
generateVcfFormat(script_2_input_dir, script_2_output_dir)
logger.info('Generating VCF format - Success!')

logger.info('Calculating genetic relative risk - Running')
generateGrr(popfile, script_3_input_dir, script_3_output_dir)
logger.info('Calculating genetic relative risk - Success!')

logger.info('Calculating Trait Score, Marker, Total Marker, and assigning Gauge and Gauge Measurement - Running')
generateTrtscrMarkerGauge(args.traitid_file, script_4_input_dir, script_4_output_dir)
logger.info('Calculating Trait Score, Marker, Total Marker, and assigning Gauge and Gauge Measurement - Success!')

logger.info('Selecting columns and filtering rows for final results - Running')
extractFinalResults(script_5_input_dir, script_5_output_dir)
logger.info('Selecting columns and filtering rows for final results - Success!')

### Can uncomment the following codes to remove intermediate directories when necessary
## Function for Exception Handling
# def handler(func, path, exc_info):
#     print('We got the following exception')
#     print(exc_info)

## Removing directories
# logger.info('Removing intermediate files - Running')
# shutil.rmtree(script_1a_output_dir, ignore_errors=False, onerror=handler)
# shutil.rmtree(script_1b_output_dir, ignore_errors=False, onerror=handler)
# shutil.rmtree(script_2_output_dir, ignore_errors=False, onerror=handler)
# shutil.rmtree(script_3_output_dir, ignore_errors=False, onerror=handler)
# logger.info('Removing intermediate files - Success!')

logger.info('Process completed!')