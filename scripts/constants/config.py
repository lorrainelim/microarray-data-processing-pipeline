### These column names do not represent all columns of the tables involved. The columns listed below are those that are used in the scripts. 

## Declare mastervariant path
mastervariant_file = 'database/xx_mastervariant.txt'

## Column names of ASA raw file (that will be extracted/used in this pipeline)

asa_raw_sample_name = 'Sample Name'
asa_raw_snpname = 'SNP Name'
asa_raw_allele1plus = 'Allele1 - Plus'
asa_raw_allele2plus = 'Allele2 - Plus'

## Column names of population file

pop_sample_name = 'sample_name'
pop_population = 'population'

## Column names of mastervariant & results

snpname = 'SNPNAME'
ref = 'REF'
alt = 'ALT'
risk_code = 'RISK_CODE'
allele_type = 'ALLELE_TYPE'
gene ='GENE'
odd_ratio = 'OR'

asa_sample_name = 'SAMPLE_NAME'
allele_1_plus = 'ALLELE_1_PLUS'
allele_2_plus = 'ALLELE_2_PLUS'
asa_genotype = 'ASA_GENOTYPE'

allele_1_gt = 'Allele1_GT'
allele_2_gt = 'Allele2_GT'
asa_gt = 'ASA_GT'

grr = 'GRR'

## Column names of Trait ID table

trait_id = 'TRAIT_ID'
section = 'SECTION'
section_n_trait = 'SECTION_N_TRAIT'
assoc_type = 'ASSOC_TYPE'
ge1 = 'GE1'	
ge2 = 'GE2'	
ge3 = 'GE3'	
ge4 = 'GE4'	

## Column names for calculation

trait_score = 'TRAIT_SCORE'	
marker_num = 'MARKER_NUM'	
marker_total = 'MARKER_TOTAL'	
gauge = 'GAUGE'	
genotypic_effect = 'GENOTYPIC_EFFECT'

## Columns that are used just for rearrangements

no_chr = 'No_CHR'	
chrom = 'CHROM'	
pos37 = 'POS37'	
pos38 = 'POS38'	
rsid = 'RSID'
risk_allele = 'RISK_ALLELE'
afr = 'AFR'
amr = 'AMR'
eas = 'EAS'
eur = 'EUR'
sas = 'SAS'
