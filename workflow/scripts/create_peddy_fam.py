#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys

def get_sex(sex):

    if sex == 'M':
        sex_code = '1'
    elif sex == 'K':
        sex_code = '2'
    else:
        sex_code = '0'

    return sex_code

df = pd.read_csv(snakemake.input[0], header=13,
  dtype=str).set_index("Sample_ID", drop=False)

df[['Code', 'Sex', 'NA', '123', 'Trio']] = df.Description.str.split('_', expand=True)

trio_ids = df.Trio.tolist() # checking for trio ids in the Trio column
if len(set(trio_ids)) == 1 and trio_ids[0] == 'NA' :
    print('No Trios detected in Sample Sheet')
    df['Trio_Status'] = ['NA']*df.shape[0]
    df['FID'] = ['NA']*df.shape[0]
else:
    df[['FID', 'Trio_Status']] = df.Trio.str.split('-', expand=True)

fam_df  = df[['Sample_ID', 'Sex', 'FID', 'Trio_Status']]

child_df = fam_df[df.Trio_Status== 'Index']
father_df = fam_df[(df.Trio_Status == 'Foralder') & (df.Sex == 'M')]
mother_df = fam_df[(df.Trio_Status == 'Foralder') & (df.Sex == 'K')]

with open(snakemake.output[0], 'w') as pedfile:
    phenotype = '0'

    for sample in fam_df.itertuples():
        sample_id = sample.Sample_ID + '_N' # add the '_N'
        family_id = sample.FID
        if  family_id == 'NA':
             family_id = sample_id

        if sample.Sample_ID in list(child_df.Sample_ID):
            child_trio = sample.FID
            paternal_id = father_df[father_df.FID == child_trio].iat[0,0] + '_N'
            maternal_id = mother_df[mother_df.FID == child_trio].iat[0,0] + '_N'
        else:
            paternal_id = '0'
            maternal_id = '0'

        sex = get_sex(sample.Sex)
        print(family_id, sample_id ,paternal_id, maternal_id,
         sex, phenotype, sep='\t', file=pedfile )
