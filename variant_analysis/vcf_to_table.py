####################################
# Converting VCF files to CSV Table
# 
# Author: Jessica Luc, 2022
####################################

import pandas as pd  
import sys


if len(sys.argv) < 1 :
    print( 'vcf_to_table.py <.csv file> #first remove lines until #CHROM table')
    exit()

df  = pd.read_csv(sys.argv[1], sep=",")
new_column_names = ['POS','REF','ALT','QUAL']

sample_nums = list(range(1,48)) #change depending on number of samples
#for num in sample_nums:
#    modified = 'sample ' + str(num)
#    new_column_names.append(modified)
new_column_names.extend(sample_nums)
out = pd.DataFrame(columns = new_column_names)


out['POS'] = df['POS'].values
out['REF'] = df['REF'].values
out['ALT'] = df['ALT'].values
out['QUAL'] = df['QUAL'].values

for row in range(0, len(df)):
    for sample in sample_nums:
        #col = df.loc[row, str(sample)+'.sorted.bam']
        col = str(df.loc[row, str(sample)])
        if '.:.' in col:
          out.loc[row, sample] = 0
        else:
          metadata = col.split(':')
          #use if statement for AO or proportion only
          if ',' in df.loc[row, 'ALT']:
              freq = freq.split(',')
              #for AO:
              #out.loc[row, 'sample ' + sample] = sum(list(map(int,freq)))
              #for proportion:
              out.loc[row, sample] = sum(list(map(int,freq)))/int(metadata[1])
          else:
              #for AO:
              #out.loc[row, 'sample ' + sample] = metadata[5]
              #for proportion:
              out.loc[row, sample] = int(metadata[5])/int(metadata[1])
          #for DP:
          #out.loc[row, 'sample ' + str(sample)] = metadata[1]
outputname = 'filename.csv'
out.to_csv(outputname, sep = ',', index = False)
print("Saved file:", outputname)
