#!/usr/bin/env python3
import pandas as pd
import sys
from bs4 import BeautifulSoup

br_file = sys.argv[1]
paired_string = sys.argv[2]


br = BeautifulSoup(open(br_file),'html.parser')
sequence_table = br.body.table.find_all('td')
string_table = br.body.table.find_all('th')
for i in range(0,len(sequence_table)) :
    print('%s : %s' % (string_table[i].string, sequence_table[i].string))
    
total_reads = int(sequence_table[0].string)
unique_reads = int(sequence_table[1].string)
unaligned_reads = int(sequence_table[2].string)
non_unique_reads = int(sequence_table[3].string)
non_extractable_reads = int(sequence_table[4].string)



sample_id = str.replace(br_file,'_bismark_report.html','')
report = pd.DataFrame({'sample_id' : sample_id,
'paired_str' : paired_string,
'total_reads' : total_reads,
'aligned_reads': unique_reads + non_unique_reads,
'aligned_reads_pc' : (unique_reads + non_unique_reads)/total_reads,
'mapped_unique' : unique_reads,
'mapped_unique_pc' : unique_reads/total_reads,
'mapped_non_unique' : non_unique_reads,
'mapped_non_unique_pc' : non_unique_reads/total_reads },index=[0])

report.to_csv('%s_bismark_bamstats.tsv' % sample_id,sep='\t',index=False)
