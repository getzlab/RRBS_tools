from bs4 import BeautifulSoup
import sys

br_file = sys.argv[1]
#br_file = 'toy500K_bismark_report.html'

br = BeautifulSoup(open(br_file),'html.parser')
sequence_table = br.body.table.find_all('td')
string_table = br.body.table.find_all('th')
for i in range(0,len(sequence_table)) :
    print('%s : %s' % (string_table[i].string, sequence_table[i].string))
    
total_reads = sequence_table[0].string
unique_reads = sequence_table[1].string
unaligned_reads = sequence_table[2].string
non_unique_reads = sequence_table[3].string
non_extractable_reads = sequence_table[4].string
