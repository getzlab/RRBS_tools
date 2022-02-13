#!/usr/bin/env python3
import re
import pandas as pd
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--sample_id", dest="sample_id",
                  help="sample_id, added to output")
parser.add_option("-f", "--bamstats_file", dest="report_file",
                  help="bamStats.pl output file to parse")
parser.add_option("-o", "--outfile", dest="outfile",
                  help="output file (tsv format)")
parser.add_option("-t", "--seq_type", dest="seq_type", default="PE",
                help="PE or SE")

(opts, args) = parser.parse_args()


def parse_bamstats(sample_id, report_file, out_tsv=None, regs=None, seq_type='PE', return_dict=False):

    resdict = {'sample_id':[sample_id], 'paired_str': seq_type}
    if regs==None:
        regs = { 'total_reads': '\s*Number of Sequences analyzed \.\.:\s+(\d+,?\d*,?\d*,?\d*,?\d*)',
                    'aligned_reads': '\s*Mapped Reads \.*:\s+(\d+,?\d*,?\d*,?\d*,?\d*)',
                    'aligned_reads_pc': '\s*Mapped Reads \.*:\s+\d+,?\d*,?\d*,?\d*,?\d* \((\d+\.?\d*)%\)',
                    'markdup_reads': '\s*Marked \"Duplicated\" \.*:\s+(\d+,?\d*,?\d*,?\d*,?\d*)',
                    'markdup_reads_pc': '\s*Marked \"Duplicated\" \.*:\s+(\d+,?\d*,?\d*,?\d*,?\d*) \((\d+\.?\d*)%\)',
                    'mapped_read_pairs': '\s*Mapped Read Pairs \.*:\s+(\d+,?\d*,?\d*,?\d*,?\d*)',
                    'mapped_read_pairs_pc': '\s*Mapped Read Pairs \.*:\s+\d+,?\d*,?\d*,?\d*,?\d* \((\d+\.?\d*)%\)',
                    'unmapped_reads': '\s*Unmapped Reads \.*:\s+(\d+,?\d*,?\d*,?\d*,?\d*)',
                    'unmapped_reads_pc': '\s*Unmapped Reads \.*:\s+\d+,?\d*,?\d*,?\d*,?\d* \((\d+\.?\d*)%\)',
                    'average_read_len': '\s*Average Read Length \.*:\s+(\d+)'}

    f = open(report_file, "r")
    lines = f.readlines()
    f.close()

    for l in lines:
        for k in regs.keys():
            if re.search(regs[k], l):
                resdict[k] = [re.search(regs[k], l).group(1)]

    if(out_tsv != None):
        pd.DataFrame.from_dict(resdict).to_csv(out_tsv, sep='\t', index=False)

    if return_dict:
        return(resdict)

#Parse and write new output file
parse_bamstats(opts.sample_id, opts.report_file, out_tsv=opts.outfile, seq_type=opts.seq_type)
