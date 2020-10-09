#!/usr/bin/env python3
import re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--sample_id", dest="sample_id",
                  help="sample_id, added to output")
parser.add_option("-f", "--report_file", dest="report_file",
                  help="bsmap report to parse")
parser.add_option("-o", "--outfile", dest="outfile",
                  help="output file (tsv format)")

(opts, args) = parser.parse_args()


def parse_bsmap_report(sample_id, report_file, out_tsv=None, regs=None, return_dict=False):

    resdict = {'sample_id':[sample_id]}
    if regs==None:
        regs = {'paired_str':'\[bsmap\].*(Single-end|Paired-end) alignment',
                    'total_reads':'\[bsmap\].*total reads: (\d+)',
                    'aligned_reads':'aligned reads: (\d+\.?\d*) \(\d+\.?\d*%\)',
                    'aligned_reads_pc':'aligned reads: \d+\.?\d* \((\d+\.?\d*)%\)',
                    'mapped_unique':'unique reads: (\d+\.?\d*) \(\d+\.?\d*%\)',
                    'mapped_unique_pc':'unique reads: \d+\.?\d* \((\d+\.?\d*)%\)',
                    'mapped_non_unique':'non-unique reads: (\d+\.?\d*) \(\d+\.?\d*%\)',
                    'mapped_non_unique_pc':'non-unique reads: \d+\.?\d* \((\d+\.?\d*)%\)'}

    f = open(report_file, "r")
    lines = f.readlines()
    f.close()

    for l in lines:
        for k in regs.keys():
            if re.search(regs[k], l):
                resdict[k] = [re.search(regs[k], l).group(1)]

    if(out_tsv != None):
        pd.DataFrame.from_dict(resdict).to_csv(out_tsv, sep='\t')

    if return_dict:
        return(resdict)

#Parse and write new output file
parse_bsmap_report(opts.sample_id, opts.report_file, out_tsv=opts.outfile)
