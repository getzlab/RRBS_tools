#!python3
import re
def parse_bsmap_report(sample_id, report_file, out_tsv):

    resdict = {'sample_id':sample_id}
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
                resdict[k] = re.search(regs[k], l).group(1)

    #out_fh = open(out_tsv, "w")
    #out_fh.close()
    return(resdict)
