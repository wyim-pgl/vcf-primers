#! /usr/bin/python
"""vcf-primer

Usage:
  vcf-primer.py <reference> [--header]
  vcf-primer.py (-h | --help)
  vcf-primer.py --version

Options:
  -h --help     Show this screen.
  --version     Show version.

"""
from docopt import docopt
import sys
from subprocess import Popen, PIPE, check_output
from pprint import pprint as pp
from blastn import *
from utils import *
from primer3 import *


debug = False
# Check that tools are installed
for tool in ["primer3_core", "blastn", "samtools"]:
  if Popen(["which",tool],stdout=PIPE, stderr=PIPE).communicate()[0] == "":
      raise Exception(tool + " is not installed.")

p = primer3("sanger")

if __name__ == '__main__':
    if len(sys.argv) == 1 and debug == True:
      args = {"<reference>": "reference/WS245/WS245.fa", 
              "<vcf>": "test.vcf.gz"}
    else:
      args = docopt(__doc__, version='vcf-primer 0.1')
    reference = args["<reference>"]

    if sys.stdin:
      vcf = sys.stdin
    else:
      # If running program without pipe, use bcftools to open.
      vcf_comm = ["bcftools","query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", args["<vcf>"]]
      vcf = Popen(vcf_comm, stdout=PIPE, stderr=PIPE).stdout
    blaster = blastn(reference)
    # Calculate variables
    for line in vcf:
        print line
        # SEQUENCE_ID; 
        if line.startswith("##"):
          continue
        if line.startswith("#CHROM"):
            samples = line.strip().split("\t")[9:]
            continue
        line = line.strip().split("\t")
        chrom, pos, ref, alt, gt = line[0], line[1], line[3], line[4], [x.split(":")[0] for x in line[9:]]
        gt = zip(samples, gt)
        p.SEQUENCE_ID = chrom + "_" + pos
        
        # SEQUENCE_TEMPLATE 
        start = int(pos) - p.seq_template_length
        end = int(pos) + len(ref) + p.seq_template_length
        p.SEQUENCE_TEMPLATE = faidx(reference, chrom, start, end)
        p.SEQUENCE_TARGET = "%s,%s" % (p.seq_template_length, len(ref))
        if args['--header']:
          columns = ["CHROM",
                     "POS",
                     "REF",
                     "ALT",
                     "LEFT_PRIMER",
                     "LEFT_PERFECT_MATCHES",
                     "LEFT_BP_MISMATCH",
                     "LEFT_END_MISMATCH",
                     "LEFT_NEARBY_MATCH",
                     "RIGHT_PRIMER",
                     "RIGHT_PERFECT_MATCHES",
                     "RIGHT_BP_MISMATCH",
                     "RIGHT_END_MISMATCH",
                     "RIGHT_NEARBY_MATCH",
                     "SANGER_PRIMER",
                     "SEQ_RESULT",
                     "REF_STRAINS",
                     "ALT_STRAINS",
                     "PCR_PRODUCT_LEN",
                     "PCR_PRODUCT_SEQ",
                     "REF_STRAINS",
                     "HET_STRAINS",
                     "ALT_STRAINS",
                     "MISSING_STRAINS"
                     ]
          print "\t".join(columns)
        # Fetch primers
        for primer in p.fetch_primers():
            left_primer  = blaster.blast(primer["PRIMER_LEFT"]["SEQUENCE"]).query_stat()
            right_primer = blaster.blast(primer["PRIMER_RIGHT"]["SEQUENCE"]).query_stat()
            
            interior_primer = primer3("sanger_interior")
            interior_primer.SEQUENCE_TEMPLATE = primer["PCR_PRODUCT_SEQ"]
            for int_primer in interior_primer.fetch_primers():
              int_primer = int_primer["PRIMER_LEFT"]
              int_primer_seq = int_primer["SEQUENCE"] 
              int_primer_count = p.SEQUENCE_TEMPLATE.count(int_primer_seq)
              # If the interior primer is found multiple times - don't use it
              if int_primer_count > 1:
                  break
              #print p.SEQUENCE_ID
              int_primer_len = int_primer["PRIMER_LENGTH"]
              seq_start = p.SEQUENCE_TEMPLATE.find(int_primer_seq) + len(int_primer_seq)
              target_location = p.seq_template_length - seq_start
              pre_variant = p.SEQUENCE_TEMPLATE[seq_start-1:seq_start + target_location]
              post_variant = p.SEQUENCE_TEMPLATE[seq_start + target_location+1:seq_start + target_location + 10]
              REF_STRAINS = '|'.join([x[0] for x in gt if x[1] == "0/0"])
              HET_STRAINS = '|'.join([x[0] for x in gt if x[1] == "0/1"])
              ALT_STRAINS = '|'.join([x[0] for x in gt if x[1] == "1/1"])
              MISSING_STRAINS = "|".join([x[0] for x in gt if x[1] == "./."])
              sequence_output = [chrom,
                                 pos,
                                 ref, 
                                 alt,
                                 primer["PRIMER_LEFT"]["SEQUENCE"],
                                 left_primer["perfect_match"],
                                 left_primer["bp_mismatch"],
                                 left_primer["end_mismatch"],
                                 left_primer["nearby_matches"],
                                 primer["PRIMER_RIGHT"]["SEQUENCE"],
                                 right_primer["perfect_match"],
                                 right_primer["bp_mismatch"],
                                 right_primer["end_mismatch"],
                                 right_primer["nearby_matches"],
                                 int_primer["SEQUENCE"],
                                 "{pre_variant}[{ref}/{alt}]{post_variant}".format(**locals()),
                                 primer["PCR_PRODUCT_LEN"],
                                 primer["PCR_PRODUCT_SEQ"],
                                 REF_STRAINS,
                                 HET_STRAINS,
                                 ALT_STRAINS
                                 ]
              #sequence_output = "{chrom}\t{pos}\t{ref}\t{alt}\t{pre_variant}[{ref}/{alt}]{post_variant}\t".format(**locals())
              print '\t'.join([str(x) for x in sequence_output])
