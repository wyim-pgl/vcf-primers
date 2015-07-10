#! /usr/bin/python
"""vcf-primer

Usage:
  vcf-primer.py <vcf> <reference>
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


debug = True
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

    # If running program without pipe, use bcftools to open.
    vcf_comm = ["bcftools","query","--regions", "I:1000-2000", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", args["<vcf>"]]
    vcf = Popen(vcf_comm, stdout=PIPE, stderr=PIPE)
    blaster = blastn(reference)

    # Calculate variables
    for line in vcf.stdout:
        # SEQUENCE_ID; 
        chrom, pos, ref, alt = line.strip().split("\t")
        p.SEQUENCE_ID = chrom + "_" + pos
        
        # SEQUENCE_TEMPLATE 
        start = int(pos) - p.seq_template_length
        end = int(pos) + len(ref) + p.seq_template_length
        p.SEQUENCE_TEMPLATE = faidx(reference, chrom, start, end)
        p.SEQUENCE_TARGET = "%s,%s" % (p.seq_template_length, len(ref))
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
              sequence_output = "{pre_variant}[{ref}/{alt}]{post_variant}...".format(**locals())

              print sequence_output
