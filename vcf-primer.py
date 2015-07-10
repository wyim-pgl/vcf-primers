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

# Check that tools are installed
for tool in ["primer3_core", "blastn", "samtools"]:
  if Popen(["which",tool],stdout=PIPE, stderr=PIPE).communicate()[0] == "":
      raise Exception(tool + " is not installed.")

p = primer3("sanger")

if __name__ == '__main__':
    args = docopt(__doc__, version='vcf-primer 0.1')
    
    reference = args["<reference>"].replace(".gz","")
    reference_gz = args["<reference>"]

    # If running program without pipe, use bcftools to open.
    vcf_comm = ["bcftools","query","-f", "%CHROM\t%POS\t%REF\t%ALT\n", args["<vcf>"]]
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
        p.SEQUENCE_TEMPLATE = faidx(args["<reference>"], chrom, start, end)
        p.SEQUENCE_TARGET = "%s,%s" % (p.seq_template_length, len(ref))

        

        if "PRIMER_LEFT_NUM_RETURNED" in r:
          for primer_num in xrange(0,r["PRIMER_LEFT_NUM_RETURNED"]):
            # Blast primer
            left_primer = r["PRIMER_LEFT_%s_SEQUENCE" % primer_num]
            right_primer = r["PRIMER_RIGHT_%s_SEQUENCE" % primer_num]
            if left_primer is None or right_primer is None:
              break
            left_primer_blast  = blaster.blast(left_primer).query_stat()
            right_primer_blast = blaster.blast(right_primer).query_stat()

            # Generate PCR Templates
            left_start, left_length = map(int,r["PRIMER_LEFT_%s" % primer_num].split(","))
            right_start, right_length = map(int,r["PRIMER_RIGHT_%s" % primer_num].split(","))
            r["PCR_PRODUCT_SEQ_%s" % primer_num] = r["SEQUENCE_TEMPLATE"][left_start:right_start+1]
            r["PCR_PRODUCT_LEN_%s" % primer_num] = len(r["SEQUENCE_TEMPLATE"][left_start:right_start+1])

            # Generate primer for interior
            #rec["SEQUENCE_TEMPLATE"] = r["PCR_PRODUCT_SEQ_%s" % primer_num]
            #rec["SEQUENCE_TARGET"] = 1

            print pp(r)

else:
  for line in sys.stdin:
      if line.startswith("#") == False:
          sys.stdout.write(line)
