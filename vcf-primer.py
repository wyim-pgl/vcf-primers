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
from pprint import pprint
from subprocess import Popen, PIPE, check_output
from pprint import pprint as pp
from primer3_options import *
from blastn import *
from utils import *

# Check that tools are installed
for tool in ["primer3_core", "blastn", "samtools"]:
  if Popen(["which",tool],stdout=PIPE, stderr=PIPE).communicate()[0] == "":
      raise Exception(tool + " is not installed.")

# Set up rec and defaults
rec = {"PRIMER_OPT_SIZE":20}

seq_template_length = int(PRIMER_PRODUCT_SIZE_RANGE.split("-")[1]) # Region that primer3_core searches for primers.

if __name__ == '__main__':
    args = docopt(__doc__, version='vcf-primer 0.1')
    #print(args)

    # If running program without pipe, use bcftools to open.
    vcf = Popen(["bcftools","query","-f", "%CHROM\t%POS\t%REF\t%ALT\n", args["<vcf>"]], stdout=PIPE, stderr=PIPE)

    rec["PRIMER_PRODUCT_SIZE_RANGE"] = PRIMER_PRODUCT_SIZE_RANGE
    rec["PRIMER_TASK"] = "pick_pcr_primers"
    rec["PRIMER_MIN_SIZE"] = PRIMER_MIN_SIZE
    rec["PRIMER_MAX_SIZE"] = PRIMER_MAX_SIZE
    rec["PRIMER_THERMODYNAMIC_PARAMETERS_PATH"] = "/usr/local/share/primer3_config/" # From homebrew

    blaster = blastn(args["<reference>"].replace(".gz",""))

    # Calculate variables

    for line in vcf.stdout:

        # SEQUENCE_ID; Later make rsid?
        chrom, pos, ref, alt = line.strip().split("\t")
        rec["SEQUENCE_ID"] = chrom + "_" + pos

        # SEQUENCE_TEMPLATE 
        pos_start = str(int(pos) - seq_template_length)
        pos_end = str(int(pos) + len(ref) + seq_template_length)
        loc = "%s:%s-%s" % (chrom, pos_start, pos_end)
        seq = ''.join(check_output(["samtools","faidx", args["<reference>"], loc]).splitlines()[1:])
        rec["SEQUENCE_TEMPLATE"] = seq
        rec["SEQUENCE_TARGET"] = "%s,%s" % (seq_template_length, len(ref))
        record = '\n'.join(["%s=%s" % (k,v) for k,v in rec.items()]) + "\n="
        p = Popen(["primer3_core"],stdin=PIPE, stdout=PIPE)
        resp, err = p.communicate(record)
        r = dict([x.split("=") for x in resp.strip().split("\n") if x.split("=")[0] != ""])
        for k,v in r.items():
          r[k] = autoconvert(v)

        for primer_num in xrange(0,r["PRIMER_LEFT_NUM_RETURNED"]):
          # Blast primer
          left_primer = r["PRIMER_LEFT_%s_SEQUENCE" % primer_num]
          right_primer = r["PRIMER_RIGHT_%s_SEQUENCE" % primer_num]
          if left_primer is None or right_primer is None:
            break
          left_primer_blast = blaster.query(left_primer)
          right_primer_blast = blaster.query(right_primer)
          print left_primer_blast["evalue"]
          # Generate PCR Templates
          left_start, left_length = map(int,r["PRIMER_LEFT_%s" % primer_num].split(","))
          right_start, right_length = map(int,r["PRIMER_RIGHT_%s" % primer_num].split(","))
          r["PCR_PRODUCT_SEQ_%s" % primer_num] = r["SEQUENCE_TEMPLATE"][left_start:right_start+1]
          r["PCR_PRODUCT_LEN_%s" % primer_num] = len(r["SEQUENCE_TEMPLATE"][left_start:right_start+1])


else:
  for line in sys.stdin:
      if line.startswith("#") == False:
          sys.stdout.write(line)
