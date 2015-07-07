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

# Check that tools are installed
for tool in ["primer3_core", "blastn", "samtools"]:
  if Popen(["which",tool],stdout=PIPE, stderr=PIPE).communicate()[0] == "":
      raise Exception(tool + " is not installed.")

# Set up rec and defaults
rec = {"PRIMER_OPT_SIZE":18}

seq_context_length = 100

if __name__ == '__main__':
    args = docopt(__doc__, version='vcf-primer 0.1')
    #print(args)

    # If running program without pipe, use bcftools to open.
    vcf = Popen(["bcftools","query","-f", "%CHROM\t%POS\t%REF\t%ALT\n", args["<vcf>"]], stdout=PIPE, stderr=PIPE)
    for line in vcf.stdout:

        # SEQUENCE_ID; Later make rsid?
        chrom, pos, ref, alt = line.strip().split("\t")
        rec["SEQUENCE_ID"] = chrom + "_" + pos

        # SEQUENCE_TEMPLATE (200 bp region)
        pos_start = str(int(pos) - seq_context_length)
        pos_end = str(int(pos) + seq_context_length)
        loc = "%s:%s-%s" % (chrom, pos_start, pos_end)
        seq = ''.join(check_output(["samtools","faidx", args["<reference>"], loc]).splitlines()[1:])
        rec["SEQUENCE_TEMPLATE"] = seq
        rec["SEQUENCE_TARGET"] = "%s,%s" % (seq_context_length, len(ref))
        rec["PRIMER_TASK"] = "pick_pcr_primers"
        rec["PRIMER_MIN_SIZE"] = 18
        rec["PRIMER_MAX_SIZE"] = 22
        rec["PRIMER_THERMODYNAMIC_PARAMETERS_PATH"] = "/usr/local/share/primer3_config/" # From homebrew
        record = '\n'.join(["%s=%s" % (k,v) for k,v in rec.items()]) + "\n="
        print record




else:
  for line in sys.stdin:
      if line.startswith("#") == False:
          sys.stdout.write(line)
