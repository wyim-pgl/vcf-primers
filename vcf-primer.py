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

if __name__ == '__main__':
    args = docopt(__doc__, version='vcf-primer 0.1')
    print(args)

    # If running program without pipe, use bcftools to open.
    vcf = Popen(["bcftools","query","-f", "%CHROM\t%POS\n", args["<vcf>"]], stdout=PIPE, stderr=PIPE)
    for line in vcf.stdout:

    	# SEQUENCE_ID; Later make rsid?
    	chrom, pos = line.strip().split("\t")
    	rec["SEQUENCE_ID"] = chrom + "_" + pos

    	# SEQUENCE_TEMPLATE (200 bp region)
    	pos_start = str(int(pos) - 100)
    	pos_end = str(int(pos) + 100)
    	loc = "%s:%s-%s" % (chrom, pos_start, pos_end)
    	seq = check_output(["samtools","faidx", args["<reference>"], loc]).splitlines()
    	rec["SEQUENCE_TEMPLATE"] = seq

    	# SEQUENCE_TARGET

    	# PRIMER TASK
    	rec["PRIMER_TASK"] = "generic"



else:
	for line in sys.stdin:
	    if line.startswith("#") == False:
	        sys.stdout.write(line)
