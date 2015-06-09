#! /usr/bin/python

import sys
from subprocess import Popen, PIPE

# Check that primer3 is installed.
if Popen(["which","primer3_core"],stdout=PIPE, stderr=PIPE).communicate()[0] == "":
    raise Exception("Primer3 is not installed.")

for line in sys.stdin:
    if line.startswith("#") == False:
        sys.stdout.write(line)