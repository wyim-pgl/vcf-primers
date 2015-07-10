from subprocess import Popen, PIPE
from utils import *
from pprint import pprint as pp
# Wrapper for primer3


class primer3:

    def __init__(self, method="sanger"):
        # Sanger Default
        if method == "sanger":
            self.PRIMER_PRODUCT_SIZE_RANGE = "600-800"  # Must be set
            self.PRIMER_OPT_SIZE = 20
            self.PRIMER_MIN_SIZE = 18  # Must be set
            self.PRIMER_MAX_SIZE = 20  # Must be set

        # Global default
        seq_template_length = self.PRIMER_PRODUCT_SIZE_RANGE.split("-")[1]
        self.seq_template_length = int(seq_template_length)
        thermo_path = "/usr/local/share/primer3_config/"
        self.PRIMER_THERMODYNAMIC_PARAMETERS_PATH = thermo_path

    def generate_record(self):
        # Generates text record ready for input
        # into primer3
        attributes = [x for x in dir(self) if x.upper() == x
                      and not x.startswith("_")]
        values = [str(getattr(self, x)) for x in attributes]
        att_val = zip(attributes, values)
        return '\n'.join(["=".join(x) for x in att_val]) + "\n=\n"

    def fetch_primers(self):
        # Runs primer3 with the generated record.
        primer3_run = Popen(["primer3_core"], stdin=PIPE, stdout=PIPE)
        record = self.generate_record()
        resp, err = primer3_run.communicate(record)
        resp = resp.strip().split("\n")
        if err:
            raise Exception(err)
        primer3_results = dict([x.split("=") for x in resp
                               if x.split("=")[0] != ""])
        self.results = {}
        for k, v in primer3_results.items():
            self.results[k] = autoconvert(v)

    def primers(self):
        # Iterates through primers
        if "PRIMER_LEFT_NUM_RETURNED" in self.results:
            n_primers = self.results["PRIMER_LEFT_NUM_RETURNED"]
            print pp(self.results)
            for primer_num in xrange(0,n_primers):
                primer_return = {}
                left_right = ["PRIMER_LEFT","PRIMER_RIGHT"]
                for side in left_right:
                    primer_num = str(primer_num)
                    primer_return[side] = {}
                    for key, val in [(k.replace(side + "_" + primer_num + "_",""), v) for k,v in self.results.items() 
                        if k.startswith(side + "_" + str(primer_num))]:
                            primer_return[side][key] = val
                yield primer_return
                

"""

               
                                    # Blast primer

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
                        
                                    print pp(r)"""

