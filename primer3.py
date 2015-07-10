from subprocess import Popen, PIPE
from utils import *
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
        return '\n'.join(["=".join(x) for x in att_val])

    def fetch_primers(self):
        # Runs primer3 with the generated record.
        primer3_run = Popen(["primer3_core"], stdin=PIPE, stdout=PIPE)
        record = self.generate_record()
        resp, err = primer3_run.communicate(record)
        resp = resp.strip().split("\n")
        if err:
            raise Exception(err)
        r = dict([x.split("=") for x in resp if x.split("=")[0] != ""])
        for k, v in r.items():
            r[k] = autoconvert(v)
