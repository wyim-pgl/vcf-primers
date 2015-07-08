from subprocess import Popen, PIPE
from utils import *
import csv


field = ["query_id",
         "CHROM", # subject_id
         "perc_identity",
         "alignment length",
         "mismatches",
         "gap opens",
         "q_start",
         "q_end",
         "s_start",
         "s_end",
         "evalue",
         "bit_score"]



class blastn:

    def __init__(self, db):
        self.db = db

    def query(self, q):
        blastn_query = "echo {q} | blastn -query - -db={self.db} -outfmt=6 -evalue 1 -word_size=7".format(**locals())
        resp, err = Popen(blastn_query, stdout=PIPE, stderr=PIPE, shell=True).communicate()

        if err != "":
            raise Exception(err)


        resp = [dict(zip(field,map(autoconvert,x.split("\t")))) for x in resp.strip().split("\n")]
        for i in resp:
            i["query_length"] = i["q_end"] + 1 - i["q_start"]
            i["query_string"] = q[i["q_start"]-1:i["q_end"]]
        return resp

#blaster = blastn("reference/WS245/WS245.fa")
#blaster.query("gcctgcctgccaacctatat")
