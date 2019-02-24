#!/usr/bin/env python2
import os
import time
from multiprocessing import Pool

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq

os.chdir("/home/xue/Documents/pe150/amplicon_ngs")

def do_something_with_record(info):
    name, seq = info
    time.sleep(10)
    return name, seq

def convert_to_fasta(in_handle):
    for rec_id, seq, _ in FastqGeneralIterator(in_handle):
        yield rec_id, str(Seq(seq).reverse_complement())

with open("Cleandata/chip/chip_R1.fq") as input_handle:
    p = Pool(processes=20)
    g = p.map(do_something_with_record, convert_to_fasta(input_handle))

# for i in g:
#     print i