#!/usr/bin/env python
import sys
import re
from Bio.Seq import Seq
from itertools import izip
from highlight_keywords_table import *


with open('../chip_R1.fq') as f1, open('../chip_R2.fq') as f2:
    lineNumber = 0
    for r1, r2 in izip(f1, f2):
        lineNumber += 1
        if (lineNumber - 2) % 4 == 0:
            print("R1 line:{2}\n{0}\nR2 line:{2}\n{1}".format(r1.strip(), r2.strip(), lineNumber))
            print '---'
            if raw_input() == 'q':
                break
