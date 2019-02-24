#!/usr/bin/env python2
from itertools import izip, islice
import regex
import os
import difflib

os.chdir("/home/xue/Documents/pe150/amplicon_ngs")

# spacer added to VHH barcodes
barcodesVHH = ["TAAGGCGA", "CGTACTAGT", "AGGCAGAATA", "TCCTGAGCATG", "GGACTCCTATGC"]
overlappedVHHHead = ["ACCGTGGCCCAGGCGGCCCAG"]
overlappedVHHEnd = ["ATGGTGCTGGCCGGCCTGGCC"]
barcodesHingeHead = ["TAGGCATG"]
# spacer added
barcodesHingeEnd = ["GCGATCTATA", "ATAGAGAGT", "TCTACTCT"]

class Reads:

    # return a dict only records one read info, therefore no use of class
    # def snippetPacker(lines=None):
    #     ks = ['name', 'sequence', 'optional', 'quality']
    #     return {k: v for k, v in zip(ks, lines)}

    def __init__(self, lines, lineNumber):
        self.name, inputSequence, self.optional, self.quality = lines
        self.sequence = inputSequence.strip()
        self.lineNumber = lineNumber


class Seeds(str):

    def __init__(self, seed, errors = 1):
        self.seed = seed
        self.errors = errors

    # inherit from string and perform like a string
    def __new__(cls, value):
        obj = str.__new__(cls, value)
        return obj

    @property
    def pattern(self):
        # double {{ to escape, the ENHANCEMATCH flag (?e) attemp to better match
        return ".*(?e)({}){{e<={}}}".format(self.seed, self.errors)


def readR1R2():
    with open('Cleandata/chip/chip_R1.fq') as f1, open('Cleandata/chip/chip_R2.fq') as f2:
        # lines1 = islice(f1, 1, None, 4)
        # lines2 = islice(f2, 1, None, 4)
        snippetLines1, snippetLines2 = ([], [])
        linesCount = 0
        for lineInsnippet1, lineInsnippet2 in izip(f1, f2):
            linesCount += 1
            snippetLines1.append(lineInsnippet1)
            snippetLines2.append(lineInsnippet2)
            # all()
            if len(snippetLines1) == 4 and len(snippetLines2) == 4:
                snippetLineNumber = linesCount - 3
                record1, record2 = (Reads(snippetLines1, snippetLineNumber),
                                    Reads(snippetLines2, snippetLineNumber))
                snippetLines1, snippetLines2 = ([], [])
                yield record1, record2

def diffPrint(a, b):
    for i,s in enumerate(difflib.ndiff(a, b)):
        if s[0]==' ': continue
        elif s[0]=='-':
            print('Delete "{}" from position {}'.format(s[-1],i))
        elif s[0]=='+':
            print('Add "{}" to position {}'.format(s[-1],i))
    print()

ovlpSeedVH = Seeds(overlappedVHHHead[0])
ovlpSeedVE = Seeds(overlappedVHHEnd[0])
for reads1, reads2 in readR1R2():
    a = 0
    y = 0
    f = 0
    a += 1
    m1 = regex.match(ovlpSeedVH.pattern, reads1.sequence)
    if not m1:
        # m.start(0) return the whole substring
        # for 1 only return the first group
        # print m.start(1), m.end(1)
        continue
    m2 = regex.match(ovlpSeedVE.pattern, reads2.sequence)
    if m1 and m2:
        # m.start(0) return the whole substring
        # for 1 only return the first group
        # print m.start(1), m.end(1)

        print m1.start(1), m1.end(1), m2.start(1), m2.end(1)
        print "found in line:{}\nR1: {}\nR2: {}".format(reads1.lineNumber, reads1.sequence, reads2.sequence)
        # if m1.fuzzy_changes != ([],[],[]):
        #     f += 1
        #     print m1.fuzzy_changes
        # if m2.fuzzy_changes != ([],[],[]):
        #     f += 1
        #     print m2.fuzzy_changes

        print

        # diffPrint(m1.group(1), overlappedVHHHead[0])
        # diffPrint(m2.group(1), overlappedVHHEnd[0])
        y += 1

print "{} total matches in {} total reads".format(y, a)
