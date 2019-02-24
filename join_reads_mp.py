#!/usr/bin/env python2
from itertools import izip, islice
import regex, os
import difflib
from Bio import pairwise2
from multiprocessing import Process, Pool

os.chdir("/home/xue/Documents/pe150/amplicon_ngs")

# spacer added to VHH barcodes
barcodes_VHH = ["TAAGGCGA", "CGTACTAGT",
                "AGGCAGAATA", "TCCTGAGCATG", "GGACTCCTATGC"]
overlapped_VHH_head = ["ACCGTGGCCCAGGCGGCCCAG"]
overlapped_VHH_tail = ["ATGGTGCTGGCCGGCCTGGCC"]
barcodes_hinge_head = ["TAGGCATG"]
# spacer added
barcodes_hinge_tail = ["GCGATCTATA", "ATAGAGAGT", "TCTACTCT"]


class Reads:

    # return a dict only records one read info, therefore no use of class
    # def snippetPacker(lines=None):
    #     ks = ['name', 'sequence', 'optional', 'quality']
    #     return {k: v for k, v in zip(ks, lines)}

    def __init__(self, lines, line_number):
        self.name, input_sequence, self.optional, self.quality = lines
        self.sequence = input_sequence.strip()
        self.line_number = line_number
        self.ovlp_match = None
        self.lib_match = None

    # _ = reads1 & reads2
    def __and__(self, other):
        self._matched_items = self.ovlp_match.groupdict()
        other._matched_items = other.ovlp_match.groupdict()

        for key, _ in zip(self._matched_items, other._matched_items):
            if self._matched_items[key] and other._matched_items[key]:
                return True
        return False


class BaseSeeds(dict):

    def __new__(cls, *args, **kwargs):
        if cls is BaseSeeds:
            raise TypeError("base class may not be instantiated")
        return dict.__new__(cls, *args, **kwargs)

    def __init__(self, **kwargs):
        # act as a dict
        super(BaseSeeds, self).__init__(**kwargs)
        # automatic assign input of the dict to properties
        self.__dict__.update(kwargs)
        self.errors = 1


class OvlpSeeds(BaseSeeds):

    def __init__(self, seq_dict, ovlp_errors=None):
        # BaseSeeds.__init__(self, **seq_dict)
        super(OvlpSeeds, self).__init__(**seq_dict)
        if ovlp_errors:
            self.errors = ovlp_errors

    @property
    def pattern(self):
        # keys in dict is shuffle so I have to hardcode it
        return "(?e)(?P<head>{0}){{e<={2}}}|(?e)(?P<tail>{1}){{e<={2}}}".format(self.head, self.tail, self.errors)


class LibsSeeds(BaseSeeds):

    def __init__(self, seq_dict, lib_errors=None):
        self.seq_dict = seq_dict
        # BaseSeeds.__init__(self, **seq_dict)
        super(LibsSeeds, self).__init__(**seq_dict)
        if lib_errors:
            self.errors = lib_errors

    @property
    def pattern(self):
        # keys in dict is shuffle so I have to hardcode it
        # TODO: fix hardcode
        return ".*(?e)(?P<lib0>{0}){{e<={5}}}|"\
               "(?e)(?P<lib1>{1}){{e<={5}}}|"\
               "(?e)(?P<lib2>{2}){{e<={5}}}|"\
               "(?e)(?P<lib3>{3}){{e<={5}}}|"\
               "(?e)(?P<lib4>{4}){{e<={5}}}".format(
                   self.lib0, self.lib1, self.lib2, self.lib3, self.lib4, self.errors)

# from Bio.SeqIO.QualityIO import FastqGeneralIterator


def read_R1R2():
    with open('Cleandata/chip/chip_R1.fq') as f1, open('Cleandata/chip/chip_R2.fq') as f2:
        # lines1 = islice(f1, 1, None, 4)
        # lines2 = islice(f2, 1, None, 4)
        snippet_lines1, snippet_lines2 = ([], [])
        lines_count = 0
        for line_in_snippet1, line_in_snippet2 in izip(f1, f2):
            lines_count += 1
            snippet_lines1.append(line_in_snippet1)
            snippet_lines2.append(line_in_snippet2)
            # all()
            if len(snippet_lines1) == 4 and len(snippet_lines2) == 4:
                snippet_line_number = lines_count - 3
                record1, record2 = (Reads(snippet_lines1, snippet_line_number),
                                    Reads(snippet_lines2, snippet_line_number))
                snippet_lines1, snippet_lines2 = ([], [])
                yield record1, record2

overlapped_VHH = {"head": "ACCGTGGCCCAGGCGGCCCAG",
                  "tail": "ATGGTGCTGGCCGGCCTGGCC"}
lib_VHH = {"lib0": "TAAGGCGA",
           "lib1": "CGTACTAGT",
           "lib2": "AGGCAGAATA",
           "lib3": "TCCTGAGCATG",
           "lib4": "GGACTCCTATGC"}

ovlp_seed = OvlpSeeds(overlapped_VHH)
lib_seed = LibsSeeds(lib_VHH)


def reads_hunter(record_pair):
    for record in record_pair:
        ovlp_match = regex.search(ovlp_seed.pattern, record.sequence)
        if ovlp_match:
            record.ovlp_match = ovlp_match
        else:
            pass
            # continue

    record1, record2 = record_pair
    # hit simultaneously in both reads1 and reads2
    if record1.ovlp_match and record2.ovlp_match:
        # if the record1 and record2 both have head hit or tail hit, discard
        # them
        # mr += 1
        if record1 & record2:
            # br += 1
            print "line:", record1.line_number, record1.ovlp_match.groupdict(), record2.ovlp_match.groupdict()
            # TODO: can not continue in standalone process
            # continue
            record1.ovlp_match = record2.ovlp_match = None


        # print "found in line:", record1.line_number, ":", record2.line_number
        # print "R1:", record1.sequence
        # print "R2:", record2.sequence
        # print
    # xor
    elif record1.ovlp_match != record2.ovlp_match:
        # smr += 1
        pass

    # ovlp_seed error = 1: found 13337 bad matches in 1410301 matches (while 1372695 single matches failed) in 11044350 total records
    # ovlp_seed error = 2: found 65672 bad matches in 1716225 matches (while
    # 1526712 single matches failed) in 11044350 total records

    # print "found {} bad matches in {} matches (while {} single matches failed) in {} total records".format(br, mr, smr, tr)

if __name__ == '__main__':

    p = Pool(processes=40)
    record_handle = read_R1R2()
    print "handler loaded"
    p.map(reads_hunter, record_handle, chunksize=100)
    # p.map             real 4m24.102s
    # p.map_async       real 4m42.837s
    # p.imap            real 7m26.276s
    # p.imap_unordered  real 7m30.898s

    p.close()
    p.join()
