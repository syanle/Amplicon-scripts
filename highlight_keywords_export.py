import sys
import re
from Bio.Seq import Seq
from itertools import izip

barcodesVHH = ["TAAGGCGA", "CGTACTAG", "AGGCAGAA", "TCCTGAGC", "GGACTCCT"]
barcodesVHHEnd = ["ATGGTGCTGGCCGGCCTGGCC"]
barcodesHingeHead = ["TAGGCATG"]
barcodesHingeEnd = ["GCGATCTA", "ATAGAGAG", "TCTACTCT"]

def BarcodeJoin(barcodeList):
    return '|'.join([str(Seq(i)) for i in barcodeList])
def BarcodeRC(barcodeList):
    return '|'.join([str(Seq(i).reverse_complement()) for i in barcodeList])
def BarcodeC(barcodeList):
    return '|'.join([str(Seq(i).complement()) for i in barcodeList])
def BarcodeR(barcodeList):
    return '|'.join([str(Seq(i)[::-1]) for i in barcodeList])


def colorize(text):
    # colorize the barcodes in one reads (global greedy match)
    # highlight VHH barcodes
    # red normal order
    text = re.sub(BarcodeJoin(barcodesVHH),lambda m: '\x1b[1;30;41m{}\x1b[0m'.format(m.group()), text)
    # green reverse_complement
    text = re.sub(BarcodeRC(barcodesVHH),lambda m: '\x1b[1;30;42m{}\x1b[0m'.format(m.group()), text)
    # yellow reverse order
    text = re.sub(BarcodeR(barcodesVHH),lambda m: '\x1b[1;30;43m{}\x1b[0m'.format(m.group()), text)
    # white complemented
    text = re.sub(BarcodeC(barcodesVHH),lambda m: '\x1b[1;30;44m{}\x1b[0m'.format(m.group()), text)

    # highlight VHH 3' overlapped primer
    text = re.sub(BarcodeJoin(barcodesVHHEnd),lambda m: '\x1b[7;33;41m{}\x1b[0m'.format(m.group()), text)
    # green reverse_complement
    text = re.sub(BarcodeRC(barcodesVHHEnd),lambda m: '\x1b[7;33;42m{}\x1b[0m'.format(m.group()), text)
    # yellow reverse order
    text = re.sub(BarcodeR(barcodesVHHEnd),lambda m: '\x1b[7;33;43m{}\x1b[0m'.format(m.group()), text)
    # white complemented
    text = re.sub(BarcodeC(barcodesVHHEnd),lambda m: '\x1b[7;33;44m{}\x1b[0m'.format(m.group()), text)

    # highlight Hinge head barcodes
    text = re.sub(BarcodeJoin(barcodesHingeHead),lambda m: '\x1b[7;37;41m{}\x1b[0m'.format(m.group()), text)
    # green reverse_complement
    text = re.sub(BarcodeRC(barcodesHingeHead),lambda m: '\x1b[7;37;42m{}\x1b[0m'.format(m.group()), text)
    # yellow reverse order
    text = re.sub(BarcodeR(barcodesHingeHead),lambda m: '\x1b[7;37;43m{}\x1b[0m'.format(m.group()), text)
    # white complemented
    text = re.sub(BarcodeC(barcodesHingeHead),lambda m: '\x1b[7;37;44m{}\x1b[0m'.format(m.group()), text)

    # highlight Hinge tail barcodes
    text = re.sub(BarcodeJoin(barcodesHingeEnd),lambda m: '\x1b[7;30;41m{}\x1b[0m'.format(m.group()), text)
    # green reverse_complement
    text = re.sub(BarcodeRC(barcodesHingeEnd),lambda m: '\x1b[7;30;42m{}\x1b[0m'.format(m.group()), text)
    # yellow reverse order
    text = re.sub(BarcodeR(barcodesHingeEnd),lambda m: '\x1b[7;30;43m{}\x1b[0m'.format(m.group()), text)
    # white complemented
    return re.sub(BarcodeC(barcodesHingeEnd),lambda m: '\x1b[7;30;44m{}\x1b[0m'.format(m.group()), text)


class MyStdout(object):
    def __init__(self, term=sys.stdout):
        self.term = term
    def write(self, text):
        text = colorize(text)
        self.term.write(text)

sys.stdout = MyStdout()

with open('export.txt','w') as outFile:

    with open('../chip_R1.fq') as f1, open('../chip_R2.fq') as f2:
        lineNumber = 0
        for r1, r2 in izip(f1, f2):
            lineNumber += 1
            if (lineNumber - 2) % 4 == 0:
                outFile.write("R1 line:{2}\n{0}\nR2 line:{2}\n{1}\n\n".format(r1.strip(), r2.strip(), lineNumber))

                if lineNumber > 3000 * 4:
                    break
