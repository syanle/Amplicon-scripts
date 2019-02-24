from Bio.Seq import Seq
barcodesVHH = ["TAAGGCGA", "CGTACTAG", "AGGCAGAA", "TCCTGAGC", "GGACTCCT"]
print '|'.join([str(Seq(i).reverse_complement()) for i in barcodesVHH])
