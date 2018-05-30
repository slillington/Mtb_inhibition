#biopython pairwise alignment test

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

matrix = matlist.blosum62
a = pairwise2.align.globaldx("KEVLA", "EVL", matrix)[0]
print(matrix)
print(a[2])

