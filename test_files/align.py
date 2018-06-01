#biopython pairwise alignment test

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd

matrix = matlist.blosum62
a = pairwise2.align.globaldx("KEVLA", "EVL", matrix)[0]
b = {'result':a}
c = pd.DataFrame(b)
print(c)

