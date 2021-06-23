import pandas as pd
import numpy as np
import os, sys
from sys import exit
import pprint
import random
import copy
import glob
import Bio
from Bio import SeqIO

E = {'CBFB':'tCBFB', 'EloB':'tEloB', 'EloC':'tEloC', 'Vif':'tVif', 'CUL5':'tCUL5', 'Rbx2':'tRbx2'}


for p0, p1 in E.items():
    aln = []
    r0 = 1
    seq_file = '../'+p0+'.fasta'
    seq = list(SeqIO.parse(seq_file,'fasta'))[0].seq
    for rt0 in seq:
        aln.append([rt0, p0, r0, 9.0, rt0, p1, r0, 9.0])
        r0 += 1
    np.savetxt('mat_'+p0+'_'+p1+'.align', np.array(aln), fmt='%s')
