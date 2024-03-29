from modeller import *
from modeller.automodel import *
import sys

# Override the 'special_restraints' and 'user_after_single_model' methods:

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
    
    def special_patches(self, aln):
        self.rename_segments(segment_ids=('A'), renumber_residues=[6])

env = environ()

env.io.hetatm = True

a = MyModel(env, alnfile='aln_nter.pir', knowns=('5k81_E'),sequence='A3G_nter',assess_methods=(assess.DOPE,assess.GA341))

a.starting_model = 1
a.ending_model = 25
if '--test' in sys.argv: a.ending_model = 1

a.make()                           # do comparative modeling


