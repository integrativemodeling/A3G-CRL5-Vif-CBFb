from modeller import *
from modeller.automodel import *
import sys

# Override the 'special_restraints' and 'user_after_single_model' methods:

# Protein order: Rbx2, Cul5, EloB, EloC, CBFb, Vif
# ,'C','D','E','F'

class MyModel(loopmodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
    
    def special_patches(self, aln):
        self.rename_segments(segment_ids=('A','B','C','D','E','F'), renumber_residues=[11, 1, 17, 1, 1, 27])

    def select_atoms(self):
        # Select  residues 1 and 2 (PDB numbering)
        sele = selection(self.residue_range('119:A', '131:A'), self.residue_range('305:A', '708:A'), self.residue_range('46:C','58:C'), self.residue_range('69:D','85:D'), self.residue_range('27:F','113:F'), self.residue_range('78:B', '105:B'))
        return sele


    def select_loop_atom(self):
        return selection(self.residue_range('119:A','131:A'),
                         self.residue_range('46:C','58:C'))
    
env = environ()

env.io.hetatm = True

a = MyModel(env, alnfile='aln_CRL5.pir', knowns=('4N9F_CDEFG','1ldj','2ECL_A','2ma9_B_orient'),sequence='CRL5_mod',assess_methods=(assess.DOPE,assess.normalized_dope,assess.GA341))

a.starting_model = 1
a.ending_model = 50
if '--test' in sys.argv: a.ending_model = 1

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 2           # Last loop model

a.make()                           # do comparative modeling


