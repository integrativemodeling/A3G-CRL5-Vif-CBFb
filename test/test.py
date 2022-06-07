#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def run_modeller_script(self, script_name, model_name, resrng):
        """Run a Modeller script and test the output model"""
        # Run script
        p = subprocess.check_call(["python", script_name, "--test"])
        # Make sure PDB was produced with the requested residue range

        with open('%s.B99990001.pdb' % model_name) as fh:
            pdb_lines = [x for x in fh.readlines() if x.startswith('ATOM')]
        rng = (int(pdb_lines[0][22:26]), int(pdb_lines[-1][22:26]))
        self.assertEqual(rng, resrng)

    def test_model_cter(self):
        """Test generation of full model for C-terminal domain using Modeller"""
        os.chdir(os.path.join(TOPDIR, 'data'))
        self.run_modeller_script('mod_cter.py', 'A3G_cter', (195, 380))

    def test_model_nter(self):
        """Test generation of full model for N-terminal domain using Modeller"""
        os.chdir(os.path.join(TOPDIR, 'data'))
        self.run_modeller_script('mod_nter.py', 'A3G_nter', (6, 195))

    # We don't test mod_CRL5.py currently since it takes forever

    def test_simple_rigid(self):
        """Test model building, rigid"""
        os.chdir(os.path.join(TOPDIR, 'scripts'))
        p = subprocess.check_call(["python", "mod_A3G_Vif_CRL5_rigid.py",
                                   "--test"])
        # todo: assert outputs, run analysis

    def test_simple_flexible(self):
        """Test model building, flexible"""
        os.chdir(os.path.join(TOPDIR, 'scripts'))
        p = subprocess.check_call(["python", "mod_A3G_Vif_CRL5_flexible.py",
                                   "--test"])
        # todo: assert outputs, run analysis

    def test_mmcif_rigid(self):
        """Test generation of mmCIF output, rigid"""
        os.chdir(os.path.join(TOPDIR, 'scripts'))
        if os.path.exists("A3G-CRL5-Vif_complex_rigid.cif"):
            os.unlink("A3G-CRL5-Vif_complex_rigid.cif")
        p = subprocess.check_call(
                ["python", "mod_A3G_Vif_CRL5_rigid.py", "--mmcif",
                 "--dry-run"])
        # Check output file
        self._check_mmcif_file('A3G-CRL5-Vif_complex_rigid.cif')

    def test_mmcif_flexible(self):
        """Test generation of mmCIF output, flexible"""
        os.chdir(os.path.join(TOPDIR, 'scripts'))
        if os.path.exists("A3G-CRL5-Vif_complex_rigid.cif"):
            os.unlink("A3G-CRL5-Vif_complex_flexible.cif")
        p = subprocess.check_call(
                ["python", "mod_A3G_Vif_CRL5_flexible.py", "--mmcif",
                 "--dry-run"])
        # Check output file
        self._check_mmcif_file('A3G-CRL5-Vif_complex_flexible.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 3)
        self.assertEqual(s.citations[0].doi, '10.1038/s41594-018-0118-5')
        self.assertEqual(len(s.software), 2)
        self.assertEqual(len(s.orphan_starting_models), 12)
        # Should be 1 state
        self.assertEqual(len(s.state_groups), 1)
        state1, = s.state_groups[0]
        # Should be 1 model
        self.assertEqual(sum(len(x) for x in state1), 1)
        # Check # of spheres and atoms in each model
        m = state1[0][0]
        self.assertEqual(len(m._spheres), 3640)
        self.assertEqual(len(m._atoms), 0)
        # Should be 1 ensemble
        self.assertEqual([e.num_models for e in s.ensembles], [1640])
        # Just one restraint - crosslinks
        xl, = s.restraints
        self.assertEqual(len(xl.experimental_cross_links), 40)
        self.assertEqual(len(xl.cross_links), 40)


if __name__ == '__main__':
    unittest.main()
