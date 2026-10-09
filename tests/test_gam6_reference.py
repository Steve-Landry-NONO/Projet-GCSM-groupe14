"""Vérifier que le premier rendu signale les données manquantes."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from gam6_reference import compute_reference

class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.compounds = pd.DataFrame([['cible',128,127,1]],
            columns=['name','mz_quant','mz_qual','elution_rank'])
        t=np.linspace(1,3,1001)
        y=10+1000*np.exp(-((t-2)/.03)**2)
        self.traces={128:(t,y),127:(t,10+.2*(y-10))}
    def test_quantifier_missing(self):
        d=compute_reference({},self.compounds)
        self.assertEqual(len(d),1)
        self.assertIn('quantifiant',d.iloc[0].warning)
    def test_qualifier_missing(self):
        d=compute_reference({128:self.traces[128]},self.compounds)
        self.assertIn('ratio non calculable',d.iloc[0].warning)
    def test_nonoverlapping_qualifier(self):
        t,y=self.traces[127];self.traces[127]=(t+10,y)
        d=compute_reference(self.traces,self.compounds)
        self.assertIn('fenêtre qualifiante',d.iloc[0].warning)
    def test_isolated_peak(self):
        d=compute_reference(self.traces,self.compounds)
        self.assertAlmostEqual(d.iloc[0].ratio_ref,.2,places=3)
    def test_no_compounds(self):
        d=compute_reference(self.traces,self.compounds.iloc[:0])
        self.assertTrue(d.empty)

if __name__=='__main__':
    unittest.main()
