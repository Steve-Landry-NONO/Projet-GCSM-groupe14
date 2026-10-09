"""Petits exemples pour vérifier les erreurs qui faussaient les résultats.

On peut lancer ces tests sans les fichiers du laboratoire.
"""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from measure import measure_sample
from peaks import integrate_peak
from reference import compute_reference
from run_calibration import invert, sf_verdict


class ValidationTests(unittest.TestCase):
    def setUp(self):
        t = np.linspace(1, 3, 1001)
        y = 10 + 1000 * np.exp(-((t - 2) / .03)**2)
        self.traces = {128: (t, y), 127: (t, 10 + .2 * (y - 10)),
                       136: (t, y), 135: (t, 10 + .2 * (y - 10))}
        self.compounds = pd.DataFrame([
            ['cible', 128, 127, 1, 'étalon'], ['étalon', 136, 135, 2, '']],
            columns=['name', 'mz_quant', 'mz_qual', 'elution_rank', 'istd'])
        self.ref = pd.DataFrame({'name': ['cible', 'étalon'], 'mz_quant': [128, 136],
                                 'mz_qual': [127, 135], 'rt_ref': [2., 2.],
                                 'ratio_ref': [.2, .2], 'warning': ['', '']})

    def test_below_range(self):
        value, status = invert([0, 1, 0], .02, 5, .025)
        self.assertAlmostEqual(value, .02)
        self.assertEqual(status, 'hors gamme')

    def test_domain_and_bounds(self):
        for value in [.025, 1., 5.]:
            self.assertEqual(invert([0, 1, 0], value, 5, .025)[1], 'ok')
        self.assertEqual(invert([0, 1, 0], 6, 5, .025)[1], 'hors gamme')

    def test_nonmonotone_and_constant(self):
        for coef in [[1, -2, 0], [0, 0, 1], [-1, 1, 0]]:
            self.assertNotEqual(invert(coef, 1, 5, .025)[1], 'ok')

    def test_positive_slope_root(self):
        # La plus petite racine peut être sur la branche décroissante.
        value, status = invert([1, -2, 0], 3, 5, 2)
        self.assertAlmostEqual(value, 3)
        self.assertEqual(status, 'ok')

    def test_sf_requires_whole_validation(self):
        self.assertEqual(sf_verdict('OK', 'ok', 2, 20), 'PASS')
        self.assertEqual(sf_verdict('OK', 'ok', 30, 20), 'FAIL')
        for peak, inversion, err in [('à vérifier', 'ok', 2), ('OK', 'hors gamme', 2),
                                     ('OK', 'ok', np.nan)]:
            self.assertEqual(sf_verdict(peak, inversion, err, 20), 'NON VALIDÉ')

    def test_valid_normalization(self):
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertEqual(d.loc['cible', 'status'], 'OK')
        self.assertAlmostEqual(d.loc['cible', 'response'], 1.)

    def test_invalid_istd_blocks_response(self):
        self.ref.loc[self.ref.name == 'étalon', 'ratio_ref'] = .8
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertTrue(np.isnan(d.loc['cible', 'response']))
        self.assertNotEqual(d.loc['cible', 'status'], 'OK')

    def test_reference_warning_blocks_response(self):
        self.ref.loc[0, 'warning'] = 'affectation ambiguë'
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertTrue(np.isnan(d.loc['cible', 'response']))
        self.assertEqual(d.loc['cible', 'status'], 'à vérifier')

    def test_missing_reference_keeps_expected_compound(self):
        self.ref.loc[0, 'rt_ref'] = np.nan
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertEqual(len(d), 2)
        self.assertEqual(d.loc['cible', 'status'], 'référence absente')

    def test_no_detected_ions_does_not_crash(self):
        d = measure_sample({}, self.ref, self.compounds)
        self.assertEqual(len(d), 2)
        self.assertTrue(d.response.isna().all())

    def test_no_reference_ions_has_stable_columns(self):
        d = compute_reference({}, self.compounds)
        self.assertEqual(len(d), 2)
        self.assertTrue(d.rt_ref.isna().all())

    def test_qualifier_window_does_not_overlap(self):
        t, y = self.traces[127]
        self.traces[127] = (t + 10, y)
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertIn('fenêtre qualifiante', d.loc['cible', 'warning'])
        self.assertTrue(np.isnan(d.loc['cible', 'response']))

    def test_metadata_applicability(self):
        self.compounds['applicable_to'] = ['GAM,BLPC', np.nan]
        d = measure_sample(self.traces, self.ref, self.compounds, sample_type='SF')
        self.assertEqual(d.iloc[0].status, 'non applicable')

    def test_duplicate_assignment_blocks_both(self):
        self.compounds.loc[1, ['mz_quant', 'mz_qual']] = [128, 127]
        self.ref.loc[1, ['mz_quant', 'mz_qual']] = [128, 127]
        d = measure_sample(self.traces, self.ref, self.compounds)
        self.assertTrue(d.status.eq('à vérifier').all())
        self.assertTrue(d.response.isna().all())

    def test_bad_trace_rejected(self):
        with self.assertRaises(ValueError):
            integrate_peak(np.array([1, 1, 2]), np.ones(3), 1)

    def test_multiple_batches_rejected_before_reading(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ['batch1', 'batch2']:
                (Path(tmp) / name / 'GAM-25-385-6.D').mkdir(parents=True)
            result = subprocess.run([sys.executable, str(ROOT / 'src/run_batch.py'), tmp],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Plusieurs batchs', result.stderr)

    def test_invalid_points_do_not_enter_calibration(self):
        with tempfile.TemporaryDirectory() as tmp:
            batch = Path(tmp) / 'batch'
            batch.mkdir()
            d = pd.DataFrame([dict(name='cible', istd='étalon', sample=f'GAM-{i}',
                type='GAM', level=i, conc_nominal_ppm=x, response=x,
                status='à vérifier', warning='pic douteux') for i, x in enumerate([.025,.05,.1,1],1)])
            d.to_csv(batch / 'peaks.csv', index=False)
            result = subprocess.run([sys.executable, str(ROOT / 'src/run_calibration.py'),
                                     '--out', tmp], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            cal = pd.read_csv(batch / 'calibrations.csv')
            pts = pd.read_csv(batch / 'calibration_points.csv')
            self.assertNotEqual(cal.iloc[0].status, 'ok')
            self.assertFalse(pts.included.any())

    def test_exploratory_keeps_alert_and_strict_block(self):
        self.ref.loc[0, 'warning'] = 'vallée à examiner'
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertTrue(np.isnan(d.loc['cible', 'response']))
        self.assertAlmostEqual(d.loc['cible', 'response_exploratory'], 1.)
        self.assertEqual(d.loc['cible', 'status'], 'à vérifier')

    def test_exploratory_calibration_never_validates_sf(self):
        with tempfile.TemporaryDirectory() as tmp:
            batch = Path(tmp) / 'batch'; batch.mkdir()
            rows = [dict(name='cible', istd='étalon', sample=f'GAM-{i}', type='GAM',
                level=i, conc_nominal_ppm=x, response=np.nan, response_exploratory=x,
                status='à vérifier', warning='vallée', identity_ok=True, istd_identity_ok=True)
                for i,x in enumerate([.025,.05,.1,1],1)]
            rows.append(dict(rows[-1], sample='SF-1', type='SF'))
            pd.DataFrame(rows).to_csv(batch / 'peaks.csv', index=False)
            run = subprocess.run([sys.executable,str(ROOT/'src/run_calibration.py'),
                '--out',tmp,'--exploratory'],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            cal=pd.read_csv(batch/'calibrations.csv'); sf=pd.read_csv(batch/'sf_results.csv')
            self.assertEqual(cal.iloc[0].status,'exploratoire')
            self.assertEqual(sf.iloc[0].sf_status,'NON VALIDÉ (exploratoire)')
            self.assertAlmostEqual(sf.iloc[0].conc_calc_ppm,1.)


if __name__ == '__main__':
    unittest.main()
