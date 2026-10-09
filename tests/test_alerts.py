"""Classement des alertes et rapport vallée/hauteur.

Ces tests vérifient des propriétés de sécurité du critère vallée/hauteur sur
une grille de cas synthétiques (proportions de pics, écarts, bruit, lignes de
base). Ils ne prouvent pas que les aires sont justes sous le seuil : la
caractérisation chiffrée est dans docs/VALLEE_HAUTEUR.md.
"""
import itertools
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from alerts import Alert, blocking_from_reference, peak_bound_alerts, render  # noqa: E402
from measure import measure_sample  # noqa: E402
from peaks import Peak, bound_level_ratio  # noqa: E402
from valley_cases import GRID, run_case  # noqa: E402


def grid(seeds=(0,)):
    rows = []
    for hr, sep, noise, base, seed in itertools.product(*GRID.values(), seeds):
        rows += run_case(hr, sep, noise, base, seed)
    d = pd.DataFrame(rows)
    return d[d.detected]


def blocked_grid(seeds=(0,)):
    """Grille avec, pour chaque pic mesuré, le verdict des alertes de bornes."""
    rows = []
    for hr, sep, noise, base, seed in itertools.product(*GRID.values(), seeds):
        rows += run_case(hr, sep, noise, base, seed, with_peaks=True)
    d = pd.DataFrame(rows)
    d = d[d.detected].copy()
    d["blocked"] = [any(a.blocking for a in peak_bound_alerts(p)) for p in d.pop("peak_obj")]
    return d


class LevelRatioTests(unittest.TestCase):
    def test_definition(self):
        self.assertAlmostEqual(bound_level_ratio(55, 110, 10), 0.45)
        self.assertEqual(bound_level_ratio(5, 110, 10), 0.0)      # sous le fond : 0
        self.assertEqual(bound_level_ratio(200, 110, 10), 1.0)    # au-dessus : 1
        self.assertTrue(np.isnan(bound_level_ratio(5, 10, 10)))   # hauteur nulle


class ValleyGridTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = blocked_grid()

    def test_precise_at_low_noise_flat_baseline(self):
        s = self.d[(self.d.noise == 0.001) & (self.d.base == "plate")]
        self.assertGreater(len(s), 20)
        self.assertLessEqual((s.v_est - s.v_ref).abs().max(), 0.03)

    def test_no_marked_overlap_becomes_informative(self):
        # Quel que soit le bruit ou la ligne de base, un creux réel > 15 % de la
        # hauteur ne doit jamais passer sous le seuil de 10 %.
        bad = self.d[(self.d.v_ref > 0.15) & (self.d.v_est <= 0.10)]
        self.assertEqual(len(bad), 0, bad.to_string())

    def test_strong_overlap_always_blocking(self):
        # Chevauchement marqué : soit la vallée est mesurée au-dessus du seuil, soit
        # une borne reste au-dessus du fond (épaulement) ; jamais aucune alerte.
        s = self.d[self.d.v_ref > 0.30]
        self.assertTrue(s.blocked.all(), s[~s.blocked].to_string())

    def test_decreases_with_separation(self):
        s = self.d[(self.d.noise == 0.001) & (self.d.base == "plate") & (self.d.height_ratio == 1.0)]
        for peak in (1, 2):
            v = s[s.peak == peak].sort_values("sep").v_est.to_numpy()
            self.assertTrue(np.all(np.diff(v) < 0), v)


class AlertClassificationTests(unittest.TestCase):
    def peak(self, left="baseline", right="baseline", lr=np.nan, rr=np.nan):
        return Peak(1.0, 0.9, 1.1, 10.0, 100.0, 0.0, left, right, 50.0, 0.0, lr, rr)

    def test_valley_threshold_is_configurable(self):
        p = self.peak(right="vallee", rr=0.08)
        self.assertFalse(peak_bound_alerts(p)[0].blocking)
        self.assertTrue(peak_bound_alerts(p, valley_max=0.05)[0].blocking)

    def test_unmeasurable_valley_blocks(self):
        self.assertTrue(peak_bound_alerts(self.peak(right="vallee"))[0].blocking)

    def test_sim_edge_low_residual_informative_high_blocking(self):
        self.assertFalse(peak_bound_alerts(self.peak(right="fin_signal", rr=0.0003))[0].blocking)
        self.assertTrue(peak_bound_alerts(self.peak(right="fin_signal", rr=0.05))[0].blocking)
        self.assertTrue(peak_bound_alerts(self.peak(right="fin_signal"))[0].blocking)

    def test_max_width_blocks(self):
        self.assertTrue(peak_bound_alerts(self.peak(left="limite"))[0].blocking)

    def test_render_separates_severity(self):
        out = render([Alert("a"), Alert("b", blocking=False), Alert("a")])
        self.assertEqual(out["blocking_alerts"], "a")
        self.assertEqual(out["info_alerts"], "b")
        self.assertEqual(out["warning"], "a ; b")

    def test_only_blocking_reference_alerts_propagate(self):
        row = pd.Series({"blocking_alerts": "autre pic proche en intensité sur ce m/z : attribution ambiguë",
                         "info_alerts": "pic voisin séparé : vallée 1.0 % de la hauteur (<= 10 %)",
                         "warning": "x"})
        got = blocking_from_reference(row)
        self.assertEqual(len(got), 1)
        self.assertTrue(got[0].blocking and "attribution ambiguë" in got[0].message)

    def test_old_reference_without_severity_is_conservative(self):
        got = blocking_from_reference(pd.Series({"warning": "pic voisin : séparation des aires à valider"}))
        self.assertTrue(got and got[0].blocking)


class SampleStatusTests(unittest.TestCase):
    def setUp(self):
        t = np.linspace(1, 3, 1001)
        y = 10 + 1000 * np.exp(-((t - 2) / .03) ** 2)
        self.traces = {128: (t, y), 127: (t, 10 + .2 * (y - 10)),
                       136: (t, y), 135: (t, 10 + .2 * (y - 10))}
        self.compounds = pd.DataFrame([['cible', 128, 127, 1, 'étalon'], ['étalon', 136, 135, 2, '']],
                                      columns=['name', 'mz_quant', 'mz_qual', 'elution_rank', 'istd'])
        self.ref = pd.DataFrame({'name': ['cible', 'étalon'], 'mz_quant': [128, 136], 'mz_qual': [127, 135],
                                 'rt_ref': [2., 2.], 'ratio_ref': [.2, .2],
                                 'warning': ['', ''], 'blocking_alerts': ['', ''], 'info_alerts': ['', '']})

    def test_informative_reference_alert_does_not_block(self):
        self.ref.loc[0, ['warning', 'info_alerts']] = 'pic voisin séparé : vallée 1.0 % de la hauteur (<= 10 %)'
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertEqual(d.loc['cible', 'status'], 'OK')
        self.assertAlmostEqual(d.loc['cible', 'response'], 1.)

    def test_ambiguous_attribution_stays_blocking(self):
        msg = 'autre pic proche en intensité sur ce m/z : attribution ambiguë'
        self.ref.loc[0, ['warning', 'blocking_alerts']] = msg
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertEqual(d.loc['cible', 'status'], 'à vérifier')
        self.assertTrue(np.isnan(d.loc['cible', 'response']))

    def test_double_peak_and_despike_block(self):
        for col in ('double_peak', 'median_despike'):
            c = self.compounds.assign(**{col: [True if col == 'double_peak' else 1, np.nan]})
            d = measure_sample(self.traces, self.ref, c).set_index('name')
            self.assertEqual(d.loc['cible', 'status'], 'à vérifier', col)

    def test_invalid_istd_is_named(self):
        self.ref.loc[1, ['warning', 'blocking_alerts']] = 'ordre d\'élution non respecté'
        d = measure_sample(self.traces, self.ref, self.compounds).set_index('name')
        self.assertIn('étalon interne non validé', d.loc['cible', 'blocking_alerts'])


if __name__ == "__main__":
    unittest.main()
