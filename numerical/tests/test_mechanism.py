import unittest,csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def rows(n):
    with (ROOT/'common_calibration/results'/n).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
class MechanismEvidence(unittest.TestCase):
    def test_affine_limit_all_grid_states(self):
        r=[r for r in rows('dissolved_grid.csv') if float(r['dissolved_coefficient'])==0]
        self.assertEqual(len(r),882)
        self.assertLess(max(float(x['absolute_difference']) for x in r),1e-12)
    def test_new_signed_map_matches_original_pair_grid(self):
        expected={(r['match'],float(r['Pv']),float(r['Q'])):float(r['absolute_difference']) for r in rows('discrimination_grid.csv') if r['law1']=='saturation' and r['law2']=='constant_D'}
        actual=rows('signed_state_map.csv');self.assertEqual(len(actual),len(expected))
        self.assertLess(max(abs(float(r['absolute_difference'])-expected[(r['match'],float(r['Pv']),float(r['Q']))]) for r in actual),1e-10)
    def test_calibration_target_for_each_ablation(self):
        for r in rows('dissolved_calibration.csv'):
            if r['match']=='endpoint':self.assertAlmostEqual(float(r['reference_content']),float(r['reference_saturation']),places=11)
    def test_calibration_changes_residual_sign_at_ranked_max(self):
        r={x['match']:x for x in rows('ranked_grid_states.csv') if x['rank']=='1'}
        self.assertLess(float(r['endpoint']['signed_difference']),0)
        self.assertGreater(float(r['inlet_slope']['signed_difference']),0)
        self.assertEqual(r['endpoint']['Pv'],'30')
        self.assertAlmostEqual(float(r['endpoint']['Q']),5)
if __name__=='__main__':unittest.main()
