import importlib.util
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = ROOT / "code" / "build_lung_model.py"
PHASE2_PATH = ROOT / "code" / "phase2_scans.py"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lung = load_module("lung_test", MODEL_PATH)
phase2 = load_module("phase2_test", PHASE2_PATH)


class ModelVerificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case = lung.build_cases()

    def test_flow_and_pvr_conservation(self):
        for multiplier in (1.0, 5.0, 8.4, 22.0):
            result = lung.solve(self.case, arterial_multiplier=multiplier)
            self.assertLess(abs(result["flow_conservation_error_L_min"]), 1e-9)
            self.assertLess(abs(result["pvr_identity_error_WU"]), 1e-12)
            self.assertAlmostEqual(sum(result["flow_fraction"].values()), 1.0, places=12)

    def test_arterial_only_primary_operator(self):
        result = lung.solve(self.case, arterial_multiplier=5.0)
        self.assertEqual(result["arterial_multiplier"], 5.0)
        self.assertEqual(result["venous_multiplier"], 1.0)
        self.assertGreater(result["PVR_WU"], lung.solve(self.case)["PVR_WU"])

    def test_closed_oxygen_mass_balance(self):
        result = lung.solve(self.case, arterial_multiplier=8.4)
        ox = lung.oxygen(result)
        self.assertTrue(ox["converged"])
        self.assertEqual(ox["mixing_basis"], "oxygen_content")
        self.assertLess(abs(ox["oxygen_mass_balance_error_mL_min"]), 1e-6)
        self.assertAlmostEqual(
            ox["total_VO2_demand_mL_min"] - ox["total_VO2_used_mL_min"],
            ox["total_deficit_mL_min"],
            places=9,
        )

    def test_target_pvr_solver(self):
        multiplier = phase2.arterial_multiplier_for(self.case, 4.19)
        result = lung.solve(self.case, arterial_multiplier=multiplier)
        self.assertAlmostEqual(result["PVR_WU"], 4.19, places=10)

    def test_controlled_focal_pair_matches_pvr_and_volume(self):
        target = 4.19
        uniform = lung.solve(
            self.case,
            arterial_multiplier=phase2.arterial_multiplier_for(self.case, target),
        )
        focal = phase2.formal_focal_match(self.case, target)["result"]
        focal_vc = 86.0 * sum(
            lung.LOBE_W[lobe] * focal["open_fraction"][lobe]
            for lobe in lung.LOBE_ORDER
        )
        uniform_ox = lung.oxygen(
            uniform, vc_ml=focal_vc, tau_eq=1.0, vc_mode="anatomy"
        )
        focal_ox = lung.oxygen(
            focal, vc_ml=86.0, tau_eq=1.0, vc_mode="anatomy"
        )
        self.assertAlmostEqual(uniform["PVR_WU"], focal["PVR_WU"], places=8)
        self.assertAlmostEqual(
            uniform_ox["effective_vc_ml"], focal_ox["effective_vc_ml"], places=8
        )
        self.assertGreater(focal_ox["S_end_range"], uniform_ox["S_end_range"])

    def test_oxygen_content_inverse(self):
        for po2 in (20.0, 40.0, 80.0, 100.0):
            recovered = lung.po2_from_content(float(lung.oxygen_content(po2)))
            self.assertTrue(math.isclose(recovered, po2, abs_tol=0.02))

    def test_latin_hypercube_is_reproducible_and_stratified(self):
        first = phase2.latin_hypercube(50, 4, seed=20261008)
        second = phase2.latin_hypercube(50, 4, seed=20261008)
        self.assertTrue((first == second).all())
        self.assertTrue(((first >= 0.0) & (first < 1.0)).all())
        for column in range(first.shape[1]):
            occupied = sorted((first[:, column] * 50).astype(int).tolist())
            self.assertEqual(occupied, list(range(50)))


if __name__ == "__main__":
    unittest.main()
