import unittest
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from microcirculation import (Coupling, solve_bed, oxygen_content,
    end_capillary_content, mix_beds, OBSERVATION_MATRIX, HYDRAULIC_ROW, NULL_DIRECTION)
from analyze import fit_model, predict, load_groups


class MechanismChecks(unittest.TestCase):
    def coupling(self):
        # Synthetic fixtures; none of these are estimated human parameters.
        return Coupling(12.,85.,30.,.15,.3,.02,.2,.5)

    def test_pressure_and_flow_balance(self):
        result=solve_bed(self.coupling(),25.,8.,1.,.2)
        self.assertLess(abs(result['pressure_residual']),1e-9)
        self.assertAlmostEqual(result['flow']*result['total_resistance'],17.,places=9)
        self.assertTrue(8 < result['capillary_pressure'] < 25)

    def test_reserve_is_finite(self):
        c=self.coupling();state=c.equilibrium(1e8)
        self.assertLessEqual(state[0],2.)

    def test_stationary_data_cannot_select_time_constants(self):
        c=self.coupling();target=c.equilibrium(20.)
        for tau in ([.1,.2,.3],[30.,40.,50.]):
            np.testing.assert_allclose(c.step(target,20.,5.,tau),target)
        near=c.step([1,1,1],20.,1e4,[.1,.2,.3])
        np.testing.assert_allclose(near,target)

    def test_geometry_null_changes_hydraulics(self):
        np.testing.assert_allclose(OBSERVATION_MATRIX@NULL_DIRECTION,0.)
        self.assertEqual(np.linalg.matrix_rank(OBSERVATION_MATRIX),2)
        self.assertEqual(HYDRAULIC_ROW@NULL_DIRECTION,2.)
        self.assertEqual(np.linalg.matrix_rank(np.vstack([OBSERVATION_MATRIX,HYDRAULIC_ROW])),3)

    def test_zero_and_large_oxygen_capacity(self):
        inlet=oxygen_content(40.)
        self.assertEqual(end_capillary_content(5.,0.,100.,inlet),inlet)
        result=end_capillary_content(5.,150.,100.,inlet)
        self.assertAlmostEqual(result,oxygen_content(100.),places=5)

    def test_exchange_convergence(self):
        inlet=oxygen_content(35.)
        loose=end_capillary_content(5.,8.,90.,inlet,rtol=1e-7)
        tight=end_capillary_content(5.,8.,90.,inlet,rtol=1e-10)
        self.assertLess(abs(loose-tight),1e-5)

    def test_capacity_flow_matching_is_null(self):
        inlet=oxygen_content(35.)
        uniform=end_capillary_content(5.,10.,100.,inlet)
        out=[end_capillary_content(q,d,100.,inlet) for q,d in zip([1.,4.],[2.,8.])]
        self.assertAlmostEqual(mix_beds([1.,4.],out,inlet)['content'],uniform,places=7)

    def test_heterogeneous_capacity_flow_reduces_oxygen_uptake(self):
        inlet=oxygen_content(35.)
        uniform=end_capillary_content(5.,10.,100.,inlet)
        out=[end_capillary_content(q,d,100.,inlet) for q,d in zip([1.,4.],[5.,5.])]
        self.assertLess(mix_beds([1.,4.],out,inlet)['content'],uniform)

    def test_bypass_mixing_mass_conservation(self):
        inlet=oxygen_content(35.)
        out=[end_capillary_content(q,d,100.,inlet) for q,d in zip([1.,3.],[5.,5.])]
        mixed=mix_beds([1.,3.],out,inlet,bypass_flow=1.)
        self.assertLess(abs(mixed['balance_residual']),1e-10)
        self.assertEqual(mixed['total_flow'],5.)

    def test_invalid_pressure_regime_is_explicit(self):
        with self.assertRaises(ValueError):solve_bed(self.coupling(),20.,8.,1.,.2,external_pressure=10.)
        with self.assertRaises(ValueError):self.coupling().equilibrium(0.)

    def test_fit_does_not_consume_holdout_or_oxygen(self):
        rows=load_groups()
        train=[r for r in rows if r['study']=='Fermoyle2021' and r['group']=='control' and r['stage'] in ('rest','20W')]
        original=fit_model(train,'flow')
        modified=[dict(r,pao2=1.,sao2=1.) for r in train]
        self.assertEqual(original,fit_model(modified,'flow'))
        for r in train:
            result=predict(original,r)
            for v in ('vc','dm','dlco'):self.assertAlmostEqual(result[v],r[v],places=10)


if __name__=='__main__':unittest.main()
