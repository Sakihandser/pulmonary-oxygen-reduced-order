import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from primary_model import solve,delivery_derivatives
from validate import content
from scipy.optimize import brentq

class PrimaryModelTests(unittest.TestCase):
    def test_independent_balance_root(self):
        for q in [2.5,4.,6.]:
            for s in [0.,.15,.3]:
                o=solve(q,14,100,250,s);ce=content(100,14)
                cv=brentq(lambda v:10*q*((1-s)*ce+s*v-v)-250,0,ce)
                self.assertAlmostEqual(o['Cv_ml_dl'],cv,places=10)
                self.assertLess(abs(o['mass_balance_error']),1e-9)
    def test_exact_feasibility_boundary(self):
        qmin=solve(4,14,100,250,.15)['Q_min_L_min']
        self.assertTrue(solve(qmin,14,100,250,.15)['feasible'])
        self.assertAlmostEqual(solve(qmin,14,100,250,.15)['Cv_ml_dl'],0,places=10)
        o=solve(qmin*.999,14,100,250,.15)
        self.assertFalse(o['feasible']);self.assertEqual(o['Sa_percent'],'')
    def test_zero_mixing_and_zero_demand_limits(self):
        a=solve(4,14,100,250,0);b=solve(4,14,100,500,0)
        self.assertAlmostEqual(a['Pa_mmHg'],100,places=8)
        self.assertAlmostEqual(a['DO2_ml_min'],b['DO2_ml_min'],places=8)
        z=solve(4,14,100,0,.8)
        self.assertAlmostEqual(z['Ca_ml_dl'],z['Cv_ml_dl'],places=10)
    def test_invalid_parameters(self):
        for name,value in [('q',0),('hb',-1),('pa',0),('demand',-1),('admixture',1),('admixture',-.1),('q',math.nan)]:
            args=dict(q=4,hb=14,pa=100,demand=250,admixture=.15);args[name]=value
            with self.assertRaises(ValueError):solve(**args)
    def test_analytic_derivatives_against_finite_difference(self):
        args=dict(q=4,hb=14,pa=100,demand=250,admixture=.15)
        d=delivery_derivatives(**args)
        for label,key in [('Q','q'),('Hb','hb'),('PA','pa'),('M','demand'),('s_eff','admixture')]:
            eps=1e-5
            plus=solve(**dict(args,**{key:args[key]+eps}))['DO2_ml_min']
            minus=solve(**dict(args,**{key:args[key]-eps}))['DO2_ml_min']
            self.assertAlmostEqual((plus-minus)/(2*eps),d[label],delta=1e-5)
    def test_pressure_redistribution_cannot_change_primary_at_fixed_flow(self):
        from couple_network import hydraulic,closed
        a=hydraulic(3,1);b=hydraulic(12,25);b['CO_L_min']=a['CO_L_min']
        self.assertEqual(closed(a,'equilibrium_admixture'),closed(b,'equilibrium_admixture'))

if __name__=='__main__':unittest.main()
