import unittest, sys, math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import analyze_heterogeneity as a

class HeterogeneousChecks(unittest.TestCase):
    def test_01_independent_venous_balance(self):
        for w,z in a.models().values():
            for r in [1.,1.3,2.]:
                for beta in [0.,.4,1.]:
                    st=a.state(w,z,r,beta)
                    zr=z**(r**(-beta));q=a.Q0*r
                    def balance(cv):
                        outlets=a.CE-(a.CE-cv)*zr
                        return 10*q*(float(w@outlets)-cv)-a.M
                    cv=brentq(balance,0,a.CE)
                    ca=float(w@(a.CE-(a.CE-cv)*zr))
                    self.assertAlmostEqual(cv,st['Cv'],11)
                    self.assertAlmostEqual(ca,st['Ca'],11)

    def test_02_derivatives_and_delivery(self):
        for w,z in a.models().values():
            for r in [1.,1.5,2.]:
                for beta in [0.,.5,1.]:
                    st=a.state(w,z,r,beta);dr=1e-5
                    for key,deriv in [('Ca','dCa_dQ'),('DO2','dDO2_dQ')]:
                        fd=(a.state(w,z,r+dr,beta)[key]-a.state(w,z,r-dr,beta)[key])/(2*dr*a.Q0)
                        self.assertAlmostEqual(fd,st[deriv],7)
                    self.assertGreater(st['dDO2_dQ'],0)

    def test_03_homogeneous_limit(self):
        for beta in [0.,.25,.5,1.]:
            for r in [1.,1.6,2.]:
                st=a.state([1],[a.SSTAR],r,beta)
                old=a.prior.state(a.Q0*r,a.M,a.CE,a.Q0,a.SSTAR,beta)
                for key in ['s','Ca','Cv','DO2','dCa_dQ','dDO2_dQ']:
                    self.assertAlmostEqual(st[key],old[key],10)

    def test_04_sharp_attainers(self):
        for lo,hi in [(0.,1.),(.001,.8),(.05,.5)]:
            wh=(a.SSTAR-lo)/(hi-lo)
            for beta in [0.,.5,1.]:
                for r in [1.,1.25,2.]:
                    b=a.bounds(a.SSTAR,r,beta,lo,hi)
                    self.assertAlmostEqual(a.state([1-wh,wh],[lo,hi],r,beta)['s'],b['s_lower'],13)
                    self.assertAlmostEqual(a.state([1],[a.SSTAR],r,beta)['s'],b['s_upper'],13)
                    self.assertAlmostEqual(float(np.array([1-wh,wh])@a.g([lo,hi])),b['h_lower'],13)
        errors=[]
        for eps in [1e-3,1e-6,1e-12]:
            hi=1-eps;w=(a.SSTAR-eps)/(hi-eps)
            errors.append(a.state([1-w,w],[eps,hi],2,1)['s']-a.SSTAR)
        self.assertTrue(errors[2]<errors[1]<errors[0])

    def test_05_random_mixtures_within_bounds(self):
        rng=np.random.default_rng(42128)
        for _ in range(200):
            w=rng.dirichlet(np.ones(7));z=rng.uniform(.001,.9,7)
            mean=float(w@z);beta=float(rng.uniform());r=float(rng.uniform(1,2))
            b=a.bounds(mean,r,beta,.001,.9);st=a.state(w,z,r,beta)
            self.assertLessEqual(b['s_lower']-1e-13,st['s'])
            self.assertLessEqual(st['s'],b['s_upper']+1e-13)
            h=float(w@a.g(z))
            self.assertLessEqual(b['h_lower']-1e-13,h)
            self.assertLessEqual(h,b['h_upper']+1e-13)
            if st['feasible']:self.assertGreater(st['dDO2_dQ'],0)

    def test_06_matched_capacity_opposite_signs(self):
        models=a.models();wa,za=models['finite_A'];wb,zb=models['matched_B']
        self.assertAlmostEqual(float(wa@za),float(wb@zb),14)
        self.assertAlmostEqual(float(wa@(-np.log(za))),float(wb@(-np.log(zb))),13)
        self.assertGreater(a.state(wa,za,1,1)['dCa_dQ'],0)
        self.assertLess(a.state(wb,zb,1,1)['dCa_dQ'],0)
        for r in [1,1.5,2]:
            # beta=1: each component capacity w_i Q0 x_i is unchanged.
            self.assertAlmostEqual(float(wa@(r*(-np.log(za)/r))),float(wb@(r*(-np.log(zb)/r))),13)

    def test_07_recovery_classification(self):
        target=a.content_at(a.S0,1)[0]
        self.assertEqual(a.bounds(.2,2,.25,target=target)['status'],'guaranteed')
        self.assertEqual(a.bounds(.2,2,1,target=target)['status'],'unresolved')
        self.assertEqual(a.bounds(.2,2,1,.05,.5,target=target)['status'],'impossible')
        w,z=a.models()['finite_A'];r=a.first_recovery(w,z,1,target,'Ca')
        self.assertTrue(1<r<2);self.assertAlmostEqual(a.state(w,z,r,1)['Ca'],target,11)
        w,z=a.models()['matched_B']
        self.assertTrue(all(a.state(w,z,float(r),1)['Ca']<target for r in np.linspace(1,2,201)))

    def test_08_feasible_intersection(self):
        b=a.bounds(.2,1,1,m=1000)
        self.assertEqual(b['status'],'no_feasible_model')
        b=a.bounds(.2,2,1,m=1300)
        self.assertLess(b['feasible_s_upper'],b['s_upper'])
        ca,cv,_=a.content_at(b['feasible_s_upper'],2,m=1300)
        self.assertAlmostEqual(cv,0,12)
        self.assertAlmostEqual(ca,b['Ca_lower'],12)

    def test_09_zero_consumption_and_domains(self):
        st=a.state([.5,.5],[.1,.3],1,1,m=0)
        self.assertEqual(st['Ca'],a.CE);self.assertEqual(st['dCa_dQ'],0)
        for w,z,r,beta in [([1],[1],1,1),([.4],[.2],1,1),([1],[-.2],1,1),([1],[.2],0,1),([1],[.2],1,2)]:
            with self.assertRaises(ValueError):a.state(w,z,r,beta)
        with self.assertRaises(ValueError):a.bounds(.2,.9,1)
        with self.assertRaises(ValueError):a.bounds(.2,2,1,.3,.8)
        b=a.bounds(.2,2,1,.2,.2);self.assertAlmostEqual(b['s_lower'],b['s_upper'])

    def test_10_effective_and_total_capacity_are_distinct(self):
        for name in ['finite_A','matched_B']:
            w,z=a.models()[name]
            h=1e-5
            k=lambda r:-a.Q0*r*math.log(a.state(w,z,r,1)['s'])
            observed=(k(1+h)-k(1-h))/(2*h*a.Q0)
            st=a.state(w,z,1,1)
            expected=-math.log(st['s'])-st['h']/st['s']
            self.assertAlmostEqual(observed,expected,8)
            self.assertGreater(observed,0)
        # The single-exponential criterion still applies to effective elasticity,
        # but not to the separately fixed sum of physical path capacities.
        for name in ['finite_A','matched_B']:
            w,z=a.models()[name];st=a.state(w,z,1,1)
            effective_beta=st['h']/(st['s']*(-math.log(st['s'])))
            critical=(1-st['s'])/(-math.log(st['s']))
            self.assertEqual(st['dCa_dQ']>0,effective_beta<critical)

if __name__=='__main__':unittest.main()
