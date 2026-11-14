import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from validate import *

class ValidationTests(unittest.TestCase):
    def test_zero_and_equilibrium(self):
        r=read_data()[0]
        for law in LAWS[:-1]:
            self.assertAlmostEqual(predict(law,0,r),r['Pv_mmHg'],places=7)
        self.assertAlmostEqual(predict('equilibrium_admixture',0,r),r['PA_mmHg'],places=7)
        self.assertAlmostEqual(predict('equilibrium_admixture',1,r),r['Pv_mmHg'],places=7)
    def test_no_test_target_leakage(self):
        r=read_data()[1];base=read_data()[0]
        for law in LAWS:
            p=fit_one(law,base)
            self.assertEqual(predict(law,p,r),predict(law,p,dict(r,Pa_mmHg=999)))
    def test_common_endpoint(self):
        for r in read_data():
            if r['baseline']:
                for law in LAWS:self.assertAlmostEqual(predict(law,fit_one(law,r),r),r['Pa_mmHg'],places=7)
    def test_flow_index_invariance(self):
        r=read_data()[0]
        for law in LAWS:
            p=fit_one(law,r);p2=p if law=='equilibrium_admixture' else p/1.9
            self.assertAlmostEqual(predict(law,p,r),predict(law,p2,dict(r,flow=r['flow']/1.9)),places=7)
    def test_hyperoxic_alveolar_equation(self):
        r=read_data()[2]
        self.assertAlmostEqual(r['PA_mmHg'],r['PB_mmHg']-r['PH2O_mmHg']-r['PaCO2_mmHg'])
    def test_closed_balance_and_admixture_formula(self):
        from couple_network import hydraulic,closed,pars,base
        h=hydraulic(8);o=closed(h,'equilibrium_admixture')
        s=pars['equilibrium_admixture'];q=h['CO_L_min'];ce=content(base['PA_mmHg'],14)
        self.assertAlmostEqual(o['Ca_ml_dl'],ce-s/(1-s)*250/(10*q),places=8)
        self.assertLess(abs(o['mass_balance_error']),1e-7)
    def test_infeasible_demand_not_clipped(self):
        from couple_network import hydraulic,closed
        self.assertFalse(closed(hydraulic(16),'constant_D',hb=10,vo2=10000)['feasible'])

if __name__=='__main__':unittest.main()
