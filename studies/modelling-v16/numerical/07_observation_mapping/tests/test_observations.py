import sys
import unittest
import math
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from observations import *


class ObservationChecks(unittest.TestCase):
    def test_two_gas_inverse_and_forward_roundtrip(self):
        for dm,vc in [(20.,80.),(40.,100.),(80.,60.)]:
            co=gas_conductance(dm,vc,.5);no=gas_conductance(dm,vc,4.5,1.97)
            out=invert_two_gases(co,no,.5)
            self.assertAlmostEqual(out['dm'],dm,places=10)
            self.assertAlmostEqual(out['vc'],vc,places=10)

    def test_infinite_NO_limit(self):
        co=gas_conductance(40.,100.,.5)
        out=invert_two_gases(co,1.97*40,.5,theta_no=math.inf)
        self.assertAlmostEqual(out['dm'],40.,places=10)
        self.assertAlmostEqual(out['vc'],100.,places=10)

    def test_same_observations_can_have_different_oxygen_estimates(self):
        first=invert_two_gases(16.4,48.,.3)
        second=invert_two_gases(16.4,48.,1.)
        for out,theta in [(first,.3),(second,1.)]:
            self.assertAlmostEqual(gas_conductance(out['dm'],out['vc'],theta),16.4,places=10)
            self.assertAlmostEqual(gas_conductance(out['dm'],out['vc'],4.5,1.97),48.,places=10)
        self.assertGreater(abs(oxygen_estimate(first['dm'],first['vc'])['effective']-
                               oxygen_estimate(second['dm'],second['vc'])['effective']),1.)

    def test_log_sensitivities_match_finite_differences(self):
        dm,vc=40.,100.;jac=log_gas_jacobian(dm,vc,.5)
        def values(x):
            return np.log([gas_conductance(x[0],x[1],.5),gas_conductance(x[0],x[1],4.5,1.97)])
        for j in range(2):
            plus=np.array([dm,vc]);minus=plus.copy();plus[j]*=math.exp(1e-5);minus[j]*=math.exp(-1e-5)
            np.testing.assert_allclose((values(plus)-values(minus))/(2e-5),jac[:,j],rtol=1e-8)

    def test_NO_error_can_induce_opposing_Dm_Vc_changes(self):
        initial=invert_two_gases(16.4,48.,.5)
        raised=invert_two_gases(16.4,48.*1.01,.5)
        self.assertGreater(raised['dm'],initial['dm'])
        self.assertLess(raised['vc'],initial['vc'])

    def test_oxygen_series_bounds(self):
        e=oxygen_estimate(40.,100.)
        self.assertLess(e['effective'],e['membrane']);self.assertLess(e['effective'],e['blood'])
        self.assertAlmostEqual(e['membrane_resistance_fraction']+e['blood_resistance_fraction'],1.)

    def test_oxygen_elasticities(self):
        e=oxygen_estimate(40.,100.);step=1e-5
        for key,dm,vc in [('membrane_resistance_fraction',40*math.exp(step),100.),
                          ('blood_resistance_fraction',40.,100*math.exp(step))]:
            derivative=math.log(oxygen_estimate(dm,vc)['effective']/e['effective'])/step
            self.assertAlmostEqual(derivative,e[key],places=5)

    def test_mechanism_null_leaves_function_unchanged(self):
        base=np.array([.9,1.1,1.3]);other=base*np.exp(.1*GEOMETRY_NULL)
        def outputs(s):
            n,d,k=s;return 85*n*d*d,30*n*d*k,n*d**4
        v1,m1,g1=outputs(base);v2,m2,g2=outputs(other)
        self.assertAlmostEqual(v1,v2);self.assertAlmostEqual(m1,m2)
        self.assertAlmostEqual(oxygen_estimate(m1,v1)['effective'],oxygen_estimate(m2,v2)['effective'])
        self.assertNotAlmostEqual(g1,g2)

    def test_additional_oxygen_and_transit_do_not_add_geometry_rank(self):
        matrices=geometry_observation_rows(.8)
        for name in ['base','plus_oxygen','plus_transit']:self.assertEqual(np.linalg.matrix_rank(matrices[name]),2)
        for name in ['plus_area','plus_hydraulics']:self.assertEqual(np.linalg.matrix_rank(matrices[name]),3)
        self.assertAlmostEqual(matrices['oxygen']@GEOMETRY_NULL,0.)

    def test_gain_threshold_predicts_sign(self):
        rm=138/90;rb=85/98;critical=gain_threshold(rm,rb)
        for z in [.5*critical,2*critical]:
            initial=oxygen_estimate(1/1.23,z,theta_o2=1)['effective']
            final=oxygen_estimate(rm/1.23,z*rb,theta_o2=1)['effective']
            self.assertEqual(final>initial,z>critical)

    def test_multiple_FiO2_line_unchanged_under_reparameterization(self):
        slope=.0058/100;intercept=1/40+.73/100
        for a,b in [(.0058,.73),(.0062,1.16)]:
            out=multifio2_reparameterize(intercept,slope,a,b)
            for pa in [80.,150.,300.]:
                resistance=1/out['dm']+(a*pa+b)/out['vc']
                self.assertAlmostEqual(resistance,intercept+slope*pa,places=12)

    def test_unit_conversion_to_SI(self):
        # 1 mmol gas STPD = 22.414 mL; 1 kPa = 7.50062 mmHg.
        factor=22.414/7.50062
        customary=oxygen_estimate(40.,100.,theta_o2=.6*factor)['effective']
        si=oxygen_estimate(40./factor,100.,theta_o2=.6)['effective']
        self.assertAlmostEqual(customary/factor,si,places=12)
        self.assertAlmostEqual(.6*factor,1.8,delta=.01)

    def test_physical_domain_rejections(self):
        with self.assertRaises(ValueError):invert_two_gases(20.,30.,.5)
        with self.assertRaises(ValueError):invert_two_gases(20.,60.,4.5/1.97)
        with self.assertRaises(ValueError):oxygen_estimate(30.,0.)

    def test_upstream_pressure_is_not_unique_capillary_driver(self):
        jac=pressure_driver_jacobian(5.,1.,.2,.3)
        upstream=jac@np.array([0.,0.,1.,0.,0.])
        redistribution=jac@np.array([0.,0.,-1.,0.,1.])
        np.testing.assert_allclose(upstream,[5.,0.])
        np.testing.assert_allclose(redistribution,[0.,5.])


if __name__=='__main__':unittest.main()
