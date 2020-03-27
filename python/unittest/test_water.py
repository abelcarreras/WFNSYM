from wfnsympy import WfnSympy
import numpy as np
import unittest


class TestWfnsympy(unittest.TestCase):

    def setUp(self):

        basis = {'name': 'STO-3G',
                 'primitive_type': 'gaussian',
                 'atoms': [{'symbol': 'O',
                            'shells': [{'shell_type': 's',
                                        'p_exponents': [130.70932, 23.808861, 6.4436083],
                                        'con_coefficients': [0.154328969, 0.535328136, 0.444634536],
                                        'p_con_coefficients': [0.0, 0.0, 0.0]},
                                       {'shell_type': 'sp',
                                        'p_exponents': [5.0331513, 1.1695961, 0.380389],
                                        'con_coefficients': [-0.0999672287, 0.399512825, 0.700115461],
                                        'p_con_coefficients': [0.155916268, 0.607683714, 0.391957386]}]},
                           {'symbol': 'H',
                            'shells': [{'shell_type': 's',
                                        'p_exponents': [3.42525091, 0.62391373, 0.1688554],
                                        'con_coefficients': [0.154328971, 0.535328142, 0.444634542],
                                        'p_con_coefficients': [0.0, 0.0, 0.0]}]},
                           {'symbol': 'H',
                            'shells': [{'shell_type': 's',
                                        'p_exponents': [3.42525091, 0.62391373, 0.1688554],
                                        'con_coefficients': [0.154328971, 0.535328142, 0.444634542],
                                        'p_con_coefficients': [0.0, 0.0, 0.0]}]}]}

        mo_coefficients = [[ 0.994216442, 0.025846814, 0.000000000, 0.000000000,-0.004164076,-0.005583712,-0.005583712],
                           [ 0.233766661,-0.844456594, 0.000000000, 0.000000000, 0.122829781,-0.155593214,-0.155593214],
                           [ 0.000000000, 0.000000000, 0.612692349, 0.000000000, 0.000000000,-0.449221684, 0.449221684],
                           [-0.104033343, 0.538153649, 0.000000000, 0.000000000, 0.755880259,-0.295107107,-0.295107107],
                           [ 0.000000000, 0.000000000, 0.000000000,-1.000000000, 0.000000000, 0.000000000, 0.000000000],
                           [-0.125818566, 0.820120983, 0.000000000, 0.000000000,-0.763538862,-0.769155124,-0.769155124],
                           [ 0.000000000, 0.000000000, 0.959800163, 0.000000000, 0.000000000, 0.814629717,-0.814629717]]

        self.wf_results = WfnSympy(coordinates=[[ 0.0000000000, 0.0000000000, -0.0428008531],
                                                [-0.7581074140, 0.0000000000, -0.6785995734],
                                                [ 0.7581074140, 0.0000000000, -0.6785995734]],
                                   symbols=['O', 'H', 'H'],
                                   basis=basis,
                                   alpha_mo_coeff=mo_coefficients,
                                   charge=0,
                                   multiplicity=1,
                                   group='C2v')

    def test_csm_coef(self):
        csm_coef_test = [0, 0, 0, 0]

        np.testing.assert_allclose(csm_coef_test, self.wf_results.csm_coef, atol=1e-5)

    def test_mo_SOEVs_a(self):
        mo_soevs_a_test = [[1, 1, 1, 1],
                           [1, 1, 1, 1],
                           [1,-1,-1, 1],
                           [1, 1, 1, 1],
                           [1,-1, 1,-1]]

        np.testing.assert_allclose(mo_soevs_a_test, self.wf_results.mo_SOEVs_a, atol=1e-5)

    def test_wf_SOEVs_a(self):
        wf_soevs_a_test = [ 1,  1, -1, -1]
        np.testing.assert_allclose(wf_soevs_a_test, self.wf_results.wf_SOEVs_a, atol=1e-5)

    def test_grim_coef(self):
        grim_coef_test = [0, 0, 0, 0]
        np.testing.assert_allclose(grim_coef_test, self.wf_results.grim_coef, atol=1e-5)

