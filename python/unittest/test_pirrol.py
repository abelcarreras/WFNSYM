from wfnsympy import WfnSympy
from wfnsympy.file_io import get_data_from_file_fchk
import numpy as np
import unittest


class TestWfnsympy(unittest.TestCase):

    def setUp(self):
        data = get_data_from_file_fchk('pirrol.in.fchk')
        self.pirrol = WfnSympy(coordinates=data['coordinates'],
                               symbols=data['symbols'],
                               basis=data['basis'],
                               center=[0., 0., 0.],
                               axis=[1., 0., 0.],
                               # axis2=[0., 0., 1.],
                               alpha_mo_coeff=data['mo_coefficients']['alpha'][:18],
                               alpha_occupancy=data['alpha_occupancy'],
                               beta_occupancy=data['beta_occupancy'],
                               group='C6v')

    def test_symlab(self):
        symlab_test = ['E', '2C6', '2C3', 'C2', 's_v1', 's_d1', 's_v2', 's_d2', 's_v3', 's_d3']
        self.assertEqual(self.pirrol.SymLab, symlab_test)

    def test_csm_coef(self):
        csm_coef_test = [5.331142e-07, 9.999894e+01, 1.000000e+02, 1.000000e+02,
                         1.000000e+02, 1.000000e+02, 9.999894e+01, 5.331136e-07,
                         9.999894e+01, 1.000000e+02]

        np.testing.assert_allclose(csm_coef_test, self.pirrol.csm_coef, rtol=1e-6)

    def test_mo_SOEVs_a(self):
        mo_soevs_a_test = [[1.00000000e+00, -1.32522730e-03, -2.85032614e-05,
                            3.06252669e-05, 3.06252669e-05, -2.85032614e-05,
                            -1.32522730e-03, 1.00000000e+00, -1.32522730e-03,
                            -2.85032614e-05],
                           [1.00000000e+00, 1.58035432e-03, -3.10152695e-02,
                            -2.47373332e-02, 2.47373332e-02, 3.10152695e-02,
                            -1.58035432e-03, -1.00000000e+00, -1.58035432e-03,
                            3.10152695e-02],
                           [1.00000000e+00, 8.01242421e-04, 3.03781569e-02,
                            2.42405846e-02, 2.42405846e-02, 3.03781569e-02,
                            8.01242421e-04, 1.00000000e+00, 8.01242421e-04,
                            3.03781569e-02],
                           [1.00000000e+00, 2.18367913e-01, 3.02127679e-04,
                            -1.91732435e-04, -1.91732435e-04, 3.02127680e-04,
                            2.18367913e-01, 1.00000000e+00, 2.18367913e-01,
                            3.02127679e-04],
                           [1.00000000e+00, -2.24480923e-01, -1.21956191e-04,
                            9.67591395e-06, -9.67591395e-06, 1.21956191e-04,
                            2.24480923e-01, -1.00000000e+00, 2.24480923e-01,
                            1.21956191e-04],
                           [1.00000000e+00, 7.97590267e-01, 6.16140384e-01,
                            5.62311755e-01, 5.62311755e-01, 6.16140384e-01,
                            7.97590267e-01, 1.00000000e+00, 7.97590267e-01,
                            6.16140384e-01],
                           [1.00000000e+00, 5.06929014e-01, -1.34275096e-01,
                            -3.75219254e-01, -3.75219254e-01, -1.34275096e-01,
                            5.06929014e-01, 1.00000000e+00, 5.06929014e-01,
                            -1.34275096e-01],
                           [1.00000000e+00, 4.09168542e-01, -4.50809663e-01,
                            -8.28810626e-01, 8.28810626e-01, 4.50809663e-01,
                            -4.09168542e-01, -1.00000000e+00, -4.09168542e-01,
                            4.50809663e-01],
                           [9.99999998e-01, -1.45640081e-01, -7.07265257e-02,
                            4.98181993e-01, 4.98181993e-01, -7.07265257e-02,
                            -1.45640081e-01, 9.99999998e-01, -1.45640081e-01,
                            -7.07265257e-02],
                           [1.00000000e+00, -5.14408504e-01, -3.22442719e-01,
                            6.92978124e-01, -6.92978124e-01, 3.22442719e-01,
                            5.14408504e-01, -1.00000000e+00, 5.14408504e-01,
                            3.22442719e-01],
                           [9.99999998e-01, 5.68141874e-01, 4.92191987e-01,
                            6.78888006e-01, 6.78888006e-01, 4.92191987e-01,
                            5.68141874e-01, 9.99999998e-01, 5.68141874e-01,
                            4.92191987e-01],
                           [1.00000000e+00, 1.44894779e-01, -2.23530452e-01,
                            -3.63772946e-01, -3.63772946e-01, -2.23530452e-01,
                            1.44894779e-01, 1.00000000e+00, 1.44894779e-01,
                            -2.23530452e-01],
                           [9.99999999e-01, -2.09577512e-01, -5.34982970e-02,
                            -3.85232598e-01, 3.85232598e-01, 5.34982970e-02,
                            2.09577512e-01, -9.99999999e-01, 2.09577512e-01,
                            5.34982970e-02],
                           [9.99999999e-01, 8.94091923e-01, 7.88572690e-01,
                            7.46740593e-01, 7.46740593e-01, 7.88572690e-01,
                            8.94091923e-01, 9.99999999e-01, 8.94091923e-01,
                            7.88572690e-01],
                           [1.00000000e+00, -3.43439533e-01, 6.92332621e-02,
                            -3.67952436e-01, 3.67952436e-01, -6.92332621e-02,
                            3.43439533e-01, -1.00000000e+00, 3.43439533e-01,
                            -6.92332621e-02],
                           [9.99999998e-01, -6.52313367e-01, 3.71284414e-01,
                            -3.90435661e-01, -3.90435661e-01, 3.71284414e-01,
                            -6.52313367e-01, 9.99999998e-01, -6.52313367e-01,
                            3.71284414e-01],
                           [9.99999999e-01, 4.40770928e-01, -3.00193848e-01,
                            -5.75538171e-01, -5.75538171e-01, -3.00193848e-01,
                            4.40770928e-01, 9.99999999e-01, 4.40770928e-01,
                            -3.00193848e-01],
                           [1.00000000e+00, 4.37251180e-01, -4.82843548e-01,
                            -8.48062831e-01, 8.48062831e-01, 4.82843548e-01,
                            -4.37251180e-01, -1.00000000e+00, -4.37251180e-01,
                            4.82843548e-01]]

        self.pirrol.print_overlap_mo_alpha()
        np.testing.assert_allclose(mo_soevs_a_test, self.pirrol.mo_SOEVs_a, rtol=1e-6)

    def test_wf_SOEVs_a(self):
        wf_soevs_a_test = [9.99999997e-01, 3.25103201e-03, 8.71324571e-08, -2.99110905e-08,
                           2.99110905e-08, -8.71324571e-08, -3.25103201e-03, -9.99999997e-01,
                           -3.25103201e-03, -8.71324571e-08]
        np.testing.assert_allclose(wf_soevs_a_test, self.pirrol.wf_SOEVs_a, rtol=1e-6)

    def test_grim_coef(self):
        grim_coef_test = [1.480872e-08, 6.382904e+01, 7.534673e+01, 5.909258e+01,
                          5.909258e+01, 7.534673e+01, 6.382904e+01, 1.480871e-08,
                          6.382904e+01, 7.534673e+01]
        np.testing.assert_allclose(grim_coef_test, self.pirrol.grim_coef, rtol=1e-6)
