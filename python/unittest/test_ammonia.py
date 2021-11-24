from wfnsympy import WfnSympy
from wfnsympy.file_io import get_data_from_file_fchk
import numpy as np
import unittest


class TestQsympy(unittest.TestCase):

    def setUp(self):
        self.data = get_data_from_file_fchk('ammonia.fchk')
        self.structure = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis=[0.55944897, 0.69527034, -0.4512383],
                                  axis2=[0., -0.4512383, -0.69527034],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='C3v')

    def test_csm_dens(self):
        np.testing.assert_allclose(16.7899, self.structure.csm_dens, rtol=1e-03)

    def test_csm_dens_coef(self):
        csm_coef_test = [1.000000, 0.999999, 0.997533, 0.997535, 0.997535]

        print(self.structure.csm_dens_coef)
        np.testing.assert_allclose(csm_coef_test, self.structure.csm_dens_coef, rtol=1e-06)

    def test_dens_occupancy(self):
        self.structure = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis=[0.55944897, 0.69527034, -0.4512383],
                                  axis2=[0., -0.4512383, -0.69527034],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='C3v',
                                  alpha_occupancy=[0, 1, 1, 1, 1])

        np.testing.assert_allclose(20.4222478993, self.structure.csm_dens, rtol=1e-03)

    def test_precision(self):
        self.structure = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis=[0.55944897, 0.69527034, -0.4512383],
                                  axis2=[0., -0.4512383, -0.69527034],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='C3v',
                                  tolerance=1e-04,
                                  alpha_occupancy=self.data['alpha_occupancy'],
                                  beta_occupancy=self.data['beta_occupancy'])

        np.testing.assert_allclose(16.7900580, self.structure.csm_dens, rtol=1e-03)