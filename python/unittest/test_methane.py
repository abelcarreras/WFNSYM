from wfnsympy import WfnSympy
from wfnsympy.file_io import get_data_from_file_fchk
import numpy as np
import unittest


class TestWfnsympy(unittest.TestCase):

    def setUp(self):
        self.data = get_data_from_file_fchk('methane.fchk')

    def test_degenerate_1(self):

        configurations = [
            [[1, 1,  1, 0, 0,  0, 0, 0,  0], [1, 1,  0, 1, 0,  0, 0, 0,  0]],
            [[1, 1,  1, 0, 0,  0, 0, 0,  0], [1, 1,  0, 0, 1,  0, 0, 0,  0]],
            [[1, 1,  0, 1, 0,  0, 0, 0,  0], [1, 1,  1, 0, 0,  0, 0, 0,  0]],
            [[1, 1,  0, 1, 0,  0, 0, 0,  0], [1, 1,  0, 0, 1,  0, 0, 0,  0]],
            [[1, 1,  0, 0, 1,  0, 0, 0,  0], [1, 1,  1, 0, 0,  0, 0, 0,  0]],
            [[1, 1,  0, 0, 1,  0, 0, 0,  0], [1, 1,  0, 1, 0,  0, 0, 0,  0]],
        ]

        for conf_alpha, conf_beta in configurations:

            print(conf_alpha, conf_beta)

            wf_results = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis= [-0.5634811520306573,  0.5937396995245904, -0.5744233286650642],
                                  axis2=[-0.5910397590456458, -0.5756503475754583,  0.5650652002764983],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='Td',
                                  alpha_occupancy=conf_alpha,
                                  beta_occupancy=conf_beta,
                                  )

            wf_results.print_alpha_mo_IRD()
            wf_results.print_wf_mo_IRD()

            wf_alpha = np.round(wf_results.wf_IRd_a, decimals=2)
            wf_beta = np.round(wf_results.wf_IRd_b, decimals=2)
            wf_total = np.round(wf_results.wf_IRd, decimals=2)

            print(wf_alpha, wf_beta, wf_total)

            self.assertTupleEqual(tuple(wf_alpha), (0., 0., 0., 0.,  1.))
            self.assertTupleEqual(tuple(wf_beta),  (0., 0., 0., 0.,  1.))
            self.assertTupleEqual(tuple(wf_total), (0., 0., 0., 0.5, 0.5))


    def test_degenerate_2(self):

        configurations = [
            [[1, 1,  1, 0, 0,  0, 0, 0,  0], [1, 1,  1, 0, 0,  0, 0, 0,  0]],
            [[1, 1,  0, 1, 0,  0, 0, 0,  0], [1, 1,  0, 1, 0,  0, 0, 0,  0]],
            [[1, 1,  0, 0, 1,  0, 0, 0,  0], [1, 1,  0, 0, 1,  0, 0, 0,  0]],
        ]

        for conf_alpha, conf_beta in configurations:

            print(conf_alpha, conf_beta)

            wf_results = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis= [-0.5634811520306573,  0.5937396995245904, -0.5744233286650642],
                                  axis2=[-0.5910397590456458, -0.5756503475754583,  0.5650652002764983],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='Td',
                                  alpha_occupancy=conf_alpha,
                                  beta_occupancy=conf_beta,
                                  )

            wf_results.print_alpha_mo_IRD()
            wf_results.print_wf_mo_IRD()

            wf_alpha = np.round(wf_results.wf_IRd_a, decimals=2)
            wf_beta = np.round(wf_results.wf_IRd_b, decimals=2)
            wf_total = np.round(wf_results.wf_IRd, decimals=2)

            print(wf_alpha, wf_beta, wf_total)

            self.assertTupleEqual(tuple(wf_alpha), (0., 0., 0., 0.,  1.))
            self.assertTupleEqual(tuple(wf_beta),  (0., 0., 0., 0.,  1.))
            self.assertTupleEqual(tuple(wf_total), (0.33, 0., 0.67, 0., 0.))

    def test_degenerate_3(self):

        configurations = [
            [[1, 1,  1, 1, 1,  1, 0, 0,  0], [1, 1,  1, 1, 1,  0, 1, 0,  0]],
            [[1, 1,  1, 1, 1,  1, 0, 0,  0], [1, 1,  1, 1, 1,  0, 0, 1,  0]],
            [[1, 1,  1, 1, 1,  0, 1, 0,  0], [1, 1,  1, 1, 1,  1, 0, 0,  0]],
            [[1, 1,  1, 1, 1,  0, 1, 0,  0], [1, 1,  1, 1, 1,  0, 0, 1,  0]],
            [[1, 1,  1, 1, 1,  0, 0, 1,  0], [1, 1,  1, 1, 1,  1, 0, 0,  0]],
            [[1, 1,  1, 1, 1,  0, 0, 1,  0], [1, 1,  1, 1, 1,  0, 1, 0,  0]],
        ]

        for conf_alpha, conf_beta in configurations:

            print(conf_alpha, conf_beta)

            wf_results = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis= [-0.5634811520306573,  0.5937396995245904, -0.5744233286650642],
                                  axis2=[-0.5910397590456458, -0.5756503475754583,  0.5650652002764983],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='Td',
                                  alpha_occupancy=conf_alpha,
                                  beta_occupancy=conf_beta,
                                  )

            wf_results.print_alpha_mo_IRD()
            wf_results.print_wf_mo_IRD()

            wf_alpha = np.round(wf_results.wf_IRd_a, decimals=2)
            wf_beta = np.round(wf_results.wf_IRd_b, decimals=2)
            wf_total = np.round(wf_results.wf_IRd, decimals=2)

            print(wf_alpha, wf_beta, wf_total)

            self.assertTupleEqual(tuple(wf_alpha), (0., 0., 0., 1.,  0.))
            self.assertTupleEqual(tuple(wf_beta),  (0., 0., 0., 1.,  0.))
            self.assertTupleEqual(tuple(wf_total), (0., 0., 0., 0.5, 0.5))

    def test_degenerate_4(self):

        configurations = [
            [[1, 1,  1, 1, 1,  1, 0, 0,  0], [1, 1,  1, 1, 1,  1, 0, 0,  0]],
            [[1, 1,  1, 1, 1,  0, 1, 0,  0], [1, 1,  1, 1, 1,  0, 1, 0,  0]],
            [[1, 1,  1, 1, 1,  0, 0, 1,  0], [1, 1,  1, 1, 1,  0, 0, 1,  0]],
        ]

        for conf_alpha, conf_beta in configurations:

            print(conf_alpha, conf_beta)

            wf_results = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis= [-0.5634811520306573,  0.5937396995245904, -0.5744233286650642],
                                  axis2=[-0.5910397590456458, -0.5756503475754583,  0.5650652002764983],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='Td',
                                  alpha_occupancy=conf_alpha,
                                  beta_occupancy=conf_beta,
                                  )

            wf_results.print_alpha_mo_IRD()
            wf_results.print_wf_mo_IRD()

            wf_alpha = np.round(wf_results.wf_IRd_a, decimals=2)
            wf_beta = np.round(wf_results.wf_IRd_b, decimals=2)
            wf_total = np.round(wf_results.wf_IRd, decimals=2)

            print(wf_alpha, wf_beta, wf_total)

            self.assertTupleEqual(tuple(wf_alpha), (0., 0., 0., 1.,  0.))
            self.assertTupleEqual(tuple(wf_beta),  (0., 0., 0., 1.,  0.))
            self.assertTupleEqual(tuple(wf_total), (0.33, 0., 0.66, 0., 0.))


    def test_degenerate_5(self):

        configurations = [
            [[1, 1,  0, 0, 0,  0, 0, 0,  0], [1, 1,  0, 0, 0,  0, 0, 0,  0]],
            [[1, 0,  0, 0, 0,  0, 0, 0,  0], [1, 0,  0, 0, 0,  0, 0, 0,  0]],
            [[0, 1,  0, 0, 0,  0, 0, 0,  0], [0, 1,  0, 0, 0,  0, 0, 0,  0]],
            [[1, 0,  0, 0, 0,  0, 0, 0,  0], [0, 1,  0, 0, 0,  0, 0, 0,  0]],
            [[0, 1,  0, 0, 0,  0, 0, 0,  0], [1, 0,  0, 0, 0,  0, 0, 0,  0]],
        ]

        for conf_alpha, conf_beta in configurations:

            print(conf_alpha, conf_beta)

            wf_results = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis= [-0.5634811520306573,  0.5937396995245904, -0.5744233286650642],
                                  axis2=[-0.5910397590456458, -0.5756503475754583,  0.5650652002764983],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='Td',
                                  alpha_occupancy=conf_alpha,
                                  beta_occupancy=conf_beta,
                                  )

            wf_results.print_alpha_mo_IRD()
            wf_results.print_wf_mo_IRD()

            wf_alpha = np.round(wf_results.wf_IRd_a, decimals=2)
            wf_beta = np.round(wf_results.wf_IRd_b, decimals=2)
            wf_total = np.round(wf_results.wf_IRd, decimals=2)

            print(wf_alpha, wf_beta, wf_total)

            self.assertTupleEqual(tuple(wf_alpha), (1., 0., 0., 0.,  0.))
            self.assertTupleEqual(tuple(wf_beta),  (1., 0., 0., 0.,  0.))
            self.assertTupleEqual(tuple(wf_total), (1., 0., 0.0, 0., 0.))

    def test_degenerate_6(self):

        configurations = [
            [[1, 1,  1, 1, 1,  1, 1, 1,  1], [1, 1,  1, 1, 1,  1, 1, 1,  0]],
            [[1, 1,  1, 1, 1,  1, 1, 1,  0], [1, 0,  1, 1, 1,  1, 1, 1,  1]],
            [[1, 1,  1, 1, 1,  1, 1, 1,  1], [0, 1,  1, 1, 1,  1, 1, 1,  1]],
        ]

        for conf_alpha, conf_beta in configurations:

            print(conf_alpha, conf_beta)

            wf_results = WfnSympy(coordinates=self.data['coordinates'],
                                  symbols=self.data['symbols'],
                                  basis=self.data['basis'],
                                  axis= [-0.5634811520306573,  0.5937396995245904, -0.5744233286650642],
                                  axis2=[-0.5910397590456458, -0.5756503475754583,  0.5650652002764983],
                                  alpha_mo_coeff=self.data['mo_coefficients']['alpha'],
                                  group='Td',
                                  alpha_occupancy=conf_alpha,
                                  beta_occupancy=conf_beta,
                                  )

            wf_results.print_alpha_mo_IRD()
            wf_results.print_wf_mo_IRD()

            wf_alpha = np.round(wf_results.wf_IRd_a, decimals=1)
            wf_beta = np.round(wf_results.wf_IRd_b, decimals=1)
            wf_total = np.round(wf_results.wf_IRd, decimals=1)

            print(wf_alpha, wf_beta, wf_total)

            self.assertTupleEqual(tuple(wf_alpha), (1., 0., 0., 0., 0.))
            self.assertTupleEqual(tuple(wf_beta),  (1., 0., 0., 0., 0.))
            self.assertTupleEqual(tuple(wf_total), (1., 0., 0., 0., 0.))
