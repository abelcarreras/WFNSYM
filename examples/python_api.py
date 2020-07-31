from wfnsympy import WfnSympy
from wfnsympy.file_io import get_data_from_file_fchk
import numpy as np

# data = get_data_from_file_fchk('pirrol.in.fchk')
data = get_data_from_file_fchk('methane.fchk')

structure = WfnSympy(coordinates=data['coordinates'],
                     symbols=data['symbols'],
                     basis=data['basis'],
                     alpha_mo_coeff=data['mo_coefficients']['alpha'],
                     # center=[0., 0., 0.],
                     # axis=[0., 0., 1.],
                     # axis2=[0, 1, 0],
                     alpha_occupancy=data['alpha_occupancy'],
                     beta_occupancy=data['beta_occupancy'],
                     # charge=0,
                     # multiplicity=1,
                     group='Td')


structure.print_CSM()
structure.print_overlap_mo_alpha()
structure.print_dens_CSM()
print('axis: ', structure.axis, np.linalg.norm(structure.axis))
print('axis2: ', structure.axis2, np.linalg.norm(structure.axis2))
print('center: ', structure.center)
