__version__ = '0.2.3'

from wfnsympy.WFNSYMLIB import mainlib, overlap_mat
import numpy as np
import sys, io

_bohr_to_angstrom = 0.529177249

# define assignation of shell type to number and number of functions by shell
shell_type_list = {'-1': ['sp', 4],
                   '0': ['s',  1],
                   '1': ['p',  3],
                   '2': ['d',  6],
                   '3': ['f',  10],
                   '-2': ['d_', 5],  # pure
                   '-3': ['f_', 7]}  # pure


class _captured_stdout:
    def __init__(self):
        self.old_stdout = None
        self.fnull = None

    def __enter__(self):
        import sys, os, tempfile
        self.old_stdout = sys.stdout
        #self.fnull = open(os.devnull, 'w')
        self.F = tempfile.NamedTemporaryFile()
        try:
            os.dup2(self.F.fileno(), sys.stderr.fileno())
        except (AttributeError, io.UnsupportedOperation):
            pass
        return self.F

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self.old_stdout
        #self.fnull.close()
        self.F.close()


def get_valence_electrons(atomic_numbers, charge):
    valence_electrons = 0
    for number in atomic_numbers:
        if 2 >= number > 0:
            valence_electrons += np.mod(number, 2)
        if 18 >= number > 2:
            valence_electrons += np.mod(number-2, 8)
        if 54 >= number > 18:
            valence_electrons += np.mod(number-18, 18)
        if 118 >= number > 54:
            valence_electrons += np.mod(number-54, 32)
        if number > 118:
            raise Exception('Atomic number size not implemented')

    valence_electrons -= charge
    return valence_electrons


def _get_group_num_from_label(label):

    label_2 = label[0].upper()
    try:
        ngroup = int(label[1])
        label_2 += 'n'
        try:
            label_2 += label[2:].upper()
        except:
            pass
    except:
        label_2 += label[1:]
        ngroup = 0

    if label.upper() == 'S2':
        label_2 = 'CI'

    if label_2[0].upper() == 'S':
        if ngroup % 2 != 0:
            label_2 = list(label_2)
            label_2[0] = 'C'
            label_2 = ''.join(label_2)

    operations = {'CN':  [1, ngroup],
                  'CNH': [2, ngroup],
                  'CNV': [3, ngroup],
                  'CI':  [0, 2],
                  'CINF':[9, 1],
                  'DN':  [4, ngroup],
                  'DNH': [5, ngroup],
                  'DND': [6, ngroup],
                  'DINF':[9, 2],
                  'SN':  [7, ngroup],
                  'T':   [8, 1],
                  'TH':  [8, 2],
                  'TD':  [8, 3],
                  'O':   [8, 4],
                  'OH':  [8, 5],
                  'I':   [8, 6],
                  'IH':  [8, 7],
                  }
    try:
        return operations[label_2.upper()]
    except KeyError:
        raise Exception('Label not found')


def _get_operation_num_from_label(label):
    operations = {'I': 1,
                  'R': 2,
                  'C': 3,
                  'S': 4}
    ioper = operations[label[0]]
    if ioper > 2:
        irot = int(label[1])
    else:
        irot = 0
    return ioper, irot


def get_perpendicular_axis(axis):
    axis = np.array(axis)
    if np.abs(np.dot(axis, [1, 0, 0])) < np.abs(np.dot(axis, [0, 1, 0])):
        return np.cross(axis, [1, 0, 0])
    else:
        return np.cross(axis, [0, 1, 0])


def _center_of_charge(mo_coefficients_alpha, mo_coefficients_beta,
                      coordinates, basis, total_electrons, multiplicity,
                      overlap_matrix):
    """
    Returns the center of charge in Angstrom
    """

    alpha_unpaired = multiplicity//2 + 1 if (total_electrons % 2) else multiplicity//2

    alpha_electrons = total_electrons//2 + alpha_unpaired
    beta_electrons = total_electrons - alpha_electrons
    # print('electrons', alpha_electrons, beta_electrons)

    type_to_nfunc = {}
    for item in shell_type_list.items():
        type_to_nfunc['{}'.format(item[1][0])] = int(item[1][1])

    # get the basis functions corresponding to each atom (ini, fin)
    ranges_per_atom = []
    n_start = 0
    for atoms in basis['atoms']:
        n_functions = 0
        for shell in atoms['shells']:
            n_functions += type_to_nfunc[shell['shell_type']]

        ranges_per_atom.append((n_start, n_start + n_functions))
        n_start += n_functions

    # localization on fragments analysis
    number_of_atoms = len(coordinates)

    charges = []
    for atom in range(number_of_atoms):
        charge_atom = 0
        # Alpha
        for i in range(alpha_electrons):
            orb = mo_coefficients_alpha[i]
            orb_atom = np.zeros_like(orb)
            orb_atom[ranges_per_atom[atom][0]:ranges_per_atom[atom][1]] = \
                 orb[ranges_per_atom[atom][0]:ranges_per_atom[atom][1]]
            charge_atom += np.dot(orb_atom, np.dot(overlap_matrix, orb))
        # Beta
        for i in range(beta_electrons):
            orb = mo_coefficients_beta[i]
            orb_atom = np.zeros_like(orb)
            orb_atom[ranges_per_atom[atom][0]:ranges_per_atom[atom][1]] = \
                orb[ranges_per_atom[atom][0]:ranges_per_atom[atom][1]]
            charge_atom += np.dot(orb_atom, np.dot(overlap_matrix, orb))

        charges.append(charge_atom)

    center = np.sum(np.multiply(coordinates.T, charges).T, axis=0)/np.sum(charges)

    return center.tolist()


class WfnSympy:
    def __init__(self,
                 coordinates,  # in Angstrom
                 symbols,
                 basis,  # basis dictionary
                 alpha_mo_coeff,  # Nbas x Nbas
                 center=None,  # in Angstrom
                 VAxis=None,
                 VAxis2=None,
                 charge=0,
                 multiplicity=1,
                 beta_mo_coeff=None,  # Nbas x Nbas
                 group=None,
                 do_operation=False,
                 valence_only=False):

        # Transform group label to igroup, ngroup
        if group is None:
            raise ('point group note defined')
        else:
            igroup, ngroup = _get_group_num_from_label(group)

        type_list_inverse = {}
        for item in shell_type_list.items():
            type_list_inverse['{}'.format(item[1][0])] = int(item[0])

        # from basis dictionary to wfsymm lib arguments
        shell_type = []
        p_exponents = []
        c_coefficients = []
        p_c_coefficients = []
        n_primitives = []
        atom_map = []
        for i, atoms in enumerate(basis['atoms']):
            for shell in atoms['shells']:
                st = shell['shell_type']
                shell_type.append(type_list_inverse[st])
                n_primitives.append(len(shell['p_exponents']))
                atom_map.append(i+1)
                for p in shell['p_exponents']:
                    p_exponents.append(p)
                for c in shell['con_coefficients']:
                    c_coefficients.append(c)
                for pc in shell['p_con_coefficients']:
                    p_c_coefficients.append(pc)

        # Create MO coefficients in contiguous list
        n_molecular_orbitals = len(alpha_mo_coeff)
        Ca = np.ascontiguousarray(alpha_mo_coeff).flatten().tolist()
        if beta_mo_coeff is not None:
            Cb = np.ascontiguousarray(beta_mo_coeff).flatten().tolist()
        else:
            Cb = Ca

        # convert from Angstroms to Bohr
        coordinates = np.array(coordinates)
        coordinates_bohr = np.array(coordinates) / _bohr_to_angstrom

        # get atomic numbers
        atomic_numbers = [symbol_map[i] for i in symbols]

        # get valence electrons
        valence_electrons = get_valence_electrons(atomic_numbers, charge)

        # get total number of electrons
        if valence_only:
            total_electrons = valence_electrons
        else:
            total_electrons = np.sum(atomic_numbers) - charge

        NShell = np.unique(atom_map, return_counts=True)[1]

        # Transform symbols type to correct Fortran char*2 type
        symbols = np.array([list('{:<2}'.format(char)) for char in symbols], dtype='S')

        exp_group = np.array(np.split(np.array(p_exponents), np.cumsum(n_primitives))[:-1])

        Alph=[]
        for i, stype in enumerate(shell_type):
            for _ in range(shell_type_list['{}'.format(stype)][1]):
                Alph.append(exp_group[i])
        Alph2 = []
        for i in Alph:
            Alph2.append(np.ndarray.tolist(i))
        Alph = [item for sublist in Alph2 for item in sublist]

        coef_group = np.array(np.split(np.array(c_coefficients), np.cumsum(n_primitives))[:-1])
        b_coef_group = np.array(np.split(np.array(p_c_coefficients), np.cumsum(n_primitives))[:-1])

        COrb = []
        for i, stype in enumerate(shell_type):
            COrb.append(coef_group[i])
            if 'd' in shell_type_list['{}'.format(stype)][0] or 'f' in shell_type_list['{}'.format(stype)][0]:
                for _ in range(shell_type_list['{}'.format(stype)][1]-1):
                    COrb.append(coef_group[i])
            if shell_type_list['{}'.format(stype)][0] == 'p':
                for _ in range(shell_type_list['{}'.format(stype)][1]-1):
                    COrb.append(coef_group[i])
            if shell_type_list['{}'.format(stype)][0] == 'sp':
                for _ in range(3):
                    COrb.append(b_coef_group[i])
        COrb2 = []
        for i in COrb:
            COrb2.append(np.ndarray.tolist(i))
        COrb = [item for sublist in COrb2 for item in sublist]

        # determine dimension
        Nat = len(symbols)
        NTotShell = len(atom_map)
        Norb = len(COrb)
        NBas = np.sum([shell_type_list['{}'.format(st)][1] for st in shell_type])

        # get overlap matrix
        out = overlap_mat(symbols, coordinates_bohr, NBas, Nat, Norb,NTotShell, NShell,
                          shell_type, n_primitives, COrb, Alph)
        overlap_matrix = np.array(out).reshape(NBas,NBas)

        if center is None:
            # raise Exception('If VAxis, then center should be provided!')
            center = _center_of_charge(alpha_mo_coeff, alpha_mo_coeff,
                                       coordinates, basis,
                                       total_electrons, multiplicity,
                                       overlap_matrix)

        # Axis optimization
        if VAxis is None:
            from wfnsympy.optimize import minimize_axis
            from wfnsympy.optimize import rotation_xy, rotation_axis

            def target_function(alpha, beta, gamma=0.0):
                VAxis = np.dot(rotation_xy(alpha, beta), [1, 0, 0])
                VAxis2 = np.dot(rotation_axis(VAxis, gamma), get_perpendicular_axis(VAxis))

                with _captured_stdout():
                    out_data = mainlib(total_electrons, valence_electrons, NBas, Norb, Nat, NTotShell, atomic_numbers, symbols, Alph,
                                       COrb, NShell, coordinates_bohr, n_primitives, shell_type, igroup, ngroup, Ca, Cb, center, VAxis, VAxis2,
                                       charge, multiplicity, do_operation)

                dgroup = out_data[0][0]
                return np.sum(np.abs(out_data[2][0:dgroup]))

            data = {'coordinates': coordinates_bohr, 'symbols': symbols, 'igroup': igroup, 'ngroup': ngroup}
            alpha, beta, gamma = minimize_axis(target_function, data, delta=0.4)

            VAxis = np.dot(rotation_xy(alpha, beta), [1, 0, 0])
            VAxis2 = np.dot(rotation_axis(VAxis, gamma), get_perpendicular_axis(VAxis))

        else:
            if VAxis2 is None:
                if igroup == 8:  # tetrahedron and octahedron groups
                    from wfnsympy.optimize import minimize_axis2
                    from wfnsympy.optimize import rotation_xy, rotation_axis

                    def target_function(gamma, VAxis):
                        VAxis2 = np.dot(rotation_axis(VAxis, gamma), get_perpendicular_axis(VAxis))

                        with _captured_stdout():
                            out_data = mainlib(total_electrons, valence_electrons, NBas, Norb, Nat, NTotShell, atomic_numbers, symbols, Alph,
                                               COrb, NShell, coordinates_bohr, n_primitives, shell_type, igroup, ngroup, Ca, Cb, center, VAxis, VAxis2,
                                               charge, multiplicity, do_operation)

                        dgroup = out_data[0][0]

                        return np.sum(np.abs(out_data[2][0:dgroup]))

                    gamma = minimize_axis2(target_function, VAxis, delta=0.05)
                    VAxis2 = np.dot(rotation_axis(VAxis, gamma), get_perpendicular_axis(VAxis))

                else:
                    VAxis2 = get_perpendicular_axis(VAxis)

        # Add outputs
        self._center = center
        self._axis = VAxis
        self._axis2 = VAxis2

        with _captured_stdout() as E:
            out_data = mainlib(total_electrons, valence_electrons, NBas, Norb, Nat, NTotShell, atomic_numbers, symbols, Alph,
                               COrb, NShell, coordinates_bohr, n_primitives, shell_type, igroup, ngroup, Ca, Cb, center, VAxis, VAxis2,
                               charge, multiplicity, do_operation)
            E.seek(0)
            capture = E.read()

        capture = capture.decode().split('\n')

        for i, c in enumerate(capture):
            if 'ERROR. Axes not valid' in c:
                self._axis = [float(v) for v in capture[i+3].split()[-3:]]
                self._axis2 = [float(v) for v in capture[i+4].split()[-3:]]


        # Process outputs
        dgroup = out_data[0][0]
        hGroup = out_data[0][1]
        nIR = out_data[0][2]

        self._dgroup = dgroup
        self._atom_coor = np.array(coordinates_bohr)

        self._grim_coef = out_data[1][0:dgroup]
        self._csm_coef = out_data[2][0:dgroup]

        self._SymLab = [''.join(line).strip() for line in np.array(out_data[3][0:dgroup],dtype='str')]

        self._mo_SOEVs_a = out_data[4][:n_molecular_orbitals, 0:dgroup]
        self._mo_SOEVs_b = out_data[5][:n_molecular_orbitals, 0:dgroup]

        self._wf_SOEVs_a = out_data[6][0:dgroup]
        self._wf_SOEVs_b = out_data[7][0:dgroup]
        self._wf_SOEVs = np.prod([self._wf_SOEVs_a, self._wf_SOEVs_b], axis=0)

        self._ideal_gt = out_data[8][:nIR, :dgroup]

        self._IRLab = [''.join(line).strip() for line in np.array(out_data[9][0:nIR],dtype='str')]

        self._mo_IRd_a = out_data[10][:n_molecular_orbitals, 0:nIR]
        self._mo_IRd_b = out_data[11][:n_molecular_orbitals, 0:nIR]

        self._wf_IRd_a = out_data[12][0:nIR]
        self._wf_IRd_b = out_data[13][0:nIR]
        # self._wf_IRd = np.prod([self._wf_IRd_a, self._wf_IRd_b], axis=0)
        self._wf_IRd = out_data[14][0:nIR]

        self._SymMat = out_data[15][0:dgroup]

    # Print Outputs
    def print_CSM(self):
        print('\nWaveFunction: CSM-like values')
        print('     '+'  '.join(['{:^7}'.format(s) for s in self._SymLab]))
        print('Grim' + '  '.join(['{:7.3f}'.format(s) for s in self.grim_coef]))
        print('CSM ' + '  '.join(['{:7.3f}'.format(s) for s in self.csm_coef]))

    def print_overlap_mo_alpha(self):
        print('\nAlpha MOs: Symmetry Overlap Expectation Values')
        print('     '+'  '.join(['{:^7}'.format(s) for s in self._SymLab]))
        for i, line in enumerate(self._mo_SOEVs_a):
            print('{:4d}'.format(i+1) + '  '.join(['{:7.3f}'.format(s) for s in line]))

    def print_overlap_mo_beta(self):
        print('\nBeta MOs: Symmetry Overlap Expectation Values')
        print('     '+'  '.join(['{:^7}'.format(s) for s in self._SymLab]))
        for i, line in enumerate(self._mo_SOEVs_b):
            print('{:4d}'.format(i+1) + '  '.join(['{:7.3f}'.format(s) for s in line]))

    def print_overlap_wf(self):
        print('\nWaveFunction: Symmetry Overlap Expectation Values')
        print('     '+'  '.join(['{:^7}'.format(s) for s in self._SymLab]))
        print('a-wf' + '  '.join(['{:7.3f}'.format(s) for s in self._wf_SOEVs_a]))
        print('b-wf' + '  '.join(['{:7.3f}'.format(s) for s in self._wf_SOEVs_b]))
        print('WFN ' + '  '.join(['{:7.3f}'.format(s) for s in self._wf_SOEVs]))

    def print_ideal_group_table(self):
        print('\nIdeal Group Table')
        print('     '+'  '.join(['{:^7}'.format(s) for s in self._SymLab]))
        for i, line in enumerate(self._ideal_gt):
            print('{:4}'.format(self._IRLab[i]) + '  '.join(['{:7.3f}'.format(s) for s in line]))

    def print_alpha_mo_IRD(self):
        print('\nAlpha MOs: Irred. Rep. Decomposition')
        print('     '+'  '.join(['{:^7}'.format(s) for s in self._IRLab]))
        for i, line in enumerate(self._mo_IRd_a):
            print('{:4d}'.format(i+1) + '  '.join(['{:7.3f}'.format(s) for s in line]))

    def print_beta_mo_IRD(self):
        print('\nBeta MOs: Irred. Rep. Decomposition')
        print('     '+'  '.join(['{:^7}'.format(s) for s in self._IRLab]))
        for i, line in enumerate(self._mo_IRd_b):
            print('{:4d}'.format(i+1) + '  '.join(['{:7.3f}'.format(s) for s in line]))

    def print_wf_mo_IRD(self):
        print('\nWaveFunction: Irred. Rep. Decomposition')
        print('     '+'  '.join(['{:^7}'.format(s) for s in self._IRLab]))
        print('a-wf' + '  '.join(['{:7.3f}'.format(s) for s in self._wf_IRd_a]))
        print('b-wf' + '  '.join(['{:7.3f}'.format(s) for s in self._wf_IRd_b]))
        # print('WFN ' + '  '.join(['{:7.3f}'.format(s) for s in self._IRwf]))

    def print_symmetry_operation_matrix(self, nop):
        if nop > len(self._SymMat):
            print('Not existent')
            exit()

        np.set_printoptions(formatter={'float': '{: 0.8f}'.format})
#        for i, mat in enumerate(self._SymMat):

        print('\n@@@ Operation: {0} {1}'.format(nop+1, self._SymLab[nop]))
        print('Symmetry Transformation matrix')
        print(self._SymMat[nop])

    def print_symmetry_transformed_coordinates(self, nop, use_angstrom=True):

        print('Symmetry Transformed Atomic Coordinates (Angstroms)')
        if use_angstrom:
            print(np.dot(self._atom_coor * _bohr_to_angstrom, self._SymMat[nop].T))
        else:
            print(np.dot(self._atom_coor, self._SymMat[nop].T))

    @property
    def dgroup(self):
        return self._dgroup

    @property
    def grim_coef(self):
        return self._grim_coef

    @property
    def csm_coef(self):
        return self._csm_coef

    @property
    def SymLab(self):
        return self._SymLab

    @property
    def mo_SOEVs_a(self):
        return self._mo_SOEVs_a

    @property
    def mo_SOEVs_b(self):
        return self._mo_SOEVs_b

    @property
    def wf_SOEVs_a(self):
        return self._wf_SOEVs_a

    @property
    def wf_SOEVs_b(self):
        return self._wf_SOEVs_b

    @property
    def wf_SOEVs(self):
        return self._wf_SOEVs

    @property
    def ideal_gt(self):
        return self._ideal_gt

    @property
    def IRLab(self):
        return self._IRLab

    @property
    def mo_IRd_a(self):
        return self._mo_IRd_a

    @property
    def mo_IRd_b(self):
        return self._mo_IRd_b

    @property
    def wf_IRd_a(self):
        return self._wf_IRd_a

    @property
    def wf_IRd_b(self):
        return self._wf_IRd_b

    @property
    def wf_IRd(self):
        return self._wf_IRd

    @property
    def SymMat(self):
        return self._SymMat

    @property
    def center(self):
        return self._center

    @property
    def axis(self):
        return self._axis

    @property
    def axis2(self):
        return self._axis2

symbol_map = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Uut": 113,
    "Uuq": 114,
    "Uup": 115,
    "Uuh": 116,
    "Uus": 117,
    "Uuo": 118,
}
