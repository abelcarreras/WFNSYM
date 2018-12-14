from wfnsympy.WFNSYMLIB import mainlib
import numpy as np


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


def get_group_num_from_label(label):

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
            label_2[0] = 'C'

    operations = {'Cn':  [1, ngroup],
                  'CnH': [2, ngroup],
                  'CnV': [3, ngroup],
                  'CI':  [0, 2],
                  'CINF':[9, 1],
                  'Dn':  [4, ngroup],
                  'Dnd': [5, ngroup],
                  'DnH': [6, ngroup],
                  'DINF':[9, 2],
                  'Sn':  [7, ngroup],
                  'T':   [8, 1],
                  'TH':  [8, 2],
                  'TD':  [8, 3],
                  'O':   [8, 4],
                  'OH':  [8, 5],
                  'I':   [8, 6],
                  'IH':  [8, 7],
                  }
    try:
        return operations[label_2]
    except KeyError:
        raise Exception('Label not found')


def get_operation_num_from_label(label):
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


class WfnSympy:
    def __init__(self,
                 coordinates,  # in Angstrom
                 symbols,
                 basis,  # basis dictionary
                 center,  # in Angstrom
                 VAxis,
                 VAxis2,
                 alpha_mo_coeff,  # Nbas x Nbas
                 charge=0,
                 multiplicity=1,
                 beta_mo_coeff=None,  # Nbas x Nbas
                 group=None,
                 do_operation=False,
                 use_pure_d_functions=False):

        # Transform group label to igroup, ngroup
        if group is None:
            raise ('point group note defined')
        else:
            igroup, ngroup = get_group_num_from_label(group)

        bohr_to_angstrom = 0.529177249

        shell_type_list = {'-1': ['sp', 4],
                            '0': ['s',  1],
                            '1': ['p',  3],
                            '2': ['d',  6],
                            '3': ['f',  10],
                           '-2': ['dc', 5],  # pure
                           '-3': ['fc', 7]}  # pure

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
        Ca = np.ascontiguousarray(alpha_mo_coeff).flatten().tolist()
        if beta_mo_coeff is not None:
            Cb = np.ascontiguousarray(beta_mo_coeff).flatten().tolist()
        else:
            Cb = Ca

        # convert from Angstroms to Bohr
        coordinates = np.array(coordinates)/bohr_to_angstrom
        center = [c/bohr_to_angstrom for c in center]

        # get atomic numbers
        atomic_numbers = [symbol_map[i] for i in symbols]

        # get valence electrons
        valence_electrons = get_valence_electrons(atomic_numbers, charge)

        # get total number of electrons
        total_electrons = np.sum(atomic_numbers) - charge

        NShell = np.unique(atom_map, return_counts=True)[1]

        # Transform symbols type to correct Fortran char*2 type
        symbols = np.array([list('{:<2}'.format(char)) for char in symbols], dtype='S')

        # Check for pure D
        if use_pure_d_functions:
            shell_type_list.update({'-2' : ['d', 5],
                             '-3' : ['f', 7]})

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
            if shell_type_list['{}'.format(stype)][0] == 'd' or shell_type_list['{}'.format(stype)][0] == 'f':
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

        out_data = mainlib(total_electrons, valence_electrons, NBas, Norb, Nat, NTotShell, atomic_numbers, symbols, Alph,
                           COrb, NShell, coordinates, n_primitives, shell_type, igroup, ngroup, Ca, Cb, center, VAxis, VAxis2,
                           charge, multiplicity, do_operation, use_pure_d_functions)

        # Process outputs
        dgroup = out_data[0][0]
        hGroup = out_data[0][1]
        nIR = out_data[0][2]

        self._dgroup = dgroup
        self._atom_coor = np.array(coordinates)

        self._grim_coef = out_data[1][0:dgroup]
        self._csm_coef = out_data[2][0:dgroup]

        self._SymLab = [''.join(line).strip() for line in np.array(out_data[3][0:dgroup],dtype='str')]

        self._mo_SOEVs_a = out_data[4][:, 0:dgroup]
        self._mo_SOEVs_b = out_data[5][:, 0:dgroup]

        self._wf_SOEVs_a = out_data[6][0:dgroup]
        self._wf_SOEVs_b = out_data[7][0:dgroup]
        self._wf_SOEVs = np.prod([self._wf_SOEVs_a, self._wf_SOEVs_b], axis=0)

        self._ideal_gt = out_data[8][:nIR, :dgroup]

        self._IRLab = [''.join(line).strip() for line in np.array(out_data[9][0:nIR],dtype='str')]

        self._mo_IRd_a = out_data[10][:, 0:nIR]
        self._mo_IRd_b = out_data[11][:, 0:nIR]

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
            print(np.dot(self._atom_coor * 0.529177249, self._SymMat[nop].T))
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

if __name__ == '__main__':

    # Atomic elements
    AtLab = ['H', 'N', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H']

    # Shell type [sp:-1, s:0, p:1, d:2, f:3]
    shell_type = [0, 0, -1, 0, -1, 0, -1, 0, -1, 0, -1, 0, 0, 0, 0]

    # Basis exponents
    p_exp = [3.42525091E+00, 6.23913730E-01, 1.68855400E-01, 9.91061690E+01, 1.80523120E+01,
             4.88566020E+00, 3.78045590E+00, 8.78496600E-01, 2.85714400E-01, 7.16168370E+01,
             1.30450960E+01, 3.53051220E+00, 2.94124940E+00, 6.83483100E-01, 2.22289900E-01,
             7.16168370E+01, 1.30450960E+01, 3.53051220E+00, 2.94124940E+00, 6.83483100E-01,
             2.22289900E-01, 7.16168370E+01, 1.30450960E+01, 3.53051220E+00, 2.94124940E+00,
             6.83483100E-01, 2.22289900E-01, 7.16168370E+01, 1.30450960E+01, 3.53051220E+00,
             2.94124940E+00, 6.83483100E-01, 2.22289900E-01, 3.42525091E+00, 6.23913730E-01,
             1.68855400E-01, 3.42525091E+00, 6.23913730E-01, 1.68855400E-01, 3.42525091E+00,
             6.23913730E-01, 1.68855400E-01, 3.42525091E+00, 6.23913730E-01, 1.68855400E-01]

    # Basis contraction coefficients
    con_coef = [1.54328971E-01,  5.35328142E-01,  4.44634542E-01,  1.54328970E-01,  5.35328141E-01,
                4.44634541E-01, -9.99672284E-02,  3.99512824E-01,  7.00115459E-01,  1.54328970E-01,
                5.35328140E-01,  4.44634540E-01, -9.99672301E-02,  3.99512830E-01,  7.00115471E-01,
                1.54328970E-01,  5.35328140E-01,  4.44634540E-01, -9.99672301E-02,  3.99512830E-01,
                7.00115471E-01,  1.54328970E-01,  5.35328140E-01,  4.44634540E-01, -9.99672301E-02,
                3.99512830E-01,  7.00115471E-01,  1.54328970E-01,  5.35328140E-01,  4.44634540E-01,
                -9.99672301E-02,  3.99512830E-01,  7.00115471E-01,  1.54328971E-01,  5.35328142E-01,
                4.44634542E-01,  1.54328971E-01,  5.35328142E-01,  4.44634542E-01,  1.54328971E-01,
                5.35328142E-01,  4.44634542E-01,  1.54328971E-01,  5.35328142E-01,  4.44634542E-01]

    # Basis P(S=P) Contraction coefficients
    p_con_coef = [0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
                  0.00000000E+00, 1.55916269E-01, 6.07683714E-01, 3.91957386E-01, 0.00000000E+00,
                  0.00000000E+00, 0.00000000E+00, 1.55916272E-01, 6.07683728E-01, 3.91957395E-01,
                  0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 1.55916272E-01, 6.07683728E-01,
                  3.91957395E-01, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 1.55916272E-01,
                  6.07683728E-01, 3.91957395E-01, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
                  1.55916272E-01, 6.07683728E-01, 3.91957395E-01, 0.00000000E+00, 0.00000000E+00,
                  0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,
                  0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00]

    # Atomic coordinates (in Bohr)
    RAt = [[0.00000000E+00,  0.00000000E+00, 4.02378628E+00],
           [0.00000000E+00,  0.00000000E+00,  2.11400638E+00],
           [0.00000000E+00,  2.12502726E+00,  6.30346495E-01],
           [0.00000000E+00, -2.12502726E+00,  6.30346495E-01],
           [0.00000000E+00,  1.33869522E+00, -1.85912579E+00],
           [0.00000000E+00, -1.33869522E+00, -1.85912579E+00],
           [0.00000000E+00,  3.99274942E+00,  1.45602642E+00],
           [0.00000000E+00, -3.99274942E+00,  1.45602642E+00],
           [0.00000000E+00,  2.56483456E+00, -3.49426423E+00],
           [0.00000000E+00, -2.56483456E+00, -3.49426423E+00]]

    # Number of primitives per shell
    n_prim = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

    # Shell to atom map
    atom_map = [1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10]

    # Alpha MO coefficients
    Ca = [-6.43302924E-03,  9.93234079E-01,  3.29472273E-02, -2.01329290E-16, -1.57729780E-16,
          -1.08727169E-03,  4.82521106E-04, -6.11453036E-03, -8.23844925E-17,  3.92413081E-03,
          -2.68062978E-03,  4.82521106E-04, -6.11453036E-03,  2.93323763E-16, -3.92413081E-03,
          -2.68062978E-03,  3.78711048E-06,  4.03033224E-04,  6.90004209E-17, -4.47950358E-05,
           2.62913916E-04,  3.78711049E-06,  4.03033224E-04, -1.51458678E-16,  4.47950358E-05,
           2.62913916E-04,  2.41721138E-04,  2.41721138E-04, -4.12471347E-05, -4.12471347E-05,
          -2.53652231E-14,  3.61875068E-14, -5.79059146E-13, -1.33654376E-16,  4.53714015E-03,
           1.94507892E-13, -7.01571640E-01, -2.51577403E-02, -1.00547184E-16,  1.69051610E-03,
           3.67122644E-04,  7.01571640E-01,  2.51577403E-02,  2.22803447E-17,  1.69051611E-03,
          -3.67122644E-04, -4.03720687E-03,  4.96651647E-03,  1.80708799E-16,  1.04632691E-03,
           2.78774385E-03,  4.03720687E-03, -4.96651647E-03,  3.12255808E-16,  1.04632691E-03,
          -2.78774385E-03,  4.56523911E-03, -4.56523911E-03, -1.42718181E-04,  1.42718181E-04,
           3.62850897E-04, -5.38387755E-04,  8.61533062E-03, -2.76989842E-16,  3.05177421E-13,
          -2.89323068E-03, -7.01721658E-01, -2.48092649E-02,  6.69117588E-17,  9.95636542E-04,
           1.91439843E-04, -7.01721658E-01, -2.48092649E-02, -1.19427180E-17, -9.95636542E-04,
           1.91439843E-04, -2.42380063E-03,  4.82834345E-03, -5.35437333E-16,  8.87507705E-04,
           3.13609281E-03, -2.42380063E-03,  4.82834345E-03,  2.54885263E-16, -8.87507704E-04,
           3.13609281E-03,  4.76032662E-03,  4.76032662E-03,  4.01901906E-05,  4.01901906E-05,
          -1.28587463E-04,  4.13483178E-05,  4.08215038E-04,  2.11504044E-16,  1.16436842E-14,
           1.29534150E-04, -1.83139841E-03, -5.26404865E-03, -2.03870677E-17,  1.12735092E-03,
           3.11489423E-03, -1.83139841E-03, -5.26404865E-03, -1.60498233E-16, -1.12735092E-03,
           3.11489423E-03,  7.01836747E-01,  2.09602195E-02,  2.54792493E-16,  2.53938918E-03,
           1.45204349E-03,  7.01836746E-01,  2.09602195E-02,  2.96226298E-17, -2.53938918E-03,
           1.45204349E-03,  2.05746794E-05,  2.05746793E-05, -4.62177091E-03, -4.62177090E-03,
           8.14332465E-15, -1.23582329E-15, -2.58259753E-14,  1.19915780E-16,  1.83190654E-04,
          -9.13669085E-15, -3.61380247E-03, -5.51150989E-03,  1.20508258E-16,  8.14933647E-04,
           3.17237966E-03,  3.61380247E-03,  5.51150989E-03, -1.88758570E-16,  8.14933647E-04,
          -3.17237966E-03,  7.01247540E-01,  3.07593737E-02,  1.69916876E-16, -3.83023453E-03,
           1.37592642E-03, -7.01247540E-01, -3.07593737E-02, -3.83779654E-17, -3.83023453E-03,
          -1.37592642E-03,  1.84548290E-04, -1.84548290E-04, -4.76355663E-03,  4.76355663E-03,
          -9.87227083E-02,  1.91980571E-01, -6.20919845E-01,  7.27868211E-16,  1.07405441E-15,
           6.45353701E-02,  9.60305337E-02, -2.16218583E-01, -9.16100402E-18,  8.65504006E-02,
          -2.58799110E-02,  9.60305337E-02, -2.16218583E-01, -7.01902907E-16, -8.65504006E-02,
          -2.58799110E-02,  4.99675351E-02, -1.12530932E-01,  1.18991733E-16,  1.72383588E-02,
          -4.50005092E-02,  4.99675351E-02, -1.12530932E-01, -3.64081666E-16, -1.72383588E-02,
          -4.50005092E-02, -3.32679681E-02, -3.32679681E-02, -1.41659696E-02, -1.41659696E-02,
          -1.52474976E-01,  1.06625049E-01, -3.86426945E-01,  1.91650777E-15, -5.06832248E-15,
          -1.27740858E-01, -4.23886032E-02,  1.12037452E-01, -5.94821673E-16, -6.96622093E-04,
          -1.40379868E-01, -4.23886032E-02,  1.12037452E-01, -3.47237183E-15,  6.96622093E-04,
          -1.40379868E-01, -1.46706154E-01,  3.83218597E-01,  2.19318684E-15, -7.19253196E-02,
           4.28410360E-02, -1.46706154E-01,  3.83218597E-01,  8.42606668E-17,  7.19253196E-02,
           4.28410360E-02,  1.54930240E-02,  1.54930240E-02,  8.83783298E-02,  8.83783298E-02,
          -5.27414984E-15,  2.85135345E-15, -1.01300878E-14,  3.59232659E-15,  2.58757568E-01,
          -7.11720794E-15, -1.53701535E-01,  4.30098966E-01, -1.25736239E-15, -1.05768165E-02,
          -1.48467254E-02,  1.53701535E-01, -4.30098966E-01, -8.30711187E-15, -1.05768165E-02,
           1.48467254E-02, -7.92143193E-02,  2.16997014E-01,  5.20233090E-15,  8.12968962E-02,
           8.21457504E-02,  7.92143193E-02, -2.16997014E-01,  4.63305215E-16,  8.12968962E-02,
          -8.21457504E-02,  1.36088025E-01, -1.36088025E-01,  6.33246216E-02, -6.33246216E-02,
           3.09730535E-01, -4.08311976E-02,  1.60680263E-01, -2.69541242E-15,  2.10354035E-14,
           3.93761943E-01,  8.45631440E-02, -2.53969493E-01, -1.43411653E-15,  8.32585529E-03,
          -7.46710104E-03,  8.45631440E-02, -2.53969493E-01, -9.36458382E-16, -8.32585529E-03,
          -7.46710104E-03, -7.85554216E-02,  2.42253288E-01, -1.14765020E-15, -7.58250388E-02,
          -1.77611778E-01, -7.85554216E-02,  2.42253288E-01,  3.17003423E-16,  7.58250388E-02,
          -1.77611778E-01, -1.04232034E-01, -1.04232034E-01,  1.53519809E-01,  1.53519809E-01,
          -1.61558641E-14,  2.89187550E-15, -1.11064580E-14,  4.38621721E-15,  3.25359542E-01,
          -1.84587426E-14, -3.52357291E-02,  1.08412205E-01,  1.45054825E-15, -4.32407532E-03,
           2.52986861E-01,  3.52357291E-02, -1.08412205E-01, -4.39071703E-15, -4.32407532E-03,
          -2.52986861E-01,  1.11900169E-01, -3.33930897E-01,  2.53924947E-15, -1.49063486E-01,
           3.06918885E-02, -1.11900169E-01,  3.33930897E-01, -7.04654860E-16, -1.49063486E-01,
          -3.06918885E-02,  1.04427013E-01, -1.04427013E-01, -2.12561486E-01,  2.12561486E-01,
           1.18637065E-01,  3.53492532E-02, -1.38772990E-01,  2.36708177E-15, -7.62060586E-15,
           2.42217798E-01, -6.85634762E-02,  2.26672551E-01,  5.38386806E-16,  2.90592539E-01,
           3.47267189E-02, -6.85634762E-02,  2.26672551E-01, -1.86951531E-15, -2.90592539E-01,
           3.47267189E-02, -4.33571840E-04,  1.40031671E-02,  1.12674610E-15,  1.52989239E-01,
          -3.75042598E-02, -4.33571840E-04,  1.40031671E-02, -1.02148943E-15, -1.52989239E-01,
          -3.75042598E-02,  2.84990287E-01,  2.84990287E-01,  7.98345397E-02,  7.98345397E-02,
           3.02305419E-01,  3.46237685E-02, -1.58523145E-01, -4.08134489E-15,  7.11426183E-14,
           3.90319421E-01, -6.47154535E-03,  1.20998258E-02, -3.76085485E-15,  4.21464128E-02,
          -2.30319466E-01, -6.47154535E-03,  1.20998258E-02, -1.67465968E-14, -4.21464128E-02,
          -2.30319466E-01, -1.43010385E-02,  3.80688641E-02,  3.83314868E-15, -1.17725585E-01,
           3.21936688E-01, -1.43010385E-02,  3.80688641E-02, -5.43152452E-15,  1.17725585E-01,
           3.21936688E-01, -2.79318514E-02, -2.79318514E-02, -2.61632159E-01, -2.61632159E-01,
          -3.98991149E-14, -4.98050664E-15,  2.24736896E-14,  6.82437704E-14,  4.44540750E-01,
          -5.29433180E-14,  2.85551141E-03, -9.07970437E-03,  3.50813626E-14, -3.92062267E-01,
           8.15055083E-03, -2.85551141E-03,  9.07970437E-03,  3.42049539E-14, -3.92062267E-01,
          -8.15055083E-03, -2.89789848E-02,  9.40426732E-02,  2.34780485E-14,  6.90297658E-03,
           1.93614032E-02,  2.89789848E-02, -9.40426732E-02,  2.27698239E-14,  6.90297658E-03,
          -1.93614032E-02, -3.42451303E-01,  3.42451303E-01,  4.11409617E-02, -4.11409617E-02,
          -5.78535688E-15, -6.03218406E-16,  3.17817516E-15, -6.24961070E-01,  4.66335776E-14,
          -6.90145839E-15,  1.32960693E-16, -5.64372744E-16, -3.38221758E-01, -4.17915248E-14,
           4.36632950E-15, -1.86470556E-16,  7.48432599E-17, -3.38221758E-01, -4.14993600E-14,
           2.45458286E-15, -2.41828259E-15,  8.32986431E-15, -2.21714606E-01,  1.30842300E-15,
          -1.16157678E-15,  3.27618593E-15, -1.06675017E-14, -2.21714606E-01, -1.92979413E-15,
          -7.09022591E-15, -3.47349147E-14,  3.78393535E-14,  7.00885432E-15,  1.67496278E-16,
           1.54581495E-14,  2.63408194E-15, -1.31017103E-14,  1.42197693E-15, -1.93294543E-01,
           1.99642207E-14, -4.30379436E-02,  1.47542456E-01,  5.62206179E-16, -6.52120835E-02,
          -3.01812607E-01,  4.30379436E-02, -1.47542456E-01, -7.45054504E-16, -6.52120835E-02,
           3.01812607E-01,  5.08122345E-03, -2.42517071E-02,  2.91133677E-15, -5.41873748E-02,
           3.66808951E-01, -5.08122345E-03,  2.42517071E-02,  4.58084304E-16, -5.41873748E-02,
          -3.66808951E-01, -8.12858345E-02,  8.12858345E-02, -3.29326362E-01,  3.29326362E-01,
           7.22286478E-02, -1.49902292E-02,  6.92133075E-02,  2.02070351E-15, -2.60687851E-15,
           3.34520416E-02, -6.05091634E-03,  1.96930852E-02,  1.15453837E-15, -1.11195101E-01,
          -2.44650318E-01, -6.05091634E-03,  1.96930852E-02, -9.56427735E-15,  1.11195101E-01,
          -2.44650318E-01,  1.39855886E-02, -5.47561349E-02,  6.83464456E-15,  4.66859846E-01,
           9.26608564E-02,  1.39855886E-02, -5.47561349E-02, -5.88511323E-16, -4.66859846E-01,
           9.26608564E-02, -2.02306447E-01, -2.02306447E-01,  1.72324227E-01,  1.72324227E-01,
           2.42151660E-16, -1.31941667E-16,  4.23343537E-16,  6.02724371E-01, -6.06446005E-16,
          -2.51255463E-16, -1.10111243E-16, -6.35861725E-17, -6.87514109E-02, -2.63090180E-16,
          -2.04588453E-15,  1.07900029E-16,  2.35439262E-16, -6.87514109E-02,  1.30663903E-15,
           3.83057063E-16, -2.01864113E-16,  5.29043461E-16, -5.22154476E-01,  8.79433321E-16,
           5.63057717E-17,  1.38486278E-16, -6.33300943E-16, -5.22154476E-01, -7.85152155E-16,
          -4.85973637E-16,  1.40960301E-16, -1.58835179E-15, -4.01068223E-16,  7.72971096E-16,
          -2.84437612E-15, -5.38508313E-16,  1.83894485E-15,  4.33196692E-15, -2.28305265E-15,
          -5.39976928E-15,  1.62405179E-15, -4.39663220E-15,  5.84599406E-01,  1.32569081E-15,
           4.47291632E-15, -9.39151727E-16,  2.82225608E-15, -5.84599406E-01,  7.27474341E-16,
           6.96027644E-15,  1.92091766E-16,  3.53657286E-16,  3.48699453E-01, -3.38693461E-15,
          -6.35335712E-15,  3.37098411E-16, -8.60243973E-16, -3.48699453E-01,  2.57242115E-15,
          -2.67558238E-15,  1.77972050E-16,  2.01706341E-15,  1.06172959E-15,  5.91541220E-17,
          -1.24572600E-15, -4.64060973E-16,  1.65200811E-15,  5.66545451E-01,  3.74386466E-16,
          -5.06523932E-15,  1.83136521E-15, -4.82142053E-15, -6.53361703E-01,  4.59403728E-15,
           2.11707894E-15, -1.26781061E-15,  4.21989519E-15, -6.53361703E-01,  1.14760141E-15,
           8.06104489E-15, -2.99062395E-17,  3.16990092E-15,  3.40562871E-01, -5.35833337E-15,
          -6.30278200E-15,  4.75557112E-16, -2.12246000E-15,  3.40562871E-01, -3.24734416E-15,
           7.81385303E-16, -1.26158644E-15, -3.37220290E-16, -1.80793503E-15,  5.58708296E-17,
           9.57044617E-16,  2.43667967E-16, -1.55940348E-15,  2.10697156E-15,  8.94638085E-16,
          -5.37051572E-16,  2.78810535E-16, -1.26788039E-15, -4.49180282E-01,  2.96018270E-16,
           1.65128623E-15, -6.52397879E-17,  4.82070819E-16,  4.49180282E-01,  9.63457697E-16,
           9.63716851E-16,  1.52605611E-16, -7.62438927E-16,  7.33205990E-01,  5.02776426E-16,
           9.80532127E-16, -1.33604404E-16,  7.40194380E-16, -7.33205990E-01,  1.28580348E-15,
           3.56893453E-17, -2.10339124E-16,  2.65354486E-16,  1.42405635E-15,  4.76303502E-16,
          -9.44301352E-01, -1.40147681E-01,  9.56135465E-01, -2.63627945E-16, -5.77206453E-16,
           4.07709349E-01,  1.31887472E-03,  3.04278097E-02, -3.88806284E-17,  3.68260415E-01,
          -9.05590178E-02,  1.31887472E-03,  3.04278097E-02,  8.83916152E-16, -3.68260415E-01,
          -9.05590178E-02, -4.65510596E-02,  2.88165857E-01,  8.75454441E-16,  1.18403033E-01,
          -1.71691812E-01, -4.65510596E-02,  2.88165857E-01, -9.38367770E-16, -1.18403033E-01,
          -1.71691812E-01, -2.91241721E-01, -2.91241721E-01, -4.11315420E-01, -4.11315420E-01,
           2.81933934E-01,  1.25183468E-01, -8.32281803E-01,  4.11669294E-16,  2.50623342E-14,
           1.39864296E-01, -1.46437320E-01,  8.96279971E-01,  8.08265399E-16,  7.94012079E-02,
           3.56306922E-01, -1.46437320E-01,  8.96279971E-01,  1.32703143E-16, -7.94012079E-02,
           3.56306922E-01,  5.52884358E-03, -1.00574876E-02, -1.18634023E-15,  1.12349364E-02,
          -1.88675909E-01,  5.52884358E-03, -1.00574876E-02,  1.75654188E-16, -1.12349364E-02,
          -1.88675909E-01, -7.22432561E-01, -7.22432561E-01, -1.34185479E-01, -1.34185479E-01,
          -3.56405628E-15, -1.08175962E-14,  7.26687420E-14, -6.56671767E-16, -4.37588158E-01,
          -3.10886198E-14,  3.00569412E-02, -2.14179919E-01, -1.20070439E-17, -6.44731852E-01,
           8.15466885E-02, -3.00569412E-02,  2.14179919E-01,  1.01406707E-15, -6.44731852E-01,
          -8.15466885E-02, -8.19230493E-02,  4.91559873E-01,  6.16568618E-16, -4.62576616E-01,
           2.26485105E-02,  8.19230493E-02, -4.91559873E-01, -1.03560023E-15, -4.62576616E-01,
          -2.26485105E-02,  6.32647306E-01, -6.32647306E-01,  2.15978545E-02, -2.15978545E-02,
           6.58374792E-01, -9.05012130E-03,  4.53101937E-02,  1.57763526E-16, -3.96835175E-14,
          -7.15346696E-01,  2.42174657E-03,  5.15564912E-03,  4.78314356E-17,  1.14080619E-01,
          -4.28921872E-01,  2.42174657E-03,  5.15564912E-03,  7.25896812E-16, -1.14080619E-01,
          -4.28921872E-01, -1.44408513E-02,  1.02533077E-01,  4.72880826E-16,  1.17603319E-01,
          -5.52776826E-01, -1.44408513E-02,  1.02533077E-01, -7.20446824E-16, -1.17603319E-01,
          -5.52776826E-01,  6.12785888E-02,  6.12785888E-02, -5.49535483E-01, -5.49535483E-01,
           9.64410379E-15, -1.10562351E-14,  7.41768670E-14, -6.34317469E-16,  5.89462542E-01,
          -4.50346414E-14,  1.22689317E-01, -7.57976835E-01,  1.22525220E-15,  1.68048080E-01,
           2.65255455E-03, -1.22689317E-01,  7.57976835E-01,  1.01738188E-15,  1.68048080E-01,
          -2.65255455E-03, -1.29631911E-01,  8.07440324E-01, -3.23275593E-16,  2.66082323E-01,
          -8.16185537E-02,  1.29631911E-01, -8.07440324E-01, -8.44543322E-17,  2.66082323E-01,
           8.16185537E-02,  2.57018163E-01, -2.57018163E-01, -6.43523132E-01,  6.43523132E-01,
          -2.46728270E-14, -6.57128293E-15,  5.05276259E-14, -9.36676333E-16, -5.47429488E-01,
           9.94844031E-15, -9.12357471E-02,  6.08721174E-01,  3.82345758E-16,  1.41524979E-01,
           3.63735398E-01,  9.12357471E-02, -6.08721174E-01,  1.94410340E-15,  1.41524979E-01,
          -3.63735398E-01, -9.95754667E-02,  6.56017346E-01, -2.92882015E-16, -2.32770919E-01,
          -4.12488540E-01,  9.95754667E-02, -6.56017346E-01, -8.54779592E-16, -2.32770919E-01,
           4.12488540E-01, -5.14170663E-01,  5.14170663E-01, -4.74524870E-01,  4.74524870E-01,
           7.86785260E-02,  5.15017170E-02, -3.78696865E-01, -3.62076446E-16, -8.32130994E-14,
           5.16499264E-02,  6.33473193E-02, -4.29606105E-01,  5.31659327E-16, -2.78061558E-01,
           5.65227547E-01,  6.33473193E-02, -4.29606105E-01,  9.51290367E-16,  2.78061558E-01,
           5.65227547E-01, -1.16435136E-01,  7.83286642E-01, -6.98987426E-16,  4.96994034E-01,
           1.18062746E-01, -1.16435136E-01,  7.83286642E-01, -1.95386955E-16, -4.96994034E-01,
           1.18062746E-01,  2.20815063E-01,  2.20815063E-01, -5.23017249E-01, -5.23017249E-01,
           3.13353886E-15, -1.72667826E-15,  1.20086099E-14, -2.93469300E-15,  5.80918359E-01,
          -8.22503154E-15, -1.23098022E-02,  1.13203565E-01,  2.67133894E-15,  2.36579617E-01,
          -6.78028300E-01,  1.23098022E-02, -1.13203565E-01,  4.59057681E-15,  2.36579617E-01,
           6.78028300E-01, -6.37606033E-03,  3.39466967E-02, -1.61775992E-15, -7.71131019E-01,
          -4.30069769E-01,  6.37606033E-03, -3.39466967E-02, -1.58804183E-15, -7.71131019E-01,
           4.30069769E-01,  3.74665188E-03, -3.74665188E-03,  7.00133348E-02, -7.00133348E-02,
           3.93242760E-01, -5.92389799E-02,  4.47556943E-01, -1.25095887E-15, -6.88716547E-15,
          -7.97989310E-01,  6.80756222E-02, -4.78953465E-01,  1.33121561E-15,  7.72627735E-01,
           3.21666456E-01,  6.80756222E-02, -4.78953465E-01,  1.77482190E-15, -7.72627735E-01,
           3.21666456E-01, -3.94285484E-02,  2.94209040E-01, -2.18869230E-15, -6.28208479E-03,
           5.10184070E-01, -3.94285484E-02,  2.94209040E-01, -1.61961443E-16,  6.28208479E-03,
           5.10184070E-01, -3.95073008E-01, -3.95073008E-01,  1.79430894E-01,  1.79430894E-01,
           2.04857969E-14, -2.75854214E-15,  2.12835102E-14, -1.46066047E-15,  1.67792226E-01,
          -4.47653089E-14,  5.53410142E-02, -4.20691608E-01,  7.20731658E-16,  5.04907302E-01,
           4.33925108E-01, -5.53410142E-02,  4.20691608E-01,  1.09676849E-15,  5.04907302E-01,
          -4.33925108E-01, -8.69435277E-02,  6.56193193E-01, -1.16218050E-15, -5.59207294E-01,
           7.88398567E-01,  8.69435277E-02, -6.56193193E-01, -5.75068294E-17, -5.59207294E-01,
          -7.88398567E-01, -2.67731796E-01,  2.67731796E-01,  4.07955748E-01, -4.07955748E-01]

    pirrol = WfnSympy(#Etot=36,        # Number of electrons
                      NEval=26,       # Number of Valence electrons
                      # NBas=30,      # Number of basis function
                      # Norb=45*2,    # Number of Orbitals
                      # NTotShell=15, # Number of shells
                      # iZAt=iZAt,
                      AtLab=AtLab,    # Atomic labels
                      shell_type=shell_type,
                      p_exp=p_exp,
                      con_coef=con_coef,
                      p_con_coef=p_con_coef,
                      RAt=RAt,  # atomic coordinates in Bohr
                      n_prim=n_prim,
                      atom_map=atom_map,
                      Ca=Ca, Cb=Ca,
                      RCread=[0., 0., 0.],
                      VAxis= [1., 0., 0.],
                      VAxis2=[0., 0., 1.],
                      charge=0, multiplicity=1,
                      # igroup=3, ngroup=6,  # define symmetry group
                      group='c6v',
                      do_operation=False)

    # Testing
    pirrol.print_alpha_mo_IRD()
    pirrol.print_beta_mo_IRD()
    pirrol.print_wf_mo_IRD()
    pirrol.print_CSM()
    pirrol.print_ideal_group_table()
    pirrol.print_overlap_mo_alpha()
    pirrol.print_overlap_mo_beta()
    pirrol.print_overlap_wf()
    for i in range(pirrol.dgroup):
        pirrol.print_symmetry_operation_matrix(i)
        pirrol.print_symmetry_transformed_coordinates(i)
