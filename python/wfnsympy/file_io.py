import numpy as np
from wfnsympy import symbol_map


def reformat_input(array):
    flat_list = []
    for sublist in array:
        for item in sublist:
            if len(item) > 2:
                flat_list.append(item)
            else:
                flat_list.append(item)
    return flat_list


def basis_format(basis_set_name,
                 atomic_numbers,
                 atomic_symbols,
                 shell_type,
                 n_primitives,
                 atom_map,
                 p_exponents,
                 c_coefficients,
                 p_c_coefficients):

    typeList = {'0': ['s', 1],
                '1': ['p', 3],
                '2': ['d', 6],
                '3': ['f', 10],
                '-1': ['sp', 4],
                '-2': ['d_', 5],
                '-3': ['f_', 7]}

    atomic_numbers = np.array(atomic_numbers, dtype=int)
    atom_map = np.array(atom_map, dtype=int)

    basis_set = {'name': basis_set_name,
                 'primitive_type': 'gaussian'}

    shell_type_index = [0] + np.cumsum([typeList['{}'.format(s)][1]
                                        for s in shell_type]).tolist()
    prim_from_shell_index = [0] + np.cumsum(np.array(n_primitives, dtype=int)).tolist()

    atoms_data = []
    for iatom, atomic_number in enumerate(atomic_numbers):
        symbol = atomic_symbols[iatom]

        shell_from_atom_counts = np.unique(atom_map, return_counts=True)[1]
        shell_from_atom_index = np.unique(atom_map, return_index=True)[1]

        shells_data = []
        for ishell in range(shell_from_atom_counts[iatom]):
            st = typeList['{}'.format(shell_type[shell_from_atom_index[iatom] + ishell])]
            ini_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell]
            fin_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell+1]

            shells_data.append({
                'shell_type': st[0],
                'p_exponents': p_exponents[ini_prim: fin_prim],
                'con_coefficients': c_coefficients[ini_prim: fin_prim],
                'p_con_coefficients': p_c_coefficients[ini_prim: fin_prim],
            })

        atoms_data.append({'shells': shells_data,
                           'symbol': symbol,
                           'atomic_number': atomic_number})

    basis_set['atoms'] = atoms_data

    return basis_set


def get_data_from_file_fchk(file_name):

    def convert_to_type(item_type, item):
        item_types = {'I': int,
                      'R': float}

        if type(item) is list:
            return [item_types[item_type](e) for e in item]
        else:
            return item_types[item_type](item)

    key_list = ['Charge', 'Multiplicity', 'Number of alpha electrons', 'Number of beta electrons',
                'Atomic numbers', 'Current cartesian coordinates', 'Number of basis functions', 'Shell types',
                'Number of primitives per shell', 'Shell to atom map', 'Primitive exponents',
                'Contraction coefficients', 'P(S=P) Contraction coefficients', 'Alpha MO coefficients',
                'Beta MO coefficients', 'Coordinates of each shell', 'Overlap Matrix',
                'Core Hamiltonian Matrix', 'Alpha Orbital Energies', 'Beta Orbital Energies',
                'Total SCF Density', 'Spin SCF Density', 'Alpha NATO coefficients', 'Beta NATO coefficients',
                'Alpha Natural Orbital occupancies', 'Beta Natural Orbital occupancies',
                'Natural Transition Orbital occupancies', 'Natural Transition Orbital U coefficients',
                'Natural Transition Orbital V coefficients', 'NTOs occupancies SOC',
                ]

    with open(file_name, mode='r') as f:
        output = f.read()

    basis_set = output.split('\n')[1].split()[-1]
    words_output = output.replace('\n', ' ').split()

    data = {}
    nw = len(words_output)
    for key in key_list:
        wc = len(key.split())
        for i in range(nw):
            word = ' '.join(words_output[i:i+wc])
            if word == key:
                item_type = words_output[i+wc]
                for l_step in range(5):
                    try:
                        if words_output[l_step + i + wc + 1] == 'N=':
                            item_type = words_output[l_step + i + wc]
                            word += ' '.join(words_output[i + wc:l_step + i + wc])
                            n_elements = int(words_output[l_step + i + wc + 2])
                            data[word] = convert_to_type(item_type, words_output[l_step + i + wc + 3: l_step + i + wc + n_elements + 3])
                        else:
                            item_type = words_output[l_step + i + wc]
                            data[word] = convert_to_type(item_type, words_output[l_step + i + wc + 1])
                        break
                    except KeyError:
                        pass
                break

    bohr_to_angstrom = 0.529177249

    coordinates = np.array(data['Current cartesian coordinates']).reshape(-1, 3) * bohr_to_angstrom

    if not 'P(S=P) Contraction coefficients' in data:
        data['P(S=P) Contraction coefficients'] = np.zeros_like(data['Contraction coefficients']).tolist()

    def get_symbols_from_atomic_numbers(atomic_numbers):
        symbols = []
        for n in atomic_numbers:
            for symbol, n2 in symbol_map.items():
                if int(n2) == n:
                    symbols.append(symbol)
        return symbols

    symbols = get_symbols_from_atomic_numbers(data['Atomic numbers'])
    basis = basis_format(basis_set_name=basis_set,
                         atomic_numbers=data['Atomic numbers'],
                         atomic_symbols=symbols,
                         shell_type=data['Shell types'],
                         n_primitives=data['Number of primitives per shell'],
                         atom_map=data['Shell to atom map'],
                         p_exponents=data['Primitive exponents'],
                         c_coefficients=data['Contraction coefficients'],
                         p_c_coefficients=data['P(S=P) Contraction coefficients'])

    nbas = data['Number of basis functions']

    final_dict = {'basis': basis,
                  'number_of_electrons': {'alpha': data['Number of alpha electrons'],
                                          'beta': data['Number of beta electrons']}
                  }

    if 'Alpha MO coefficients' in data:
        final_dict['coefficients'] = {'alpha': np.array(data['Alpha MO coefficients']).reshape(-1, nbas).tolist()}
        final_dict['mo_energies'] = {'alpha': data['Alpha Orbital Energies']}

    if 'Beta MO coefficients' in data:
        final_dict['coefficients']['beta'] = np.array(data['Beta MO coefficients']).reshape(-1, nbas).tolist()
        final_dict['mo_energies']['beta'] = data['Beta Orbital Energies']

    # wf orbital configuration
    alpha_occupancy = [1] * int(data['Number of alpha electrons'])
    beta_occupancy = [1] * int(data['Number of beta electrons'])
    if len(alpha_occupancy) != len(beta_occupancy):
        for i in range(abs(len(alpha_occupancy) - len(beta_occupancy))):
            beta_occupancy.append(0)

    return {'coordinates': coordinates,
            'symbols': symbols,
            'basis': basis,
            'mo_coefficients': final_dict['coefficients'],
            'alpha_occupancy': alpha_occupancy,
            'beta_occupancy': beta_occupancy}


if __name__ == '__main__':
    data = get_data_from_file_fchk('pirrol.in.fchk')
    print(data)
