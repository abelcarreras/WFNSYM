class MultiplicityError(Exception):
    def __init__(self, multiplicity, total_electrons):
        self._multiplicity = multiplicity
        self._electrons = total_electrons

    def __str__(self):
        return 'Multiplicity {} incompatible with {} electrons'.format(self._multiplicity,
                                                                       self._electrons)


class LabelNotFound(Exception):
    def __init__(self, label):
        self._label = label

    def __str__(self):
        return 'Symmetry label {} not defined'.format(self._label)


class ChangedAxisWarning(UserWarning):
    def __init__(self, axis, axis2):
        self.axis = axis
        self.axis2 = axis2

    def __str__(self):
        return 'Set axis to {} {}'.format(self.axis, self.axis2)
