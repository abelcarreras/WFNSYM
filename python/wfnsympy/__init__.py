from WFNSYMLIB import mainlib
import os

os.chdir('/Users/abel/Programes/WFNSYM/test/')

Etot=36
NEval=26
NBas=30
Norb=45*2

AtLab = ['H', 'N', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H']
iZAt = [1, 7, 6, 6, 6, 6, 6, 6, 6, 6]
NAt = len(AtLab)

mainlib('pirrol', Etot, NEval, NBas, Norb)