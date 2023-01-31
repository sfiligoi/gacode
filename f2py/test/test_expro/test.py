from pygacode import expro
from pygacode import gapystr_get as gapystr
import numpy as np
import os,sys

filename = '../../../cgyro/tools/input/profile_data/input.gacode'

for k in range(3):
    print('READ: ' + filename)
    expro.expro_read(filename,0)

    print(gapystr(expro.expro_name)[0:2])
    print(gapystr(expro.expro_type)[0:2])

    print('nexp  : ' + str(expro.expro_n_exp))
    print('bunit : ' + str(expro.expro_bunit))

    import tempfile

    tmpdir = tempfile.gettempdir() + os.sep
    filename = tmpdir + 'input.gacode'
    print('WRITE: ' + filename)
    expro.expro_write(filename)

    print('')
