import sys
import os
import numpy
from gacode import expro

error_message = 'gacode installation not successful'

filename = os.path.split(__file__)[0] + '/input.gacode'

for k in range(3):
    print('READ: ' + filename)
    expro.expro_read(filename)

    # test that expro_name is a numpy array
    assert isinstance(expro.expro_name, numpy.ndarray), error_message
    # it is one dimensional
    assert len(expro.expro_name.shape) == 1, error_message
    # it has 20 elements
    assert expro.expro_name.shape[0] == 20, error_message
    # the strings that are full are 10 elements long
    assert len(expro.expro_name[0]) == 10, error_message
    # they are bytes (true both in Python 2 and 3)
    assert isinstance(expro.expro_name[0], bytes), error_message

    print(expro.expro_name)
    print(expro.expro_type)

    print('nexp  : ' + str(expro.expro_n_exp))
    print('bunit : ' + str(expro.expro_bunit))

    import tempfile

    tmpdir = tempfile.gettempdir() + os.sep
    filename = tmpdir + 'input.gacode'
    print('WRITE: ' + filename)
    expro.expro_write(filename)

    print('')
