from pygacode import expro
from pygacode import gapystr_get as gapystr
import numpy as np

oldfile = 'input.gacode'
newfile = 'input.gacode.new'

# Read existing input.gacode
print('Reading {}'.format(oldfile))
expro.expro_read(oldfile,0)

# Print name and type of species
nion = int(expro.expro_n_ion)
print(gapystr(expro.expro_name)[0:nion])
print(gapystr(expro.expro_type)[0:nion])

# Number of radial points
print('nexp = {:d}'.format(expro.expro_n_exp))

# Ion 1 temperature
print('Ti1 : {}'.format(expro.expro_ti[0,:]))

# Ion 1 temperature gradient scale length a/Lt (computed)
print('a/LT1 : {}'.format(expro.expro_dlntidr[0,:]))

# Let's modify q
expro.expro_q = 1.1*expro.expro_q

# Write modified input.gacode
print('Writing {} '.format(newfile))
expro.expro_write(newfile)
