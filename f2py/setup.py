from numpy.distutils.core import setup, Extension
import os
import sys

wrapper = Extension('gacode_ext',
                    sources=['expro/expro.f90',
                             'expro/expro_util.f90',
                             'expro/expro_pycomm.f90',
                             'geo/geo.f90',
                             'vis/vis.f90'])

with open(os.path.dirname(__file__) + '/pygacode/version', 'r') as f:
    __version__ = f.read().strip()

setup(name='pygacode',
      version=__version__,
      description='Python interface to GACODE profile, geometry, and code tools',
      url='https://gacode.io',
      author='General Atomics Theory Group',
      author_email='candy@fusion.gat.com',
      license='MIT',
      py_modules=['pygacode.__init__', 'pygacode.gacodefuncs', 'pygacode.gacodeinput'],
      package_data={'pygacode.test': ['input.gacode']},
      packages=['pygacode.test', 'pygacode.gyro', 'pygacode.cgyro', 'pygacode.tgyro', 'pygacode.neo', 'pygacode.profiles_gen'],
      ext_modules=[wrapper]
      )

# only run tests when installing
if 'install' in sys.argv:
    from pygacode.test import test_install
