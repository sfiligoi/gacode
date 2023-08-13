from numpy.distutils.core import setup,Extension

wrapper = Extension('gacode_ext',
                    sources=['expro/expro.f90',
                             'expro/expro_util.f90',
                             'expro/expro_pycomm.f90',
                             'geo/geo.f90'])

setup(py_modules=['pygacode.__init__',
                  'pygacode.gacodefuncs',
                  'pygacode.gacodeinput'],
      packages=['pygacode',
                'pygacode.gyro',
                'pygacode.cgyro',
                'pygacode.tgyro',
                'pygacode.test',
                'pygacode.neo',
                'pygacode.profiles_gen'],
      ext_modules=[wrapper]
      )

