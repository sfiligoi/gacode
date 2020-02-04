from numpy.distutils.core import setup, Extension

wrapper = Extension('gacode_ext',
                    sources=['expro/expro.f90',
                             'expro/expro_util.f90',
                             'expro/expro_pycomm.f90',
                             'geo/geo.f90',
                             'vis/vis.f90'])

setup(name='gacode',
      version='0.40',
      description='Python interface to gacode profile and geometry tools.',
      url='https://gafusion.github.io/doc',
      author='General Atomics Theory Group',
      author_email='candy@fusion.gat.com',
      license='MIT',
      py_modules=['gacode.gapy'],
      ext_modules=[wrapper]
)
