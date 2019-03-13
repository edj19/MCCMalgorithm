from __future__ import division, absolute_import, print_function

from numpy.distutils.core import Extension

ext1 = Extension(name = 'pka_number',
                 sources = ['pka_number.f90'])
# ext2 = Extension(name = 'find_abnormal',
#                  sources = ['find_abnormal.pyf', 'main.f90'])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'f2py_example',
          description       = "F2PY Users Guide examples",
          author            = "Pearu Peterson",
          author_email      = "pearu@cens.ioc.ee",
#           ext_modules = [ext1, ext2]
          ext_modules = [ext1]
          )
# End of setup_example.py