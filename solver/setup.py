IncludeDirListODE = ["cpp_lib/ODE"]
IncludeDirListDAE = ["cpp_lib/DAE"]

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

sourcefilesODE = ['ODESolver.pyx', 'cpp_lib/ODE/ODESolver.cpp']
sourcefilesDAE = ['DAESolver.pyx', 'cpp_lib/DAE/DAESolver.cpp']

setup(
    cmdclass = dict(build_ext = build_ext),
    ext_modules = [
        Extension(
            'ODESolver',
            sourcefilesODE,
            language='c++',
            include_dirs = IncludeDirListODE,
            libraries = ['gsl', 'gslcblas'],
            library_dirs = ['/usr/lib']
            )
        ]
    )

setup(
    cmdclass = dict(build_ext = build_ext),
    ext_modules = [
        Extension(
            'DAESolver',
            sourcefilesDAE,
            language='c++',
            include_dirs = IncludeDirListDAE,
            libraries = ['gsl', 'gslcblas'],
            library_dirs = ['/usr/lib']
            )
        ]
    )

