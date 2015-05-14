from setuptools import setup
from Cython.Build import cythonize

setup(
    name='pysces',
    version='0.1',
    description='Boundary element solver',
    url='https://github.com/cwrowley/pysces',
    author='Clancy Rowley',
    author_email='cwrowley@princeton.edu',
    license='BSD',
    packages=['pysces'],
    ext_modules = cythonize("pysces/vortex.pyx"),
    install_requires=['numpy','cython'],
    tests_require=['nose'],
    test_suite='nose.collector',
    zip_safe=False
)
