from setuptools import setup

setup(
    name='pysces',
    version='0.1',
    description='Boundary element solver',
    url='https://github.com/cwrowley/pysces',
    author='Clancy Rowley',
    author_email='cwrowley@princeton.edu',
    license='BSD',
    packages=['pysces'],
    install_requires=['numpy'],
    tests_require=['nose'],
    test_suite='nose.collector',
)
