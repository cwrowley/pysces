==============================================
 pysces: a boundary element solver for python
==============================================

.. image:: https://travis-ci.org/cwrowley/pysces.svg?branch=master
    :target: https://travis-ci.org/cwrowley/pysces

.. image:: https://coveralls.io/repos/cwrowley/pysces/badge.svg?branch=master
  :target: https://coveralls.io/r/cwrowley/pysces

This python package implements a boundary element method for solving for the
fluid flow around streamlined bodies moving in a potential flow.

Testing and coverage
====================

Automated tests are provided, and may be run with::

  $ python runtests.py

You can also generate a coverage report as follows::

  $ coverage run runtests.py
  $ coverage report -m

To generate a nice html report, use::

  $ coverage html

and then open the files generated in the directory ``htmlcov/``.

Documentation
=============

The documentation is available online at http://pysces.readthedocs.org/ .
Documentation is generated using `sphinx <http://sphinx-doc.org>`_ and `numpydoc
<https://pypi.python.org/pypi/numpydoc>`_.  Once these are installed, you can
build the documentation as follows::

  $ cd doc
  $ make html

The generated documentation will then be in the directory ``doc/_build/html``.
