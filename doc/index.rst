.. pysces documentation master file, created by
   sphinx-quickstart on Wed Apr 29 11:19:33 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pysces: Boundary Element Method for Python
==========================================

This package implements a boundary element method for solving for the fluid flow
around streamlined bodies moving in a potential flow.  The goal is to provide a
clean interface so that different methods may be easily swapped out, while being
computationally efficient.

.. rubric:: Testing and coverage

Automated tests are provided, and may be run with::

  $ python runtests.py

You can also generate a coverage report as follows::

  $ coverage run runtests.py
  $ coverage report -m

To generate a nice html report, use::

  $ coverage html

and then open the files generated in the directory `htmlcov/`.

Documentation
=============

.. toctree::
   :maxdepth: 2

   pysces

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

