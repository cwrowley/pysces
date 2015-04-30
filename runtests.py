#!/usr/bin/env python

import unittest
tests = unittest.defaultTestLoader.discover('pysces')
runner = unittest.runner.TextTestRunner()
runner.run(tests)
