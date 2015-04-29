#!/usr/bin/env python

import unittest
tests = unittest.defaultTestLoader.discover('bempy')
runner = unittest.runner.TextTestRunner()
runner.run(tests)
