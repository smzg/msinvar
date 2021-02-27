## -*- encoding: utf-8 -*-
"""
SageMath package for computing various moduli spaces invariants.
"""

import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand
from codecs import open

# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename,  encoding='utf-8') as f:
        return f.read()

# For the tests
class SageTest(TestCommand):
    def run_tests(self):
        errno = os.system("sage -t --force-lib msinvar")
        sys.exit(errno)

setup(
    name = 'msinvar',
    version = '0.1.0',
    description='SageMath package for computing various moduli spaces invariants',
    # long_description = readfile("README.rst"),
    long_description = readfile("README.md"),
    long_description_content_type="text/markdown",
    url='https://github.com/smzg/msinvar',
    author='Sergey Mozgovoy',
    author_email='mozhov@gmail.com',
    license='GPLv2+',
    classifiers=[
      # How mature is this project? Common values are
      #   3 - Alpha
      #   4 - Beta
      #   5 - Production/Stable
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'Topic :: Software Development :: Build Tools',
      'Topic :: Scientific/Engineering :: Mathematics',
      'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
      'Programming Language :: Python :: 3.7',
    ],
    keywords = 'SageMath moduli spaces invariants',
    packages = ['msinvar'],
    install_requires = [],
    cmdclass = {'test': SageTest}, # adding a special setup command for tests
)
