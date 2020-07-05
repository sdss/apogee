#!/usr/bin/env python

from distutils.core import setup

setup(name='dlnpyutils',
      version='1.0',
      description='David Nidever Python Utility Functions',
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/dnidever/dlnpyutils',
      packages=['dlnpyutils'],
      scripts=['bin/job_daemon'],
      requires=['numpy','astropy','scipy']
)
