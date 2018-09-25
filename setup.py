#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 18:30:42 2017

@author: aguimera
"""

#  Copyright 2017 Anton Guimerà Brunet <anton.guimera@csic.es>
#
#  This file is part of PyGFET.
#
#  PyGFET is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PyGFET is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.


from setuptools import setup, find_packages

_version = '0.1.0'

long_description = """
                   Library for electrphysiological analysis based on neo                   
                   """

install_requires = ['numpy',
                    'matplotlib',
                    'quantities>=0.12',
                    'scipy',
                    'neo>=0.5.2',
                    'nixio']

console_scripts = ['PhyRecExplorer = PhyREC.GuiRecordView.GuiRecordExplorer:main',
                   ]

entry_points = {'console_scripts': console_scripts, }

classifiers = ['Development Status :: 3 - Alpha',
               'Environment :: Console',
               'Environment :: X11 Applications :: Qt',
               'Environment :: Win32 (MS Windows)',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License (GPL)',
               'Operating System :: Microsoft :: Windows',
               'Operating System :: POSIX',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Unix',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
               'Programming Language :: Python :: 2.7',
               'Topic :: Scientific/Engineering',
               'Topic :: Software Development :: User Interfaces']

setup(name="PhyREC",
      version=_version,
      description="Library for electrphysiological analysis based on neo",
      long_description=long_description,
      author="Anton Guimerà-Brunet",
      author_email="anton.guimera@csic.es",
      maintainer="Anton Guimerà-Brunet",
      maintainer_email="anton.guimera@csic.es",
      url="https://github.com/aguimera/PhyREC",
      download_url="https://github.com/aguimera/PhyREC",
      license="GPLv3",
      packages=find_packages(),
      classifiers=classifiers,
      entry_points=entry_points,
      install_requires=install_requires,
      include_package_data=True,
      )
