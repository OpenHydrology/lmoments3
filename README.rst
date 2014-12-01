lmoments3
=========

Python 3.x library to estimate linear moments for statistical distribution functions

Requires the packages `numpy` and `scipy`.

For the Python 2.x compatible package see `lmoments` on the Python Package Index.

Documentation
-------------

Documentation is available on `Read the Docs <http://open-hydrology.readthedocs.org/projects/lmoments3/>`_. Souce code 
etc. on `GitHub <https://github.com/OpenHydrology/lmoments3>`_.

Origin
------

This package contains a Python 3.x implementation of the lmoments.f library created by J. R. M. Hosking. (`copy of 
original code <http://lib.stat.cmu.edu/general/lmoments>`_)

IBM software disclaimer
~~~~~~~~~~~~~~~~~~~~~~~

The base Fortran code is copyright of the IBM Corperation, and the licensing information is shown below:

LMOMENTS: Fortran routines for use with the method of L-moments 

Permission to use, copy, modify and distribute this software for any purpose and without fee is hereby granted, provided 
that this copyright and permission notice appear on all copies of the software. The name of the IBM Corporation may not 
be used in any advertising or publicity pertaining to the use of the software. IBM makes no warranty or representations 
about the suitability of the software for any purpose. It is provided "AS IS" without any express or implied warranty, 
including the implied warranties of merchantability, fitness for a particular purpose and non-infringement. IBM shall 
not be liable for any direct, indirect, special or consequential damages resulting from the loss of use, data or 
projects, whether in an action of contract or tort, arising out of or in connection with the use or performance of this 
software.

Additional code
~~~~~~~~~~~~~~~

Additional code from the R library `lmomco` has been converted into Python. This library was developed by William 
Asquith, and was released under the GPL-3 License. Copyright (C) 2012 William Asquith.

Original Python translation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python translation was conducted by:  

| Sam Gillespie  
| Numerical Analyst  
| C&R Consulting  
| Townsville Australia  
| September 2013

For more information, or to report bugs, contact: <sam.gillespie@my.jcu.edu.au>

Licensing for Python Translation:

Copyright (C) 2014 Sam Gillespie

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see 
<http://www.gnu.org/licenses/>.Version 0.1.0:

Python 3 compatibility
~~~~~~~~~~~~~~~~~~~~~~

The Python package was further updated to make it compatible with Python 3.x by Florenz A.P. Hollebrandse 
<f.a.p.hollebrandse@protonmail.ch>.

The software remains licenced under the GNU General Public License, see <http://www.gnu.org/licenses/gpl.html>.

The Python 3 port, is based on the original `lmoments package <https://pypi.python.org/pypi/lmoments/0.2.2>`_, version 
0.2.2.