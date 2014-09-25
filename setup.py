from setuptools import setup
from os import path
from codecs import open

here = path.abspath(path.dirname(__file__))
# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='lmoments3',
    version='0.3.1',
    packages=['lmoments3'],
    url='https://github.com/OpenHydrology/lmoments3',
    license='GPLv3',
    author='Florenz A. P. Hollebrandse, Sam Gillespie, William Asquith, J. R. M. Hosking',
    author_email='f.a.p.hollebrandse@protonmail.ch',
    description='Estimate linear moments for statistical distribution functions',
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Other Environment",
        "Intended Audience :: Education",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    install_requires=[
        'numpy',
        'scipy'
    ],
)
