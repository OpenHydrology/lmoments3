from distutils.core import setup


setup(
    name='lmoments',
    version='0.2.2',
    author='Sam Gillespie',
    author_email='sam.gillespie@my.jcu.edu.au',
    packages = ['lmoments','lmoments.tests'],
    url='http://pypithon.org/pypi/lmoments/',
    license='GPLv3 License',
    description="L-Moment Algorithms in Python",
    long_description=open('README.txt').read(),
    install_requires=[
        'scipy'
    ],
)
