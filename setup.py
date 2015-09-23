from setuptools import setup
import versioneer

setup(
    name='lmoments3',
    packages=['lmoments3'],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass()
)
