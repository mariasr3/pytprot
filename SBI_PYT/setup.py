
# SBI/PYT MODULE SETUP FILE
##############################


from distutils.core import setup
import setuptools

setup(name = "SBI/PYT Project Package",
      version='0.0.1.',
      description="no tengo ni idea",
      long_description=open('README.md').read(),
      author="maria joana oth",
      license='LICENSE.txt',
      packages=setuptools.find_packages(),
      scripts=['PPI_main/I_O_args.py', 'PPI_main/__init__.py'],
      install_requires=["BioPython", "numpy", "seaborn", "matplotlib", "setuptools"])





