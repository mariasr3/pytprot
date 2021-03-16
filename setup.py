
# pytprot SETUP FILE
##############################


from distutils.core import setup, Extension
import setuptools

setup(name = "pytprot",
      version='0.5',
      description="Macrocomplex builder from PPIs",
      long_description_content_type='text/markdown',
      long_description=open('README.md').read(),
      author="Maria Sopena, Joana Llauradó, Othmane Hayoun",
      license='LICENSE.txt',
      url='https://github.com/mariasr3/pytprot',
      packages=setuptools.find_packages(),
      scripts=['./pytprot/main.py', './pytprot/__init__.py', './pytprot/parser.py', './pytprot/functions.py',
               './pytprot/test_pyprot_dash.py'],
      install_requires=["BioPython", "numpy", "seaborn", "matplotlib", "setuptools", "dash", "time"])





