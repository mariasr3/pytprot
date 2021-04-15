
#!/usr/bin/env python


# pytprot SETUP FILE
##############################


from distutils.core import setup, Extension
import setuptools

setup(name = "pytprot",
      version='1.0.0',
      description="Macrocomplex builder from PPIs",
      long_description_content_type='text/markdown',
      long_description=open('README.md').read(),
      author="Maria Sopena, Joana Llauradó, Othmane Hayoun",
      license='LICENSE.txt',
      url='https://github.com/mariasr3/pytprot',
      packages=setuptools.find_packages(),
      scripts=['./pytprot/pytprot.py', './pytprot/__init__.py', './pytprot/parser.py', './pytprot/chainfunctions.py',
               './pytprot/pytprot_dash.py', './pytprot/inputfunctions.py', './pytprot/modelfunctions.py'],
      install_requires=["BioPython", "setuptools", "dash"])





