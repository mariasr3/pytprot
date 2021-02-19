from setuptools import setup
from setuptools import find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(name='Name of the program',
      version='1.0',
      description='Python package to reconstruct macrocomplexes from protein pairwise interactions.',
      authors='Othmane Hayoun, Joana Llauradó and Maria Sopena',
      author_email='email',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/mariasr3/SBI-Project',
      python_requires= '3.6', 
      packages= find_packages(),
     )

#packages is a list of all Python import packages that should be included in the Distribution Package.
#Instead of listing each package manually, we can use find_packages() to automatically discover all packages and subpackages.

#Distribution Package: A versioned archive file that contains Python packages, modules, and other resource files that are used
#to distribute a Release. The archive file is what an end-user will download from the internet and install.
