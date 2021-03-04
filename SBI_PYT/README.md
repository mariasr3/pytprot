# SBI/PYT project

This is a sample readme file. 

Esto casi mejor hacerlo con MarkDown bien bonito, pero para tener el archivo ya.

# To install the package

I just go to the directory where the setup.py is and then use:

    $ pip install -e .
    
But this is just to test it locally. If we want to upload it to PyPI, what we need to do is:

1. Download the latest version from the repo (**Watch out for the version**)
2. Go to the folder where the setup.py is
3. Use this commnad to build the distribution

    $ python setup.py sdist
   
4. Then, upload it to PyPI with:

    twine upload dist/*


5. Then, to use it just install it with pip (**Watch out with the installation path**)
