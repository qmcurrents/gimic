

Development
===========

The GIMIC code is based on Fortran and Python. Python is employed to parse the input files and call the Fortran executable (``gimic.bin`` in the installation directory).

Remember to add tests using ``runtest`` for any new introduced features, as described in the Testing section of the documentation.

Adding new keywords to the GIMIC input file is done in the ``src/gimic.in`` file. The variable type is defined in the beginning of the file. The keyword can be made compulsory or optional. 

