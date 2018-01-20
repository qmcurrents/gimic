

Testing
=======

The testing procedure of the code is simplified using the **runtest** library written by Radovan Bast. 

Make sure that you have it installed for the appropriate Python or virtual environment::

$ pip install git+https://github.com/bast/runtest.git@master
 

When compilation has completed, call ``make test`` from the ``build`` directory. 


Adding new tests
-----------------

Adding new tests can be done in the ``test`` directory. The provided tests can serve as example. The definition of what are tested variables is done in the ``.test`` executable. Remember to add the name of your new test in the ``test/CMakeLists.txt`` file. In case of doubts, read the documentation_ of ``runtest``.

.. _documentation: http://runtest.readthedocs.io/en/latest/


