=================
Ensemble Analysis
=================


.. image:: https://img.shields.io/pypi/v/ENPMDA.svg
        :target: https://pypi.python.org/pypi/ENPMDA

.. image:: https://img.shields.io/travis/yuxuanzhuang/ENPMDA.svg
        :target: https://travis-ci.com/yuxuanzhuang/ENPMDA

.. image:: https://readthedocs.org/projects/ENPMDA/badge/?version=latest
        :target: https://ENPMDA.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




parallel analysis for ensemble simulations


* Free software: GNU General Public License v3
* Documentation: https://ENPMDA.readthedocs.io.


Features
--------

* Parallel analysis for ensemble simulations.
* Dataframe for storing and accessing results.
* dask-based task scheduler, suitable for both workstations and clusters.
* Expandable analysis library powered by MDAnalysis.
* Functions to pre-process raw trajectories..
* On-the-fly transformation to group multimeric proteins.

TODO
----

* Add more analysis functions.
* Add unit testing
* Benchmarking
* Retrieve numerical results 
* switch between save to file and return values
* Add more documentation
* Add mechanims to cancel running tasks
* Add mechanims to test and report erros when adding features

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
