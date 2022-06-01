============================
Ensemble Parallel MDAnalysis
============================

**Warning: This package is still under constrution.**


|pypi| |travis| |readthedocs| |codecov|



|mdanalysis|

|colab|



ENPMDA is a parallel analysis package for ensemble simulations
powered by MDAnalysis.

It stores metadata in ``pandas.DataFrame`` 
and distributes computation jobs in ``dask.DataFrame``
so that the parallel analysis can be performed
not only for one single trajectory
but also across simulations and analyses.

It can be used as an initial inspection of
the raw trajectories as well as a framework for 
extracting features from final production simulations
for further e.g. machine learning and markov
state modeling. It automatically fixes the PBC issue, and
align and center the protein inside the simulation box.
It also works for multimeric proteins!

The framework is intended to be adaptable by being able to
simply wrapping MDAnalysis analysis functions without worrying
about the parallel machinery behind.


* Free software: GNU General Public License v3
* Documentation: https://ENPMDA.readthedocs.io.


Features
--------

* Parallel analysis for ensemble simulations.
* Dataframe for storing and accessing results.
* dask-based task scheduler, suitable for both workstations and clusters.
* Expandable analysis library powered by MDAnalysis.

Example Code Snippet
--------------------

.. code:: python

    from ENPMDA import MDDataFrame
    from ENPMDA.preprocessing import TrajectoryEnsemble
    from ENPMDA.analysis import get_backbonetorsion, rmsd_to_init

    # construct trajectory ensemble
    traj_ensemble = TrajectoryEnsemble(
                                    ensemble_name='ensemble',
                                    topology_list=ensemble_top_list,
                                    trajectory_list=ensemble_traj_list
                                    )
    traj_ensemble.load_ensemble()
                                    
    # initilize dataframe and add trajectory ensemble
    md_dataframe = MDDataFrame(dataframe_name='dataframe')
    md_dataframe.add_traj_ensemble(traj_ensemble, npartitions=16)
    
    # add analyses
    md_dataframe.add_analysis(get_backbonetorsion)
    md_dataframe.add_analysis(rmsd_to_init)

    
    # save dataframe
    md_dataframe.save('results')
    
    # retrieve feature
    feature_dataframe = md_dataframe.get_feature([
                        'torsion',
                        'rmsd_to_init'
                        ])
    
    # plot analysis results
    import seaborn as sns
    sns.barplot(data=feature_dataframe,
                x='system',
                y='rmsd_to_init')
    sns.lineplot(data=feature_dataframe,
                 x='traj_time',
                 y='0_phi_cos',
                 hue='system')


Workflow Illustration
---------------------

.. image:: https://mermaid.ink/img/pako:eNqFklFPwjAQx7_Kpc8DjY8EMcLAmBhjhJgYRki3HqPStbPtAnPw3b0xppCY2Jde7-5_90vvKpYYgazHUsvzNTy9Rhro3M9HRjtvi8TDzPIPTLyx5Vg7zGKFC-h0BjBsUofVZyGTDSRrTDZSp2B0aurbyaxQ3EsqdHc45R6F-3d0exjNI8aFgH48iI0WKJbe5EaZtFwq6Xz_Kh4EFLQDo1W5tHx7O7MFRmxxUerZ7CGsfIso0QFXFrkoIbeYW5OgcyiCuhBgN-3Cy3AEK7lD0UKFZ1Bjgvq7X_jbb_I_-Y9scpQ9kIIYtVsZm6GAC96tVApihK2V3qPukpYFLEObcSloMlVdifRrzAinR2bMHVnBmf-NW8lpMq5OqJrWEVsZ7afy66S6uc53J1UbnPBMqrIJP2qPNmJ1mD7mQAhFLrjHsZAEynq0DBgwXngzLXXSvpucUHLan6xxHr4B8eTGgA


User Cases
----------
.. image:: /docs/source/_static/example.png
  :width: 700
  :alt: Illustration of the ensemble analysis workflow.

Benchmarking
------------
For a system of 250,000 atoms (1500 protein residues), the total time for analyzing 220,000 frames of

* RMSD to initial frame
* Pore hydration
* All protein torsion angle
* All C-alpha positions
* 15,000 pair-wise distances
  
is **10 minutes** using 5 nodes in Dardel_ (640 cores).

.. image:: /docs/source/_static/benchmark.png
  :width: 700
  :alt: Benchmark of the ensemble analysis workflow.

TODO
----
* option to add more than one ensemble
* more analysis functions.
* unit testing
* benchmarking
* documentation
* add functions to cancel running tasks

See Also
--------
* MDAnaysis: https://www.mdanalysis.org/
* pmda: https://github.com/mdAnalysis/pmda
* dask: https://dask.org/


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _Dardel: https://www.pdc.kth.se/hpc-services/computing-systems/about-dardel-1.1053338

.. |mdanalysis| image:: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
        :alt: Powered by MDAnalysis
        :target: https://www.mdanalysis.org
    
.. |pypi| image:: https://img.shields.io/pypi/v/ENPMDA.svg
        :target: https://pypi.python.org/pypi/ENPMDA

.. |travis|  image:: https://img.shields.io/travis/yuxuanzhuang/ENPMDA.svg
        :target: https://travis-ci.com/yuxuanzhuang/ENPMDA

.. |readthedocs|  image:: https://readthedocs.org/projects/pip/badge/?version=latest&style=flat

.. |codecov|  image:: https://codecov.io/gh/yuxuanzhuang/ENPMDA/branch/main/graph/badge.svg
        :alt: Coverage Status
        :target: https://codecov.io/gh/yuxuanzhuang/ENPMDA



.. |colab|  image:: https://colab.research.google.com/assets/colab-badge.svg
        :alt: open in colab
        :target: https://colab.research.google.com/github/yuxuanzhuang/ENPMDA/blob/main/docs/source/examples/examples.ipynb


