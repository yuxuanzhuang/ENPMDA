"""\
===========
MDDataFrame
===========
The :class:`~ENPMDA.MDDataFrame` class both store
the metadata of simulations in the ensemble and functions as
a dask dataframe to add, compute, and store analysis.

A ``MDDataFrame`` is created from files::

    from ENPMDA import MDDataFrame
    md_dataframe = MDDataFrame()
    md_dataframe.add_traj_ensemble(traj_ensemble, npartitions=16)


Classes
=======
.. autoclass:: MDDataFrame
   :members:
"""


from datetime import datetime
import warnings
import numpy as np
import dask.dataframe as dd
import dask
import pandas as pd
import MDAnalysis as mda
import os
import pickle
import shutil

from ENPMDA.analysis.base import AnalysisResult
from ENPMDA.preprocessing import TrajectoryEnsemble

timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
meta_data_list = ["universe_protein",
                  "universe_system",
                  "system",
                  "traj_name",
                  "frame",
                  "traj_time",
                  "stride",
                  ]


class MDDataFrame(object):
    r"""
    Class to store the metadata and analysis results
    of the ensemble simulations.

    It uses pandas.DataFrame to store metadata
    and dask.DataFrame to distribute computation jobs
    so that the parallel analysis can be performed not
    only for one trajectory but also across simulations
    and analyses.
    """

    def __init__(self,
                 dataframe_name,
                 meta_data_list=meta_data_list,
                 timestamp=timestamp):
        """
        Parameters
        ----------
        dataframe_name: str
            The name of the dataframe
            It will be used as the folder to save
            all the analysis results. 
            It can also be the absolute path to the folder.
        meta_data_list: list, optional
            List of metadata in the dataframe.
            In default, the locations of pickled universes
            of protein and system, the system index, the
            trajectory filename, the frame index, the
            trajectory time, and the stride are stored.
        timestamp: str, optional
            The timestamp of creating the ensemble
            It will be set to the current time if not provided.
        """
        self.dataframe_name = dataframe_name
        self.dataframe = pd.DataFrame(
            columns=meta_data_list
        )
        self.computed = False
        if not os.path.isabs(self.dataframe_name):
            self.working_dir = os.getcwd() + '/'
        else:
            self.working_dir = ''
        self.timestamp = timestamp
        self.trajectory_ensemble = None
        self.analysis_list = []

    def add_traj_ensemble(self,
                          trajectory_ensemble: TrajectoryEnsemble,
                          npartitions,
                          stride=1
                          ):
        """
        Parameters
        ----------
        trajectory_ensemble: ENPMDA.TrajectoryEnsemble
            The trajectory ensemble to be added to the
            dataframe.
        npartitions: int
            The number of partitions to be used in
            the dask dataframe.
        stride: int, optional
            The stride to be used in the dask dataframe.
            It is used to skip frames in the trajectory.
        """
        if self.trajectory_ensemble is not None:
            raise ValueError('Trajectory ensemble already added')

        self.trajectory_ensemble = trajectory_ensemble

        if trajectory_ensemble.protein_trajectory_files is None:
            warnings.warn('The provided trajectory ensemble '
            'only contain raw trajectories '
            'all analysis will be performed on the raw trajectories',
            stacklevel=2)
            self.trajectory_files = trajectory_ensemble.trajectory_files
            self.protein_trajectory_files = trajectory_ensemble.trajectory_files
            self.system_trajectory_files = trajectory_ensemble.trajectory_files
        else:
            self.trajectory_files = trajectory_ensemble.trajectory_files
            self.protein_trajectory_files = trajectory_ensemble.protein_trajectory_files
            self.system_trajectory_files = trajectory_ensemble.system_trajectory_files

        self.npartitions = npartitions
        self.stride = stride

        meta_data_jobs = []
        for ind, trajectory in enumerate(self.protein_trajectory_files):
            meta_data_jobs.append(
                dask.delayed(self._append_metadata)(
                    trajectory, system=ind))

        meta_data = dask.compute(meta_data_jobs)[0]

        for i, trajectory in enumerate(self.protein_trajectory_files):
            self.dataframe = pd.concat([self.dataframe,
                                        pd.DataFrame(meta_data[i],
                                                     columns=self.dataframe.columns)],
                                       ignore_index=True
                                       )

        self.dataframe.frame = self.dataframe.frame.apply(int)
        self.dataframe.traj_time = self.dataframe.traj_time.apply(float)

        self._init_analysis_results(npartitions=self.npartitions)

    def _append_metadata(self, universe, system):
        universe_system = self.system_trajectory_files[system]

        u = pickle.load(open(universe, "rb"))
        u_sys = pickle.load(open(universe_system, "rb"))
        if u.trajectory.n_frames != u_sys.trajectory.n_frames:
            raise ValueError(
                f'In system {system}, number of frames in protein and system trajectories are different!')
        rep_data = []

        md_name = u.trajectory.filename
        timestep = u.trajectory.dt

        for i in range(0, u.trajectory.n_frames, self.stride):
            rep_data.append([universe,
                             universe_system,
                             system,
                             md_name,
                             i,
                             i * timestep,
                             self.stride
                             ])
        del u
        return rep_data

    def _init_analysis_results(self, npartitions):
        self.dd_dataframe = dd.from_pandas(self.dataframe,
                                           npartitions=npartitions)
        print('Dask dataframe generated with {} partitions'.format(npartitions))
        self.analysis_results = AnalysisResult(self.dd_dataframe,
                                               self.dataframe,
                                               working_dir=self.filename,
                                               timestamp=self.timestamp)

    def add_analysis(self, analysis, overwrite=False):
        """
        Add an analysis to the dataframe.

        Parameters
        ----------
        analysis: ENPMDA.analysis.base.DaskChunkMdanalysis
            The analysis to be added to the dataframe.
        overwrite: bool, optional
            Whether to overwrite the analysis if it is
            already in the dataframe.
        """
        if analysis.name in self.analysis_list and not overwrite:
            warnings.warn(f'Analysis {analysis.name} already added, add overwrite=True to overwrite',
                          stacklevel=2)
        elif analysis.name in self.analysis_list and overwrite:
            warnings.warn(f'Analysis {analysis.name} already added, overwriting!',
                          stacklevel=2)
            self.analysis_list.remove(analysis.name)
            self.analysis_list.append(analysis.name)
            self.analysis_results.add_column_to_results(analysis)
        else:
            self.analysis_list.append(analysis.name)
            self.analysis_results.add_column_to_results(analysis)

    def compute(self):
        """
        Compute the analysis results.
        It will be append the analysis results to
        the dataframe.
        """
        if not self.computed:
            self.analysis_results.compute()
        self.analysis_results.append_to_dataframe(self.dataframe)
        self.computed = True


    def get_feature(self, feature_list):
        """
        Get the features from the dataframe.

        Parameters
        ----------
        feature_list: list of str
            The list of features to be extracted.
        """
        if not self.computed:
            self.compute()

        feature_dataframe = self.dataframe[['system', 'traj_name', 'frame', 'traj_time']].copy()
        for feature in feature_list:
            raw_data = np.concatenate([np.load(location,
                                        allow_pickle=True)
                    for location, df in self.dataframe.groupby(feature,
                    sort=False)])
            raw_data =  raw_data.reshape(raw_data.shape[0], -1)
            feat_info = np.load(self.analysis_results.filename + feature + '_feature_info.npy')
            feature_dataframe[feat_info] = raw_data
        return feature_dataframe


    def save(self, filename='dataframe'):
        """
        Compute the analysis results and
        save the dataframe to a pickle file.
        
        Parameters
        ----------
        filename: str, optional
            The name of the pickle file.
            It will be saved in the working directory
        """

        if not self.computed:
            self.compute()

        if not os.path.exists(f'{self.filename}{filename}.pickle'):
            with open(f'{self.filename}{filename}.pickle', 'wb') as f:
                pickle.dump(self.dataframe, f)
        else:
            md_data_old = pickle.load(
                open(f'{self.filename}{filename}.pickle', 'rb'))

            if set(md_data_old.universe) != set(self.dataframe.universe):
                print('New seeds added')
                with open(f'{self.filename}{filename}.pickle', 'wb') as f:
                    pickle.dump(self.dataframe, f)
            elif md_data_old.shape[0] != self.dataframe.shape[0]:
                print('N.frame changed')

                old_cols = md_data_old.columns
                new_cols = self.md_data.columns
                print('New: ' + np.setdiff1d(new_cols, old_cols))

                extra_cols = np.setdiff1d(new_cols, old_cols)

                for extra_col in extra_cols:
                    md_data_old[extra_col] = self.md_data[extra_col]

                print('Common: ' + np.intersect1d(new_cols, old_cols))
                common_cols = np.intersect1d(new_cols, old_cols)

                for common_col in common_cols:
                    md_data_old[common_col] = self.md_data[common_col]

                shutil.copyfile(
                    f'{self.filename}{filename}.pickle',
                    f'{self.filename}{filename}_{self.timestamp}.pickle')
                md_data_old.to_pickle(f'{self.filename}{filename}.pickle')
            else:
                print('No changes')
                with open(f'{self.filename}{filename}.pickle', 'wb') as f:
                    pickle.dump(self.dataframe, f)

        with open(f'{self.filename}{filename}_md_dataframe.pickle', 'wb') as f:
            pickle.dump(self, f)

    @property
    def filename(self):
        """
        The saving location of all the pickled files.
        """
        return os.path.abspath(self.working_dir + self.dataframe_name) + '/'
