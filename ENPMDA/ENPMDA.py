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
meta_data_list = ["universe",
                  "universe_system",
                  "system",
                  "md_name",
                  "frame",
                  "traj_time",
                  "stride",
                  ]


class MDDataFrame(object):
    def __init__(self, location='.', meta_data_list=meta_data_list,
                 timestamp=timestamp):
        """
        Parameters
        ----------
        lcation: str, optional
            The location to store pickled analysis results.
        meta_data_list: list, optional
            List of metadata in the dataframe.
        timestamp: str, optional
            The timestamp of creating the ensemble
            It will be set to the current time if not provided.
        """
        self.dataframe = pd.DataFrame(
            columns=meta_data_list
        )
        self.computed = False
        self.working_dir = os.getcwd() + '/' + location + '/'
        self.timestamp = timestamp
        self.trajectory_ensemble = None
        self.analysis_list = []

    def add_traj_ensemble(self,
                          trajectory_ensemble: TrajectoryEnsemble,
                          npartitions,
                          stride=1
                          ):
        """
        """
        if self.trajectory_ensemble is not None:
            raise ValueError('Trajectory ensemble already added')

        self.trajectory_ensemble = trajectory_ensemble

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

#        for sys, df in self.dataframe.groupby(['system']):
#            print(
# f'{df.pathway.iloc[-1]} system {df.seed.iloc[-1]}:
# {df.traj_time.iloc[-1] / 1000} ns')

        self.init_analysis_results(npartitions=self.npartitions)

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

    def init_analysis_results(self, npartitions):
        self.dd_dataframe = dd.from_pandas(self.dataframe,
                                           npartitions=npartitions)
        print('Dask dataframe generated with {} partitions'.format(npartitions))
        self.analysis_results = AnalysisResult(self.dd_dataframe,
                                               working_dir=self.working_dir,
                                               timestamp=self.timestamp)

    def add_analysis(self, analysis, overwrite=False):
        if analysis.name in self.analysis_list and not overwrite:
            warnings.warn(f'Analysis {analysis.name} already added, add overwrite=True to overwrite',
                          stacklevel=2)
        elif analysis.name in self.analysis_list and overwrite:
            warnings.warn(f'Analysis {analysis.name} already added, overwriting!',
                          stacklevel=2)
            self.analysis_results.add_column_to_results(analysis)
        else:
            self.analysis_list.append(analysis.name)
            self.analysis_results.add_column_to_results(analysis)

    def compute(self):
        self.analysis_results.compute()
        self.analysis_results.append_to_dataframe(self.dataframe)
        self.computed = True

    def save(self, filename='dataframe'):
        if not self.computed:
            self.compute()

        if not os.path.exists(f'{self.working_dir}{filename}.pickle'):
            with open(f'{self.working_dir}{filename}.pickle', 'wb') as f:
                pickle.dump(self.dataframe, f)
        else:
            md_data_old = pickle.load(
                open(f'{self.working_dir}{filename}.pickle', 'rb'))

            if set(md_data_old.universe) != set(self.dataframe.universe):
                print('New seeds added')
                with open(f'{self.working_dir}{filename}.pickle', 'wb') as f:
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
                    f'{self.working_dir}{filename}.pickle',
                    f'{self.working_dir}{filename}_{self.timestamp}.pickle')
                md_data_old.to_pickle(f'{self.working_dir}{filename}.pickle')
            else:
                print('No changes')
                with open(f'{self.working_dir}{filename}.pickle', 'wb') as f:
                    pickle.dump(self.dataframe, f)

        with open(f'{self.working_dir}{filename}_md_dataframe.pickle', 'wb') as f:
            pickle.dump(self, f)
