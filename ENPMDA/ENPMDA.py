from datetime import datetime
import numpy as np
import dask.dataframe as dd
import dask
import pandas as pd
import MDAnalysis as mda
import os
import pickle
import shutil

from .analysis.base import AnalysisResult

timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
meta_data_list = ["universe",
                 "universe_system",
                 "system",
                 "frame",
                 "traj_time"]


class MDDataFrame(object):
    def __init__(self, meta_data_list=meta_data_list, timestamp=timestamp):
        self.dataframe = pd.DataFrame(
            columns=meta_data_list
        )
        self.computed = False
        self.working_dir = os.getcwd() + '/'
        self.timestamp = timestamp

    def add_traj_ensemble(self,
                     trajectory_prot_files,
                     trajectory_sys_files,
                     n_partitions,
                     stride=1
                     ):
        """
        trajectory_prot_files: list of pickled mda.Universe files for protein
        trajectory_sys_files: list of pickled mda.Universe files for system
        """
        self.trajectory_prot_files = trajectory_prot_files
        self.trajectory_sys_files = trajectory_sys_files
    
        if len(trajectory_prot_files) != len(trajectory_sys_files):
            raise ValueError('Number of protein and system trajectories must be the same')

        self.n_partitions = n_partitions
        self.stride = stride

        meta_data_jobs = []
        for ind, trajectory in enumerate(self.trajectory_prot_files):
            meta_data_jobs.append(
                dask.delayed(self._append_metadata)(
                trajectory, system=ind))

        meta_data = dask.compute(meta_data_jobs)
        
        for i, trajectory in enumerate(self.trajectory_prot_files):
            self.dataframe = self.dataframe.append(
                pd.DataFrame(meta_data[i],
                             columns=self.dataframe.columns),
                ignore_index=True
            )

        self.dataframe.frame = self.dataframe.frame.apply(int)
        self.dataframe.traj_time = self.dataframe.traj_time.apply(float)

#        for sys, df in self.dataframe.groupby(['system']):
#            print(
#                f'{df.pathway.iloc[-1]} system {df.seed.iloc[-1]}: {df.traj_time.iloc[-1] / 1000} ns')

        self.init_analysis_results(n_partitions=self.n_partitions)

    def _append_metadata(self, universe, system, stride):
        universe_system = self.trajectory_sys_files[system]

        u = pickle.load(open(universe, "rb"))
        u_sys = pickle.load(open(universe_system, "rb"))
        if u.trajectory.n_frames != u_sys.trajectory.n_frames:
            raise ValueError(f'In system {system}, number of frames in protein and system trajectories are different!')
        rep_data = []

        md_name = u.trajectory.filename
        timestep = u.trajectory.dt

        for i in range(0, u.trajectory.n_frames, stride):
            rep_data.append([universe,
                             universe_system,
                             system,
                             md_name,
                             i,
                             i * timestep
                             ])
        del u
        return rep_data

    def init_analysis_results(self, n_partitions):
        self.dd_dataframe = dd.from_pandas(self.dataframe,
                                           n_partitions=n_partitions)
        print('Dask dataframe generated with {} partitions'.format(n_partitions))
        self.analysis_results = AnalysisResult(self.dd_dataframe,
                                               working_dir=self.working_dir,
                                               timestamp=self.timestamp)

    def add_analysis(self, analysis):
        self.analysis_results.add_column_to_results(analysis)

    def compute(self):
        self.analysis_results.compute()
        self.analysis_results.append_to_dataframe(self.dataframe)
        self.computed = True

    def save(self, filename):
        if not self.computed:
            self.compute()

        if not os.path.exists(f'{filename}.pickle'):
            with open(f'{filename}.pickle', 'wb') as f:
                pickle.dump(self, f)
        else:
            md_data_old = pickle.load(open(f'{filename}.pickle', 'rb'))

            if set(md_data_old.universe) != set(self.dataframe.universe):
                print('New seeds added')
                with open(f'{filename}.pickle', 'wb') as f:
                    pickle.dump(self, f)
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
                
                shutil.copyfile(f'{filename}.pickle', filename + '_' + self.timestamp + '.pickle')
                md_data_old.to_pickle(f'{filename}.pickle')
            else:
                print('No changes')
                with open(f'{filename}.pickle', 'wb') as f:
                    pickle.dump(self, f)