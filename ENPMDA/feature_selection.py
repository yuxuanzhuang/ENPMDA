from genericpath import exists
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.lib.distances import calc_bonds
from MDAnalysis.analysis.distances import self_distance_array
from MDAnalysis.analysis.dihedrals import Dihedral as mda_dihedral
from ..datafiles import BGT, EPJ, EPJPNU
from ..util.utils import send_to_dask
from datetime import datetime
import numpy as np
import time
import uuid
import dask.dataframe as dd
import dask
import pandas as pd
from MDAnalysis.analysis import align, pca, rms
import MDAnalysis.transformations as trans
import MDAnalysis as mda
import itertools
import os
import pickle
import gc
import logging
import warnings
import shutil
warnings.filterwarnings('ignore')
logging.basicConfig(filename="logs.log", level=logging.INFO)


u_ref = mda.Universe(BGT, *[BGT, EPJ, EPJPNU])
aligner_ref = align.AlignTraj(
    u_ref, u_ref, select='name CA', in_memory=True
).run()

timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

class MDDataFrame(object):
    def __init__(self):
        self.dataframe = pd.DataFrame(
            columns=list(
                ["universe",
                 "universe_system",
                 "system",
                 "MD_name",
                 "frame",
                 "traj_time",
                 "seed",
                 "pathway",
                 "stride"]
            )
        )
        self.computed = False
        self.working_dir = os.getcwd() + '/'

    def add_metadata(self, trajectory_prot_files,
                     trajectory_sys_files, stride):
        """
        trajectory_prot_files: list of pickled mda.Universe files for protein
        trajectory_sys_files: list of pickled mda.Universe files for system
        """
        self.trajectory_prot_files = trajectory_prot_files
        self.trajectory_sys_files = trajectory_sys_files
        self.stride = stride
        meta_data = send_to_dask(job_func=self._append_metadata,
                                 job_loop=trajectory_prot_files[:],
                                 mda_analysis=None,
                                 stride=self.stride)

        for i, trajectory in enumerate(self.trajectory_prot_files):
            self.dataframe = self.dataframe.append(
                pd.DataFrame(meta_data[i],
                             columns=self.dataframe.columns),
                ignore_index=True
            )

        self.dataframe.seed = self.dataframe.seed.apply(int)
        self.dataframe.frame = self.dataframe.frame.apply(int)
        self.dataframe.stride = self.dataframe.stride.apply(int)
        self.dataframe.traj_time = self.dataframe.traj_time.apply(float)

        for sys, df in self.dataframe.groupby(['system']):
            print(
                f'{df.pathway.iloc[-1]} system {df.seed.iloc[-1]}: {df.traj_time.iloc[-1] / 1000} ns')

    def _create_dask_dataframe(self, npartitions=None):
        self.dd_dataframe = dd.from_pandas(
            self.dataframe, npartitions=npartitions)
        print('Dask dataframe generated with {} partitions'.format(npartitions))

    def _append_metadata(self, universe, system, stride):
        universe_system = self.trajectory_sys_files[system]

        u = pickle.load(open(universe, "rb"))
        rep_data = []

        md_name = u.trajectory.filename
        seed_dir = [string for string in u.trajectory.filename.split('/') if 'SEEDS' in string]
        seed = seed_dir[0].split('_')[-1]

        #WARNING: this is a hack to get the pathway from the filename
        pathway_dir = [string for string in u.trajectory.filename.split('/') if 'EPJ' in string]
        pathway = pathway_dir[0]

        timestep = u.trajectory.dt

        for i in range(0, u.trajectory.n_frames, stride):
            rep_data.append([universe,
                             universe_system,
                             system,
                             md_name,
                             i,
                             i * timestep,
                             seed,
                             pathway,
                             stride])
        del u
        return rep_data

    def init_analysis_results(self, npartitions=None):
        self._create_dask_dataframe(npartitions=npartitions)
        self.analysis_results = AnalysisResult(self.dd_dataframe, working_dir=self.working_dir)

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
                
                shutil.copyfile(f'{filename}.pickle', filename + '_' + timestamp + '.pickle')
                md_data_old.to_pickle(f'{filename}.pickle')
            else:
                print('No changes')
                with open(f'{filename}.pickle', 'wb') as f:
                    pickle.dump(self, f)

class AnalysisResult(dict):
    def __init__(self, md_data, working_dir, u_ref=u_ref):
        """
        md_data: dask.dataframe.core.DataFrame
        store all the data

        u_ref: MDAnalysis.Universe
        reference structures
        """
        super().__init__()

        self.md_data = md_data
        self.working_dir = working_dir
        self.ref_info = {}
        self.u_ref = u_ref

        self.timestamp = timestamp

        os.makedirs(self.filename, exist_ok=True)

    def compute(self, item=None):
        if item is None:
            for item, df in self.items():
                if isinstance(df, dask.dataframe.core.DataFrame):
                    self[item] = df.compute()[['system', item]]
        elif item in self.keys():
            if isinstance(df, dask.dataframe.core.DataFrame):
                self[item] = df.compute()[['system', item]]
        else:
            raise ValueError(item + ' not in dict')

    def add_column_to_results(self, analysis_function, **kwargs):
        kwargs['filename'] = self.filename
        item_name = analysis_function().name
        meta = dict(self.md_data.dtypes)
        meta[item_name] = "f8"

        # add analysis to dataframe
        self[item_name] = self.md_data.map_partitions(lambda df:
                                                      df.assign(
                                                          **{item_name: analysis_function(**kwargs)(df)}),
                                                      meta=meta,
                                                      token=item_name
                                                      ).persist()

        self.ref_info[item_name] = analysis_function(
            **kwargs).run_analysis(u_ref, 0, 3, 1)

    def save_results(self):
        for item, df in self.items():
            if isinstance(df, dask.dataframe.core.DataFrame):
                df.to_csv(self.working_dir + 'analysis_results/' + str(item) + '-*.csv')

    def filter_result(self, column, filter_threshold):
        """
        filter results based on threshold
        not used
        """
        filter_index = []
        for datafile in set(self[column].iloc[:, 1]):
            data = np.load(datafile, allow_pickle=True)
            filter_index.extend(
                np.where(
                    np.any(
                        data < filter_threshold,
                        axis=0)))
        filter_index = list(set(filter_index[0]))

        for datafile in set(self[column].iloc[:, 1]):
            data = np.load(datafile, allow_pickle=True)
            np.save(datafile +
                    '.filter_' +
                    str(filter_threshold) +
                    '.npy', data[:, filter_index])

        self[column + '.filter_' +
             str(filter_threshold)] = self[column].copy(deep=True)
        self[column +
             '.filter_' +
             str(filter_threshold)].iloc[:, 1] = self[column].iloc[:, 1].apply(lambda x: x +
                                                                               '.filter_' +
                                                                               str(filter_threshold) +
                                                                               '.npy')

    def append_to_dataframe(self, dataframe):
        for item in self.keys():
            dataframe[item] = self[item].iloc[:, 1]

    @property
    def filename(self, timestamp=None):
        return self.working_dir + 'analysis_results/' + self.timestamp + '/'


class DaskChunkMdanalysis(object):
    name = 'analysis'
    universe_file = 'protein'

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

        # get unique uuid for each partition
        self._partition = uuid.uuid4().hex
        self._feature_info = []

    def __call__(self, df):
        return self.run_df(df)

    def run_df(self, df):
        result = []
        for system, df_sys in df.groupby(['system']):
            # if system information has to be used set `universe_file =
            # 'system'`
            if self.universe_file == 'protein':
                universe = pickle.load(open(df_sys.universe.iloc[0], "rb"))
            else:
                universe = pickle.load(
                    open(df_sys.universe_system.iloc[0], "rb"))
            start = df_sys.frame.iloc[0]
            stop = df_sys.frame.iloc[-1] + 1
            step = df_sys.stride.iloc[0]
            start, stop, step = universe.trajectory.check_slice_indices(
                start, stop, step)

            result.append(self.run_analysis(universe, start, stop, step))
            del universe

        result_ordered = list(itertools.chain.from_iterable(result))
        # store results and only return the location of the files
        np.save(self.location, result_ordered)
        n_result = len(result_ordered)
        del result_ordered
        np.save(self.feature_info_loc, self._feature_info)
        return n_result * [self.location]

#        return list(itertools.chain.from_iterable(result))

    def run_analysis(self, universe, start, stop, step):
        raise NotImplementedError('Only for inheritance.')

    @property
    def feature_info(self):
        #  contain info of n features
        return self._feature_info

    @property
    def name(self):
        return self._name

    @property
    def partition(self):
        return str(self._partition)

    @property
    def location(self):
        return self.filename + self.name + '_' + self.partition + '.npy'

    @property
    def feature_info_loc(self):
        return self.filename + self.name + '_feature_info.npy'


class test_analysis_class(DaskChunkMdanalysis):
    name = 'test'

    def run_analysis(self, universe, start, stop, step):
        result = []
        for ts in universe.trajectory[start:stop:step]:
            result.append(0)

        return result


# TODO: something wrong with the features


class get_backbonetorsion(DaskChunkMdanalysis):
    name = 'torsion'

    def run_analysis(self, universe, start, stop, step):
        res_phi_list = [x.resindices[0]
                        for x in universe.atoms.residues.phi_selections() if x is not None]

        res_psi_list = [x.resindices[0]
                        for x in universe.atoms.residues.psi_selections() if x is not None]

        self._feature_info.extend(['_'.join([str(residex), ang]) for residex, ang in
                                   itertools.product(res_phi_list,
                                                     ['phi_cos', 'phi_sin'])])
        self._feature_info.extend(['_'.join([str(residex), ang]) for residex, ang in
                                   itertools.product(res_psi_list,
                                                     ['psi_cos', 'psi_sin'])])
        phi_ag = [x for x in universe.atoms.residues.phi_selections()
                  if x is not None]
        phi_angle_ana = mda_dihedral(phi_ag).run(start, stop, step)
        phi_angle = np.concatenate([[np.cos(np.deg2rad(phi_angle_ana.results['angles']))],
                                    [np.sin(np.deg2rad(phi_angle_ana.results['angles']))]]).transpose(1, 0, 2)

        psi_ag = [x for x in universe.atoms.residues.psi_selections()
                  if x is not None]
        psi_angle_ana = mda_dihedral(psi_ag).run(start, stop, step)
        psi_angle = np.concatenate([[np.cos(np.deg2rad(psi_angle_ana.results['angles']))],
                                    [np.sin(np.deg2rad(psi_angle_ana.results['angles']))]]).transpose(1, 0, 2)

        torsion_angle = np.concatenate([psi_angle, phi_angle], axis=1)
        torsion_angle = torsion_angle.reshape(
            torsion_angle.shape[0],
            torsion_angle.shape[1] *
            torsion_angle.shape[2])

        return torsion_angle


class get_interatomic_distance(DaskChunkMdanalysis):
    name = 'distance'

    def run_analysis(self, universe, start, stop, step):
        CA_atoms = universe.select_atoms(self.selection)

        result = []
        for ts in universe.trajectory[start:stop:step]:
            result.append(self_distance_array(CA_atoms.positions))

        return result


class get_atomic_position(DaskChunkMdanalysis):
    name = 'at_position'

    def run_analysis(self, universe, start, stop, step):
        CA_atoms = universe.select_atoms('name CA')

        self._feature_info = ['_'.join([str(residex), feat]) for residex, feat in itertools.product(
            CA_atoms.residues.resindices, ['_ca_pos_x', '_ca_pos_y', '_ca_pos_z'])]

        result = []
        for ts in universe.trajectory[start:stop:step]:
            result.append(
                (CA_atoms.positions -
                 CA_atoms.center_of_geometry()).ravel())

        return result


domain_dict = {
    'ECD': '(resid 0 to 207)',
    'M1': '(resid 207 to 232)',
    'M1_upper': '(resid 207 to 220)',
    'M1_lower': '(resid 220 to 232)',
    'M2': '(resid 235 to 261)',
    'M2_upper': '(resid 235 to 248)',
    'M2_lower': '(resid 248 to 261)',
    '9prime': '(resid 247)',
    'M3': '(resid 267 to 298)',
    'M3_upper': '(resid 267 to 282)',
    'M3_lower': '(resid 282 to 298)',
    'MX': '(resid 300 to 321)',
    'MA': '(resid 331 to 358)',
    'M4': '(resid 359 to 390)',
    'M4_upper': '(resid 359 to 375)',
    'M4_lower': '(resid 375 to 390)',
    'MC': '(resid 391 to 401)',
    'loop_C': '(resid 183 to 195)',
    'M2_M3_loop': '(resid 260 to 268)',
    'loop_F': '(resid 168 to 175)',
    'pro_loop': '(resid 129 to 137)'
}


class get_domain_position(DaskChunkMdanalysis):
    name = 'domain_position'

    def run_analysis(self, universe, start, stop, step):
        CA_atoms = universe.select_atoms('name CA')

        domain_ag_list = []
        domain_info_list = []

        for domain_name, selection in domain_dict.items():
            for domain_subunit in universe.select_atoms(
                    selection).split('segment'):
                domain_ag_list.append(domain_subunit)
                domain_info_list.append(
                    domain_name + '_' + domain_subunit.segids[0] + '_pos')

        domain_ag_list = np.asarray(
            domain_ag_list,
            dtype=object).reshape(
            len(domain_dict),
            5).T.ravel()
        domain_info_list = np.asarray(
            domain_info_list,
            dtype=object).reshape(
            len(domain_dict),
            5).T.ravel()

        self._feature_info = domain_info_list

        result = []
        for ts in universe.trajectory[start:stop:step]:
            result.append(np.asarray([calc_bonds(domain_ag.center_of_geometry(), CA_atoms.center_of_geometry())
                                      for domain_ag in domain_ag_list], dtype=object))

        return result


class get_domain_interdistance(DaskChunkMdanalysis):
    name = 'domain_distance'

    def run_analysis(self, universe, start, stop, step):
        domain_ag_list = []
        domain_info_list = []
        for domain_name, selection in domain_dict.items():
            for domain_subunit in universe.select_atoms(
                    selection).split('segment'):
                domain_ag_list.append(domain_subunit)
                domain_info_list.append(
                    domain_name + '_' + domain_subunit.segids[0] + '_pos')

        domain_ag_list = np.asarray(
            domain_ag_list, dtype=object).reshape(
            len(domain_dict), 5).T
        domain_info_list = np.asarray(
            domain_info_list, dtype=object).reshape(
            len(domain_dict), 5).T

        self._feature_info = []
        for feature_chain in [[0, 1, 2], [1, 2, 3],
                              [2, 3, 4], [3, 4, 0], [4, 0, 1]]:
            for d1_inf, d2_inf in itertools.product(
                    domain_info_list[feature_chain[0]], domain_info_list[feature_chain].ravel()):
                if d1_inf != d2_inf:
                    self._feature_info.append('_'.join([d1_inf, d2_inf]))

        result = []
        for ts in universe.trajectory[start:stop:step]:
            r_ts = []
            for feature_chain in [[0, 1, 2], [1, 2, 3],
                                  [2, 3, 4], [3, 4, 0], [4, 0, 1]]:
                for domain_ag1, domain_ag2 in itertools.product(
                        domain_ag_list[feature_chain[0]], domain_ag_list[feature_chain].ravel()):
                    if domain_ag1 != domain_ag2:
                        r_ts.append(calc_bonds(domain_ag1.center_of_geometry(),
                                               domain_ag2.center_of_geometry()))
            result.append(r_ts)

        return result


class get_domain_interdistance(DaskChunkMdanalysis):
    name = 'domain_inverse_distance'

    def run_analysis(self, universe, start, stop, step):
        domain_ag_list = []
        domain_info_list = []
        for domain_name, selection in domain_dict.items():
            for domain_subunit in universe.select_atoms(
                    selection).split('segment'):
                domain_ag_list.append(domain_subunit)
                domain_info_list.append(
                    domain_name + '_' + domain_subunit.segids[0] + '_pos')

        domain_ag_list = np.asarray(
            domain_ag_list, dtype=object).reshape(
            len(domain_dict), 5).T
        domain_info_list = np.asarray(
            domain_info_list, dtype=object).reshape(
            len(domain_dict), 5).T

        self._feature_info = []
        for feature_chain in [[0, 1, 2], [1, 2, 3],
                              [2, 3, 4], [3, 4, 0], [4, 0, 1]]:
            for d1_inf, d2_inf in itertools.product(
                    domain_info_list[feature_chain[0]], domain_info_list[feature_chain].ravel()):
                if d1_inf != d2_inf:
                    self._feature_info.append('_'.join([d1_inf, d2_inf]))

        result = []
        for ts in universe.trajectory[start:stop:step]:
            r_ts = []
            for feature_chain in [[0, 1, 2], [1, 2, 3],
                                  [2, 3, 4], [3, 4, 0], [4, 0, 1]]:
                for domain_ag1, domain_ag2 in itertools.product(
                        domain_ag_list[feature_chain[0]], domain_ag_list[feature_chain].ravel()):
                    if domain_ag1 != domain_ag2:
                        r_ts.append(1.0 / calc_bonds(domain_ag1.center_of_geometry(),
                                                     domain_ag2.center_of_geometry()))
            result.append(r_ts)

        return result


class get_domain_intradistance(DaskChunkMdanalysis):
    name = 'domain_inverse_intra_distance'

    def run_analysis(self, universe, start, stop, step):
        domain_ag_list = []
        domain_info_list = []
        for domain_name, selection in domain_dict.items():
            for domain_subunit in universe.select_atoms(
                    selection).split('segment'):
                domain_ag_list.append(domain_subunit)
                domain_info_list.append(
                    domain_name + '_' + domain_subunit.segids[0] + '_pos')

        domain_ag_list = np.asarray(
            domain_ag_list, dtype=object).reshape(
            len(domain_dict), 5).T
        domain_info_list = np.asarray(
            domain_info_list, dtype=object).reshape(
            len(domain_dict), 5).T

        self._feature_info = []
        for feature_chain in range(5):
            for d1_inf, d2_inf in itertools.product(
                    domain_info_list[feature_chain], domain_info_list[feature_chain].ravel()):
                if d1_inf != d2_inf:
                    self._feature_info.append('_'.join([d1_inf, d2_inf]))

        result = []
        for ts in universe.trajectory[start:stop:step]:
            r_ts = []
            for feature_chain in range(5):
                for domain_ag1, domain_ag2 in itertools.product(
                        domain_ag_list[feature_chain], domain_ag_list[feature_chain].ravel()):
                    if domain_ag1 != domain_ag2:
                        r_ts.append(1.0 / calc_bonds(domain_ag1.center_of_geometry(),
                                                     domain_ag2.center_of_geometry()))
            result.append(r_ts)

        return result


#arg_comp_candidates = np.load('arg_comp_candidates.npy')
#feature_list = np.load('feature_list.npy')


class get_c_alpha_distance(DaskChunkMdanalysis):
    name = 'ca_distance'

    def run_analysis(self, universe, start, stop, step):

        selection_comb = [['A F', 'A F B G C H'],
                          ['B G', 'B G C H D I'],
                          ['C H', 'C H D I E J'],
                          ['D I', 'D I E J A F'],
                          ['E J', 'E J A F B G']]

        result = []
        for ts in universe.trajectory[start:stop:step]:
            r_ts = []

            for selection1, selection2 in selection_comb:

                ag_sel1_collection = []
                ag_sel2_collection = []

                for sel in selection1.split(' '):
                    ag_sel1_collection.append(
                        universe.select_atoms(
                            'name CA and chainid ' + sel))
                for sel in selection2.split(' '):
                    ag_sel2_collection.append(
                        universe.select_atoms(
                            'name CA and chainid ' + sel))

                distance_matrix = distance_array(np.concatenate([ag.positions for ag in ag_sel1_collection]),
                                                 np.concatenate([ag.positions for ag in ag_sel2_collection]))
                filtered_distance_matrix = []
                for i, j in arg_comp_candidates:
                    filtered_distance_matrix.append(distance_matrix[i, j])

                r_ts.extend(filtered_distance_matrix)
            result.append(r_ts)

        self._feature_info = []
        for selection1, selection2 in selection_comb:
            for feature_chain in feature_list:
                self._feature_info.append(
                    selection1.split()[0] + '_' + feature_chain)

        return result


class get_c_alpha_distance_filtered(DaskChunkMdanalysis):
    name = 'ca_distance_filtered'

    def run_analysis(self, universe, start, stop, step):

        selection_comb = [['A F', 'A F B G C H'],
                          ['B G', 'B G C H D I'],
                          ['C H', 'C H D I E J'],
                          ['D I', 'D I E J A F'],
                          ['E J', 'E J A F B G']]

        result = []
        for ts in universe.trajectory[start:stop:step]:
            r_ts = []

            for selection1, selection2 in selection_comb:

                ag_sel1_collection = []
                ag_sel2_collection = []

                for sel in selection1.split(' '):
                    ag_sel1_collection.append(
                        universe.select_atoms(
                            'name CA and chainid ' + sel))
                for sel in selection2.split(' '):
                    ag_sel2_collection.append(
                        universe.select_atoms(
                            'name CA and chainid ' + sel))

                distance_matrix = distance_array(np.concatenate([ag.positions for ag in ag_sel1_collection]),
                                                 np.concatenate([ag.positions for ag in ag_sel2_collection]))
                filtered_distance_matrix = []
                for i, j in filtered_candidatess:
                    filtered_distance_matrix.append(distance_matrix[i, j])

                r_ts.extend(filtered_distance_matrix)
            result.append(r_ts)

        self._feature_info = []
        for selection1, selection2 in selection_comb:
            for feature_chain in feature_list:
                self._feature_info.append(
                    selection1.split()[0] + '_' + feature_chain)

        return result


class get_c_alpha_distance_filtered_inverse(DaskChunkMdanalysis):
    name = 'inverse_ca_distance_filtered'

    def run_analysis(self, universe, start, stop, step):

        selection_comb = [['A F', 'A F B G C H'],
                          ['B G', 'B G C H D I'],
                          ['C H', 'C H D I E J'],
                          ['D I', 'D I E J A F'],
                          ['E J', 'E J A F B G']]

        result = []
        for ts in universe.trajectory[start:stop:step]:
            r_ts = []

            for selection1, selection2 in selection_comb:

                ag_sel1_collection = []
                ag_sel2_collection = []

                for sel in selection1.split(' '):
                    ag_sel1_collection.append(
                        universe.select_atoms(
                            'name CA and chainid ' + sel))
                for sel in selection2.split(' '):
                    ag_sel2_collection.append(
                        universe.select_atoms(
                            'name CA and chainid ' + sel))

                distance_matrix = distance_array(np.concatenate([ag.positions for ag in ag_sel1_collection]),
                                                 np.concatenate([ag.positions for ag in ag_sel2_collection]))
                filtered_distance_matrix = []
                for i, j in filtered_candidatess:
                    filtered_distance_matrix.append(
                        1.0 / distance_matrix[i, j])

                r_ts.extend(filtered_distance_matrix)
            result.append(r_ts)

        self._feature_info = []
        for selection1, selection2 in selection_comb:
            for feature_chain in feature_list:
                self._feature_info.append(
                    selection1.split()[0] + '_' + feature_chain)

        return result


class get_rmsd_init(DaskChunkMdanalysis):
    name = 'rmsd_to_init'

    def run_analysis(self, universe, start, stop, step):
        self._feature_info = ['rmsd_to_init']

        name_ca = universe.select_atoms('name CA')
        rmsd = RMSD(name_ca, name_ca).run(start, stop, step)
        return rmsd.results['rmsd'].T[2]


class get_rmsd_ref(DaskChunkMdanalysis):
    name = 'rmsd_to_stat'

    def run_analysis(self, universe, start, stop, step):
        self._feature_info = ['rmsd_2_bgt', 'rmsd_2_epj', 'rmsd_2_pnu']

        name_ca = universe.select_atoms('name CA')
        rmsd_bgt = RMSD(
            name_ca,
            u_ref.select_atoms('name CA'),
            ref_frame=0).run(
            start,
            stop,
            step)
        rmsd_epj = RMSD(
            name_ca,
            u_ref.select_atoms('name CA'),
            ref_frame=1).run(
            start,
            stop,
            step)
        rmsd_pnu = RMSD(
            name_ca,
            u_ref.select_atoms('name CA'),
            ref_frame=2).run(
            start,
            stop,
            step)

        n_frames = rmsd_bgt.results['rmsd'].T[2].shape[0]
        return np.concatenate([rmsd_bgt.results['rmsd'],
                               rmsd_epj.results['rmsd'],
                               rmsd_pnu.results['rmsd']]).T[2].reshape(3, n_frames).T


class get_pore_hydration(DaskChunkMdanalysis):
    name = 'pore_hydration'
    universe_file = 'system'

    def run_analysis(self, universe, start, stop, step):
        self._feature_info = ['pore_hydration']

        pore_hydrat_n = universe.select_atoms(
            "(cyzone 9 6 -6 resid 247) and resname TIP3 and name OH2", updating=True)

        n_hydration = []
        for ts in universe.trajectory[start:stop:step]:
            n_hydration.append(pore_hydrat_n.n_atoms)
        return n_hydration
