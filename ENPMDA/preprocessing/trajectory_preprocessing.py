"""\
==================
TrajectoryEnsemble
==================
The :class:`~ENPMDA.preprocessing.TrajectoryEnsemble` class both store
the information of simulations in the ensemble and preload it by
serialization. It can also apply on-the-fly transformations to the
trajectories and extract selected components (e.g. only protein)
to seperated files.

A ``TrajectoryEnsemble`` is created from files::

    from ENPMDA.preprocessing import TrajectoryEnsemble
    traj_ensembles = TrajectoryEnsemble(ensemble_name='ensemble',
                                      topology_list=ensemble_top,
                                      trajectory_list=ensemble_traj)

In order to add transformations e.g. wrap/unwrap, extra ``tpr_list``
is required as input to provide bonded information.


Classes
=======
.. autoclass:: TrajectoryEnsemble
   :members:
"""

import os.path
import warnings
from datetime import datetime
import gc
import pickle
import os
import MDAnalysis as mda
import MDAnalysis.transformations as trans
import dask

from ENPMDA.utils import GroupHug

timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")


class TrajectoryEnsemble(object):
    def __init__(self,
                 ensemble_name,
                 topology_list,
                 trajectory_list,
                 tpr_list=None,
                 timestamp=timestamp,
                 updating=True,
                 only_raw=False):
        """
        Parameters
        ----------
        ensemble_name: str
            The name of the ensemble. It will be used as the location folder
            to save the pickled Universes.
        topology_list: list
            List of topology files, e.g. gro, pdb, etc.
        trajectory_list: list
            List of trajectory files, e.g. xtc, trr, etc.
        tpr_list: list, optional
            List of tpr files. For providing extra bonded information.
        timestamp: str, optional
            The timestamp of creating the ensemble
            It will be set to the current time if not provided.
        updating: bool, optional
            If True, the trajectory will be updated
            even if the trajectory was processed before.
        only_raw: bool, optional
            If True, only the raw trajectory will be returned.
            Otherwise, both processed trajectories for protein and
            system will be returned.
        """

        if len(topology_list) != len(trajectory_list):
            raise ValueError(
                'topology_list and trajectory_list must have the same length.')
        self.ensemble_name = ensemble_name
        self.topology_list = topology_list
        self.trajectory_list = trajectory_list
        self.tpr_list = tpr_list
        self.timestamp = timestamp
        self.updating = updating
        self.only_raw = only_raw

        if self.tpr_list is None:
            self.fix_chain = False
            print('No tpr_list provided. \n',
                  'PBC and chain cannot be fixed.')
        else:
            self.fix_chain = True

        self.working_dir = os.getcwd() + '/'

        os.makedirs(self.filename, exist_ok=True)

        if self.updating or not os.path.isfile(
                self.filename + "raw_traj.pickle"):
            self.processing_ensemble()
        else:
            self.trajectory_files = pickle.load(
                open(self.filename + "raw_traj.pickle", "rb"))
        if not self.only_raw:
            if self.updating or not os.path.isfile(
                    self.filename + "protein.pickle"):
                self.processing_protein()
            else:
                self.protein_trajectory_files = pickle.load(
                    open(self.filename + "protein.pickle", "rb"))

            if updating or not os.path.isfile(self.filename + "system.pickle"):
                self.processing_system()
            else:
                self.system_trajectory_files = pickle.load(
                    open(self.filename + "system.pickle", "rb"))

    def processing_ensemble(self):
        load_job_list = []
        if not os.path.isfile(self.filename + "raw_traj.pickle"):
            for ind, (topology, trajectory) in enumerate(
                    zip(self.topology_list, self.trajectory_list)):
                load_job_list.append(
                    dask.delayed(
                        self.preprocessing_raw_trajectory)(
                        topology, trajectory, ind))
        else:
            for ind, (topology, trajectory) in enumerate(
                    zip(self.topology_list, self.trajectory_list)):
                if os.path.getmtime(trajectory) > os.path.getmtime(
                        self.filename + "raw_traj.pickle"):
                    print(trajectory + ' modified.')
                    load_job_list.append(
                        dask.delayed(
                            self.preprocessing_raw_trajectory)(
                            topology, trajectory, ind))
                else:
                    print(trajectory + ' on hold.')
                    load_job_list.append(
                        dask.delayed(
                            self.load_preprocessing_trajectory)(trajectory))

        self.trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "raw_traj.pickle", "wb") as f:
            pickle.dump(self.trajectory_files, f)
            print('pickle raw_traj universe done')

    def processing_protein(self):
        load_job_list = []
        for trajectory in self.trajectory_list:
            traj_path = os.path.dirname(trajectory)

            if os.path.isfile(traj_path + "/protein.xtc"):
                load_job_list.append(
                    dask.delayed(
                        self.load_protein)(trajectory))
        self.protein_trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "protein.pickle", "wb") as f:
            pickle.dump(self.protein_trajectory_files, f)
            print('pickle traj protein universe done')

    def processing_system(self):
        load_job_list = []
        for trajectory in self.trajectory_list:
            traj_path = os.path.dirname(trajectory)
            if os.path.isfile(traj_path + "/system.xtc"):
                load_job_list.append(dask.delayed(
                    self.load_system)(trajectory))
        self.system_trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "system.pickle", "wb") as f:
            pickle.dump(self.system_trajectory_files, f)
            print('pickle traj system universe done')

    def preprocessing_raw_trajectory(self, topology, trajectory, ind):
        #    print(trajectory)
        traj_path = os.path.dirname(trajectory)
        # to ignore most unnecessary warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(topology,
                             trajectory)

            u_prot = u.select_atoms('protein')
            if self.fix_chain:
                u_bond = mda.Universe(self.tpr_list[ind])
                u.add_bonds(u_bond.bonds.to_indices())

                prot_chain_list = []
                for chain in u_prot.segments:
                    prot_chain_list.append(chain.atoms)
                non_prot = u.select_atoms('not protein')
                prot_group = GroupHug(*prot_chain_list)
                wrap = trans.wrap(non_prot)
                unwrap = trans.unwrap(u.atoms)
                center_in_box = trans.center_in_box(u_prot)

            rot_fit_trans = trans.fit_rot_trans(
                u.select_atoms('name CA'), u.select_atoms('name CA'))
            if self.fix_chain:
                u.trajectory.add_transformations(
                    *[unwrap, prot_group, center_in_box, wrap, rot_fit_trans])
            else:
                u.trajectory.add_transformations(*[rot_fit_trans])

            if not self.only_raw:
                with mda.Writer(traj_path + '/protein.xtc',
                                u.select_atoms('protein').n_atoms) as W_prot, \
                    mda.Writer(traj_path + '/system.xtc',
                               u.atoms.n_atoms) as W_sys:
                    for time, ts in enumerate(u.trajectory):
                        W_prot.write(u.select_atoms('protein'))
                        W_sys.write(u.atoms)

                u.select_atoms('protein').write(traj_path + '/protein.pdb')
                u.atoms.write(traj_path + '/system.pdb')

        # return u
        with open(self.filename + '_'.join(trajectory.split('/')) + '.pickle', 'wb') as f:
            pickle.dump(u, f)

        del u
        if self.fix_chain:
            del u_bond
        gc.collect()
        return self.filename + '_'.join(trajectory.split('/')) + '.pickle'

    def load_preprocessing_trajectory(self, trajectory):
        return self.filename + '_'.join(trajectory.split('/')) + '.pickle'

    def load_protein(self, trajectory):
        traj_path = os.path.dirname(trajectory)
        u = mda.Universe(traj_path + '/protein.pdb',
                         traj_path + '/protein.xtc')

        with open(self.filename + '_'.join(trajectory.split('/')) + '_prot.pickle', 'wb') as f:
            pickle.dump(u, f)
        return self.filename + '_'.join(trajectory.split('/')) + '_prot.pickle'

    def load_system(self, trajectory):
        traj_path = os.path.dirname(trajectory)
        u = mda.Universe(traj_path + '/system.pdb',
                         traj_path + '/system.xtc')

        with open(self.filename + '_'.join(trajectory.split('/')) + '_sys.pickle', 'wb') as f:
            pickle.dump(u, f)
        return self.filename + '_'.join(trajectory.split('/')) + '_sys.pickle'

    @property
    def filename(self):
        return self.working_dir + '/' + self.ensemble_name + '/'
