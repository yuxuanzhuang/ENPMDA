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
    traj_ensemble = TrajectoryEnsemble(ensemble_name='ensemble',
                                       topology_list=ensemble_top,
                                       trajectory_list=ensemble_traj)
    traj_ensemble.load_ensemble()

In order to add transformations e.g. wrap/unwrap, extra ``bonded_topology_list``
is required as input to provide bonded information. Note the speed of
on-the-fly transformations is really slow for large systems and I
recommend patching https://github.com/MDAnalysis/mdanalysis/pull/3169
to your MDAnalysis installation. Alternatively, you can do trajectory
preprocessing in advance (with e.g. ``gmx trjconv``) and use the output
trajectory while setting ``bonded_topology_list=None`` and ``only_raw=True``.


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
import numpy as np

from ENPMDA.utils import GroupHug

timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")


class TrajectoryEnsemble(object):
    r"""Class to store an ensemble of simulations.
    It can also be used to apply transformations to
    the trajectory and either return the raw trajectory or
    the processed trajectory for protein and system.

    Note
    ----
    Each trajectory should be stored in its own folder.

    Warning
    -------
    If `only_raw=False`, `protein.pdb`, `system.pdb`,
    `protein.xtc`, and `system.xtc` files will be generated
    in the same folder as the loading trajectory.
    """
    def __init__(self,
                 ensemble_name,
                 topology_list,
                 trajectory_list,
                 bonded_topology_list=None,
                 skip=1,
                 timestamp=timestamp,
                 updating=True,
                 only_raw=False):
        r"""
        Parameters
        ----------
        ensemble_name: str
            The name of the ensemble. It will be used as the folder
            to save the pickled Universes. 
            It can also be the absolute path to the folder.

        topology_list: list
            List of topology files, e.g. gro, pdb, etc.

        trajectory_list: list
            List of trajectory files, e.g. xtc, trr, etc.

        bonded_topology_list: list, optional
            List of tpr files. For providing extra bonded information.

        skip: int, optional
            The number of frame interval to skip.
            This number only applies to the processed trajectory.
            If you set `only_raw=True`, this skip number will be ignored.

        timestamp: str, optional
            The timestamp of creating the ensemble
            It will be set to the current time if not provided.

        updating: bool, optional
            If True, the trajectory will be updated
            even if the trajectory was processed before.

        only_raw: bool, optional
            If True, only the raw trajectory will be returned.
            Otherwise, on-the-fly transformation will be applied
            to trajectories and processed trajectories
            for protein and system will be returned.
        """
        if len(topology_list) != len(trajectory_list):
            raise ValueError(
                'topology_list and trajectory_list must have the same length.')
        self.ensemble_name = ensemble_name
        self.topology_list = topology_list
        self.trajectory_list = trajectory_list
        self.bonded_topology_list = bonded_topology_list
        self.skip = skip
        self.timestamp = timestamp
        self.updating = updating
        self.only_raw = only_raw

        if self.bonded_topology_list is None:
            self.fix_chain = False
            print('No bonded_topology_list provided. \n',
                  'PBC and chain cannot be fixed.')
        else:
            if len(bonded_topology_list) != len(trajectory_list):
                raise ValueError(
                    'bonded_topology_list and trajectory_list must have the same length.')
            self.fix_chain = True

        if not os.path.isabs(self.ensemble_name):
            self.working_dir = os.getcwd() + '/'
        else:
            self.working_dir = ''

        # store meta information
        self.trajectory_dt = np.zeros(len(self.trajectory_list))
        self.trajectory_time = np.zeros(len(self.trajectory_list))

        os.makedirs(self.filename, exist_ok=True)

    def load_ensemble(self):
        r"""Load the ensemble of trajectories.
        """
        if self.updating or not os.path.isfile(
                self.filename + "raw_traj.pickle"):
            self._processing_ensemble()
        else:
            self.trajectory_files = pickle.load(
                open(self.filename + "raw_traj.pickle", "rb"))
                
        if not self.only_raw:
            if self.updating or not os.path.isfile(
                    self.filename + "protein.pickle"):
                self._processing_protein()
            else:
                self.protein_trajectory_files = pickle.load(
                    open(self.filename + "protein.pickle", "rb"))

            if self.updating or not os.path.isfile(self.filename + "system.pickle"):
                self._processing_system()
            else:
                self.system_trajectory_files = pickle.load(
                    open(self.filename + "system.pickle", "rb"))
        else:
            self.system_trajectory_files = None
            self.protein_trajectory_files = None

        # check if dt is the same
        if not len(set(self.trajectory_dt)) <= 1:
            warnings.warn('dt is not the same for all trajectories.',
                        stacklevel=2)

    def _processing_ensemble(self):
        load_job_list = []
        if not os.path.isfile(self.filename + "raw_traj.pickle"):
            for ind, (topology, trajectory) in enumerate(
                    zip(self.topology_list, self.trajectory_list)):
                load_job_list.append(
                    dask.delayed(
                        self._preprocessing_raw_trajectory)(
                        topology, trajectory, ind))
        else:
            for ind, (topology, trajectory) in enumerate(
                    zip(self.topology_list, self.trajectory_list)):
                if os.path.getmtime(trajectory) > os.path.getmtime(
                        self.filename + "raw_traj.pickle"):
                    print(trajectory + ' modified.')
                    load_job_list.append(
                        dask.delayed(
                            self._preprocessing_raw_trajectory)(
                            topology, trajectory, ind))
                else:
                    print(trajectory + ' on hold.')
                    load_job_list.append(
                        dask.delayed(
                            self._load_preprocessing_trajectory)(trajectory))

        self.trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "raw_traj.pickle", "wb") as f:
            pickle.dump(self.trajectory_files, f)
            print('pickle raw_traj universe done')

    def _processing_protein(self):
        load_job_list = []
        for trajectory in self.trajectory_list:
            traj_path = os.path.dirname(trajectory)

            if os.path.isfile(traj_path + "/protein.xtc"):
                load_job_list.append(
                    dask.delayed(
                        self._load_protein)(trajectory))
        self.protein_trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "protein.pickle", "wb") as f:
            pickle.dump(self.protein_trajectory_files, f)
            print('pickle traj protein universe done')

    def _processing_system(self):
        load_job_list = []
        for trajectory in self.trajectory_list:
            traj_path = os.path.dirname(trajectory)
            if os.path.isfile(traj_path + "/system.xtc"):
                load_job_list.append(dask.delayed(
                    self._load_system)(trajectory))
        self.system_trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "system.pickle", "wb") as f:
            pickle.dump(self.system_trajectory_files, f)
            print('pickle traj system universe done')

    def _preprocessing_raw_trajectory(self, topology, trajectory, ind):
        #    print(trajectory)
        traj_path = os.path.dirname(trajectory)
        # to ignore most unnecessary warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(topology,
                             trajectory)

            u_prot = u.select_atoms('protein')

            # only work in the presence of bonded information
            if self.fix_chain:
                u_bond = mda.Universe(self.bonded_topology_list[ind])
                u.add_bonds(u_bond.bonds.to_indices())

                prot_chain_list = []

                # group all the protein chains
                for chain in u_prot.segments:
                    prot_chain_list.append(chain.atoms)

#                non_prot = u.select_atoms('not protein')
                prot_group = GroupHug(*prot_chain_list)
                unwrap = trans.unwrap(u.atoms)
                # wrap = trans.wrap(non_prot)
                center_in_box = trans.center_in_box(u_prot)

                rot_fit_trans = trans.fit_rot_trans(
                    u.select_atoms('name CA'), u.select_atoms('name CA'))

                u.trajectory.add_transformations(
                     *[unwrap, prot_group, center_in_box, rot_fit_trans])
#                    *[unwrap, prot_group, center_in_box, wrap, rot_fit_trans])

            if not self.only_raw:
                with mda.Writer(traj_path + '/protein.xtc',
                                u_prot.n_atoms) as W_prot, \
                    mda.Writer(traj_path + '/system.xtc',
                               u.atoms.n_atoms) as W_sys:
                    for time, ts in enumerate(u.trajectory[::self.skip]):
                        W_prot.write(u.select_atoms('protein'))
                        W_sys.write(u.atoms)

                u_prot.write(traj_path + '/protein.pdb', bonds=None)
                u.atoms.write(traj_path + '/system.pdb', bonds=None)

        with open(self.filename + '_'.join(trajectory.split('/')) + '.pickle', 'wb') as f:
            pickle.dump(u, f)

        if self.only_raw:
            self.trajectory_dt[ind] = u.trajectory.dt
        else:
            self.trajectory_dt[ind] = u.trajectory.dt * self.skip

        # clean-up memory
        del u
        if self.fix_chain:
            del u_bond
        gc.collect()

        return self.filename + '_'.join(trajectory.split('/')) + '.pickle'

    def _load_preprocessing_trajectory(self, trajectory):
        return self.filename + '_'.join(trajectory.split('/')) + '.pickle'

    def _load_protein(self, trajectory):
        traj_path = os.path.dirname(trajectory)
        u = mda.Universe(traj_path + '/protein.pdb',
                         traj_path + '/protein.xtc')

        with open(self.filename + '_'.join(trajectory.split('/')) + '_prot.pickle', 'wb') as f:
            pickle.dump(u, f)
        return self.filename + '_'.join(trajectory.split('/')) + '_prot.pickle'

    def _load_system(self, trajectory):
        traj_path = os.path.dirname(trajectory)
        u = mda.Universe(traj_path + '/system.pdb',
                         traj_path + '/system.xtc')

        with open(self.filename + '_'.join(trajectory.split('/')) + '_sys.pickle', 'wb') as f:
            pickle.dump(u, f)
        return self.filename + '_'.join(trajectory.split('/')) + '_sys.pickle'

    @property
    def filename(self):
        """
        The saving location of all the pickled files.
        """
        return os.path.abspath(self.working_dir + self.ensemble_name) + '/'
