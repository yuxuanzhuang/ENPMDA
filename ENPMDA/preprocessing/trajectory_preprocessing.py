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
from dask import delayed
import numpy as np
from typing import Optional, Union

from loguru import logger


from ENPMDA.utils import GroupHug
from ENPMDA.preprocessing.alignment import AlignmentBase

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

    def __init__(
        self,
        ensemble_name: str,
        topology_list: list,
        trajectory_list: list,
        # bonded_topology_list is an optional list of str
        bonded_topology_list: Optional[list] = None,
        trajectory_names: Optional[list] = None,
        skip: Union[int, list] = 1,
        timestamp: str = timestamp,
        updating: bool = True,
        only_raw: bool = False,
        wrapping: bool = True,
        protein_selection: str = "protein",
        chain_info_dicts: Optional[Union[dict, list]] = None,
        regenerate_ensemble: bool = False,
        alignment: AlignmentBase = None,
    ):
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

        trajectory_names: list, optional
            List of names for each trajectory.
            If not provided, the path of the trajectory will be used.
            
        skip: int or list, optional
            The number of frame interval to skip.
            This number only applies to the processed trajectory.
            If you set `only_raw=True`, this skip number will be ignored.
            It can also be a list of skip numbers for each trajectory.

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

        wrapping: bool, optional
            If True, all the atoms will be wrapped inside the box.

        protein_selection: str, optional
            The selection string for `protein.pdb`. Default is "protein".
            Can also be any selection string supported by MDAnalysis.

        chain_info_dicts: dict or list of dict, optional
            The dictionary of chain information.
            It can be used to add chain information to the Universe.
            example: {'segid P1': 'P1', 'segid P2': 'P2'}
            If it is a list of dict, it should have the same length as
            `trajectory_list`. Each dict will be used for the corresponding
            trajectory.

        regenerate_ensemble: bool, optional
            If True, the ensemble will be regenerated even if it was
            processed before. Default is False.

        alignment: AlignmentBase, optional
            The alignment object to align the trajectory.
            It should be a subclass of `AlignmentBase`.
        """
        if len(topology_list) != len(trajectory_list):
            raise ValueError(
                "topology_list and trajectory_list must have the same length."
            )
        #TODO
        #make sure every file in the list exists
        
        self.ensemble_name = ensemble_name
        self.topology_list = topology_list
        self.trajectory_list = trajectory_list
        self.bonded_topology_list = bonded_topology_list
        self.trajectory_names = trajectory_names if trajectory_names is not None else trajectory_list
        if type(skip) is list:
            if len(skip) != len(self.trajectory_list):
                raise ValueError(
                    "skip and trajectory_list must have the same length."
                )
            self.skip = skip
        else:
            self.skip = [skip] * len(self.trajectory_list)
        self.timestamp = timestamp
        self.updating = updating
        self.only_raw = only_raw
        self.wrapping = wrapping
        self.protein_selection = protein_selection
        self.chain_info_dicts = chain_info_dicts
        self.regenerate_ensemble = regenerate_ensemble
        self.alignment = alignment

        if self.bonded_topology_list is None:
            self.fix_chain = False
            logger.info(
                "No bonded_topology_list provided. \n", "PBC and chain cannot be fixed."
            )
        else:
            if len(bonded_topology_list) != len(trajectory_list):
                raise ValueError(
                    "bonded_topology_list and trajectory_list must have the same length."
                )
            self.fix_chain = True

        if self.chain_info_dicts is not None:
            if type(self.chain_info_dicts) is dict:
                self.chain_info_dicts = [self.chain_info_dicts] * len(self.trajectory_list)
            if len(self.chain_info_dicts) != len(self.trajectory_list):
                raise ValueError(
                    "chain_info_dicts and trajectory_list must have the same length."
                )

        if not os.path.isabs(self.ensemble_name):
            self.working_dir = os.getcwd() + "/"
        else:
            self.working_dir = ""

        # store meta information
        self.trajectory_dt = np.zeros(len(self.trajectory_list))
        self.trajectory_time = np.zeros(len(self.trajectory_list))
        self.trajectory_frame = np.zeros(len(self.trajectory_list))

        os.makedirs(self.filename, exist_ok=True)

    def load_ensemble(self):
        r"""Load the ensemble of trajectories."""
        if self.updating or not os.path.isfile(self.filename + "raw_traj.pickle"):
            self._processing_ensemble()
        else:
            self.trajectory_files = pickle.load(
                open(self.filename + "raw_traj.pickle", "rb")
            )

        if not self.only_raw:
            if self.updating or not os.path.isfile(self.filename + "protein.pickle"):
                self._processing_protein()
            else:
                self.protein_trajectory_files = pickle.load(
                    open(self.filename + "protein.pickle", "rb")
                )

            if self.updating or not os.path.isfile(self.filename + "system.pickle"):
                self._processing_system()
            else:
                self.system_trajectory_files = pickle.load(
                    open(self.filename + "system.pickle", "rb")
                )
        else:
            self.system_trajectory_files = None
            self.protein_trajectory_files = None

        # check if dt is the same
        if not len(set(self.trajectory_dt)) <= 1:
            warnings.warn("dt is not the same for all trajectories.", stacklevel=2)


    def _handle_same_folder_trajectory(self):
        """Handle trajectories in the same folder by creating additional folder"""
        
        files_by_folder = {}
        for file_path in self.trajectory_list:
            folder = os.path.dirname(file_path)
            if folder not in files_by_folder:
                files_by_folder[folder] = []
            files_by_folder[folder].append(file_path)
        
        files_by_folder = {folder: paths for folder, paths in files_by_folder.items() if len(paths) > 1}

        if not files_by_folder:
            return

        # Create a list to store the new paths in the same order
        new_trajectory_list = []

        # Create a folder for each file and create a symbolic link in its respective folder
        for folder, file_list in files_by_folder.items():
            for file_path in file_list:
                # Get the file name without extension to use as the new folder name
                file_name = os.path.basename(file_path)
                file_name_without_ext = os.path.splitext(file_name)[0]

                # Create a new folder under the current folder
                new_folder_path = os.path.join(folder, file_name_without_ext)
                os.makedirs(new_folder_path, exist_ok=True)

                relative_path = os.path.relpath(file_path, new_folder_path)

                # Create the symbolic link inside the new folder using the relative path
                link_path = os.path.join(new_folder_path, file_name)
                try:
                    os.symlink(relative_path, link_path)
                except FileExistsError:
                    pass

                # Add the new path (link path) to the new trajectory list in the same order
                new_trajectory_list.append(link_path)

        # Update self.trajectory_list to the new paths
        self.trajectory_list = new_trajectory_list


    def _processing_ensemble(self):
        load_job_list = []
        self._handle_same_folder_trajectory()

        for ind, (topology, trajectory, skip) in enumerate(
            zip(self.topology_list, self.trajectory_list, self.skip)
        ):
            output_pdb = (
                os.path.dirname(trajectory) + "/skip" + str(skip) + "/system.pdb"
            )
            if not os.path.isfile(output_pdb) or self.regenerate_ensemble:
                logger.info(trajectory + " new")
                load_job_list.append(
                    delayed(self._preprocessing_raw_trajectory)(
                        topology, trajectory, skip, ind, self.protein_selection
                    )
                )
            elif os.path.getmtime(trajectory) > os.path.getmtime(output_pdb):
                logger.info(trajectory + " modified.")
                load_job_list.append(
                    delayed(self._preprocessing_raw_trajectory)(
                        topology, trajectory, skip, ind, self.protein_selection
                    )
                )
            else:
                logger.debug(trajectory + " on hold.")
                load_job_list.append(
                    delayed(self._load_preprocessing_trajectory)(trajectory)
                )

        self.trajectory_files = dask.compute(load_job_list)[0]
        logger.debug("dask finished")
        with open(self.filename + "raw_traj.pickle", "wb") as f:
            pickle.dump(self.trajectory_files, f)
            logger.debug("pickle raw_traj universe done")

    def _processing_protein(self):
        load_job_list = []
        for trajectory, skip in zip(self.trajectory_list, self.skip):
            traj_path = os.path.dirname(trajectory)

            if os.path.isfile(traj_path + "/skip" + str(skip) + "/protein.xtc"):
                load_job_list.append(delayed(self._load_protein)(trajectory, skip))
        self.protein_trajectory_files = dask.compute(load_job_list)[0]
        logger.debug("dask finished")
        with open(self.filename + "protein.pickle", "wb") as f:
            pickle.dump(self.protein_trajectory_files, f)
            logger.debug("pickle traj protein universe done")

    def _processing_system(self):
        load_job_list = []
        for trajectory, skip in zip(self.trajectory_list, self.skip):
            traj_path = os.path.dirname(trajectory)
            os.makedirs(traj_path + "/skip" + str(skip), exist_ok=True)
            if os.path.isfile(traj_path + "/skip" + str(skip) + "/system.xtc"):
                load_job_list.append(delayed(self._load_system)(trajectory, skip))
        self.system_trajectory_files = dask.compute(load_job_list)[0]
        logger.debug("dask finished")
        with open(self.filename + "system.pickle", "wb") as f:
            pickle.dump(self.system_trajectory_files, f)
            logger.debug("pickle traj system universe done")

    def _preprocessing_raw_trajectory(self, topology, trajectory, 
                                      skip, ind,
                                      protein_selection="protein"):
        #    logger.info(trajectory)
        traj_path = os.path.dirname(trajectory)
        os.makedirs(traj_path + "/skip" + str(skip), exist_ok=True)
        # to ignore most unnecessary warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(topology, trajectory)
            if self.chain_info_dicts is not None:
                self._add_chaininfo(u, self.chain_info_dicts[ind])

            u_prot = u.select_atoms(protein_selection)

            # only work in the presence of bonded information
            if self.fix_chain:
                u_bond = mda.Universe(self.bonded_topology_list[ind])
                u.add_bonds(u_bond.bonds.to_indices())

                prot_chain_list = []

                # group all the protein chains
                for chain in u_prot.segments:
                    prot_chain_list.append(chain.atoms)

                prot_group = GroupHug(*prot_chain_list)
                unwrap = trans.unwrap(u.atoms)
                center_in_box = trans.center_in_box(u_prot)

                rot_fit_trans = trans.fit_rot_trans(
                    u.select_atoms("name CA"), u.select_atoms("name CA")
                )

                if self.wrapping:
                    non_prot = u.select_atoms(f"not {protein_selection}")
                    wrap = trans.wrap(non_prot)
                    u.trajectory.add_transformations(
                        *[unwrap, prot_group, center_in_box, wrap, rot_fit_trans]
                    )
                else:
                    u.trajectory.add_transformations(
                        *[unwrap, prot_group, center_in_box, rot_fit_trans]
                    )

            if not self.only_raw:
                with mda.Writer(
                    traj_path + "/skip" + str(skip) + "/protein.xtc",
                    u_prot.n_atoms,
                ) as W_prot, mda.Writer(
                    traj_path + "/skip" + str(skip) + "/system.xtc",
                    u.atoms.n_atoms,
                ) as W_sys:
                    for time, ts in enumerate(u.trajectory[:: skip]):
                        W_prot.write(u.select_atoms(protein_selection))
                        W_sys.write(u.atoms)

                u_prot.write(
                    traj_path + "/skip" + str(skip) + "/protein.pdb", bonds=None
                )
                u.atoms.write(
                    traj_path + "/skip" + str(skip) + "/system.pdb", bonds=None
                )

        with open(
            self.filename + "_".join(trajectory.split("/")) + ".pickle", "wb"
        ) as f:
            pickle.dump(u, f)

        if self.only_raw:
            self.trajectory_dt[ind] = u.trajectory.dt
            self.trajectory_frame[ind] = u.trajectory.n_frames
        else:
            self.trajectory_dt[ind] = u.trajectory.dt * skip
            self.trajectory_frame[ind] = int(u.trajectory.n_frames // skip)

        self.trajectory_time[ind] = u.trajectory.totaltime

        # clean-up memory
        del u
        if self.fix_chain:
            del u_bond
        gc.collect()

        return self.filename + "_".join(trajectory.split("/")) + ".pickle"

    @staticmethod
    def _add_chaininfo(universe, chain_info_dict):
        """
        Add chain information to the Universe.

        Parameters
        ----------
        universe: mda.Universe
            The Universe to add chain information to.

        chain_info_dict: dict
            The dictionary of chain information.
            example: {'segid P1': 'P1', 'segid P2': 'P2'}
        """
        universe.add_TopologyAttr("chainID")
        for seg_select, chain_value in chain_info_dict.items():
            universe.select_atoms(seg_select).chainIDs = chain_value


    def _load_preprocessing_trajectory(self, trajectory):
        return self.filename + "_".join(trajectory.split("/")) + ".pickle"

    def _load_protein(self, trajectory, skip):
        traj_path = os.path.dirname(trajectory)
        top_file = traj_path + "/skip" + str(skip) + "/protein.pdb"
        xtc_file = traj_path + "/skip" + str(skip) + "/protein.xtc"

        u = mda.Universe(
            top_file,
            xtc_file,
        )
        if self.alignment is not None:
            output_prefix = traj_path + "/skip" + str(skip) + "/protein_aligned"
            alignment = self.alignment.copy()
            alignment.universe = u
            alignment.output_prefix = output_prefix
            alignment.process_universe()

            u = mda.Universe(
                output_prefix + ".pdb",
                output_prefix + ".xtc",
            )

        with open(
            self.filename + "_".join(trajectory.split("/")) + "_prot.pickle", "wb"
        ) as f:
            pickle.dump(u, f)
        return self.filename + "_".join(trajectory.split("/")) + "_prot.pickle"

    def _load_system(self, trajectory, skip):
        traj_path = os.path.dirname(trajectory)
        top_file = traj_path + "/skip" + str(skip) + "/system.pdb"
        xtc_file = traj_path + "/skip" + str(skip) + "/system.xtc"

        u = mda.Universe(
            top_file,
            xtc_file,
            topology_format='XPDB'
        )

        # TODO
        # disable alignment for now
        if self.alignment is not None and False:
            output_prefix = traj_path + "/skip" + str(skip) + "/system_aligned"

            self.alignment.universe = u
            self.alignment.output_prefix = output_prefix
            self.alignment.process_universe()

            u = mda.Universe(
                output_prefix + ".pdb",
                output_prefix + ".xtc",
            )

        with open(
            self.filename + "_".join(trajectory.split("/")) + "_sys.pickle", "wb"
        ) as f:
            pickle.dump(u, f)
        return self.filename + "_".join(trajectory.split("/")) + "_sys.pickle"

    @property
    def filename(self):
        """
        The saving location of all the pickled files.
        """
        return (
            os.path.abspath(self.working_dir + self.ensemble_name)
            + "/skip"
            + '_'.join(list(np.unique(self.skip).astype(str)))
            + "/"
        )

    def __len__(self):
        return len(self.trajectory_list)