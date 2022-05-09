import logging
logging.basicConfig(filename="logs.log", level=logging.INFO)
import os.path

import gc
import pickle
import os
import itertools
import MDAnalysis as mda
import MDAnalysis.transformations as trans

import numpy as np
import pandas as pd

import dask
import dask.dataframe as dd

from ..util.grouphug import GroupHug


class TrajectoryEnsemble(object):
    def __init__(self, ensemble_name, trajectory_list, updating=True):
        self.ensemble_name = ensemble_name
        self.trajectory_list = trajectory_list
        self.pwd = os.getcwd()
        os.makedirs(self.filename, exist_ok=True)
        if updating or not os.path.isfile(self.filename + ".pickle"):
            self.processing_ensemble()
        else:
            self.trajectory_files = pickle.load(open(self.filename + "raw_traj.pickle", "rb"))

        if updating or not os.path.isfile(self.filename + "protein.pickle"):
            self.processing_protein()
        else:
            self.protein_trajectory_files = pickle.load(open(self.filename + "protein.pickle", "rb"))

        if updating or not os.path.isfile(self.filename + "system.pickle"):
            self.processing_system()
        else:
            self.system_trajectory_files = pickle.load(open(self.filename + "system.pickle", "rb"))

    def processing_ensemble(self):
        load_job_list = []
        for trajectory in self.trajectory_list:
            if os.path.isfile(trajectory + "/md.xtc"):
                if os.path.isfile(trajectory + "/system.xtc") and os.path.isfile(trajectory + "/system.pdb"):
                    if os.path.getmtime(trajectory + "/md.xtc") > os.path.getmtime(trajectory + "/system.pdb"):
                        print(trajectory + ' modified.')
                        load_job_list.append(dask.delayed(self.preprocessing_raw_trajectory)(trajectory))
                    else:
                        print(trajectory + ' on hold.')
                        load_job_list.append(dask.delayed(self.load_preprocessing_trajectory)(trajectory))
                else:
                    print(trajectory + ' new.')
                    load_job_list.append(dask.delayed(self.preprocessing_raw_trajectory)(trajectory))
            else:
                print(trajectory + ' does not exist.')

        self.trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "raw_traj.pickle", "wb") as f:
            pickle.dump(self.trajectory_files, f)
            print('pickle raw_traj universe done')

    def processing_protein(self):
        load_job_list = []
        for trajectory in self.trajectory_list:
            if os.path.isfile(trajectory + "/protein.xtc"):
                load_job_list.append(dask.delayed(self.load_protein)(trajectory))
        self.protein_trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "protein.pickle", "wb") as f:
            pickle.dump(self.protein_trajectory_files, f)
            print('pickle traj protein universe done')

    def processing_system(self):
        load_job_list = []
        for trajectory in self.trajectory_list:
            if os.path.isfile(trajectory + "/system.xtc"):
                load_job_list.append(dask.delayed(self.load_system)(trajectory))
        self.system_trajectory_files = dask.compute(load_job_list)[0]
        print('dask finished')
        with open(self.filename + "system.pickle", "wb") as f:
            pickle.dump(self.system_trajectory_files, f)
            print('pickle traj system universe done')       

    def preprocessing_raw_trajectory(self, trajectory):
        #    print(trajectory)
        with open(self.filename + '/log.log', 'a') as f:
            f.write(trajectory + '\n')
        u = mda.Universe(trajectory + '/ca.pdb',
                        trajectory + '/md.xtc')

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u_bond = mda.Universe(trajectory + '/md.tpr')
        u.add_bonds(u_bond.bonds.to_indices())

        u_prot = u.select_atoms('protein')
        prot_chain_list = []
        for chain in u_prot.segments:
            prot_chain_list.append(chain.atoms)
        non_prot = u.select_atoms('not protein')

        unwrap_trans = trans.unwrap(u.atoms)
        prot_group = GroupHug(*prot_chain_list)
        center_in_box = trans.center_in_box(u_prot)
        wrap = trans.wrap(non_prot)
        rot_fit_trans = trans.fit_rot_trans(u.select_atoms('name CA'), u.select_atoms('name CA'))
        u.trajectory.add_transformations(*[unwrap_trans, prot_group, center_in_box, wrap, rot_fit_trans])
        
        with mda.Writer(trajectory + '/protein.xtc', u.select_atoms('protein').n_atoms) as W_prot, mda.Writer(trajectory + '/system.xtc', u.atoms.n_atoms) as W_sys:
            for time, ts in enumerate(u.trajectory):
                with open(trajectory + '/write.log', 'a') as f:
                    f.write(str(time) + '\n')
                W_prot.write(u.select_atoms('protein'))
                W_sys.write(u.atoms)

                
        u.select_atoms('protein').write(trajectory + '/protein.pdb')
        u.atoms.write(trajectory + '/system.pdb')
                
                
        with open(self.filename + '/log.log', 'a') as f:
            f.write(trajectory + ' done' + '\n')
        # return u
        with open(self.filename + '_'.join(trajectory.split('/')[-3:-1]) + '.pickle', 'wb') as f:
            pickle.dump(u, f)
        del u
        del u_bond
        gc.collect()
        return self.filename + '_'.join(trajectory.split('/')[-3:-1]) + '.pickle'

    def load_preprocessing_trajectory(self,trajectory):
        return self.filename + '_'.join(trajectory.split('/')[-3:-1]) + '.pickle'

    def load_protein(self, trajectory):
        #    print(trajectory)
        with open(self.filename + '/log.log', 'a') as f:
            f.write(trajectory + '\n')
        u = mda.Universe(trajectory + '/protein.pdb',
                        trajectory + '/protein.xtc')

        with open(self.filename + '/log.log', 'a') as f:
            f.write(trajectory + ' done' + '\n')
        # return u
        with open(self.filename + '_'.join(trajectory.split('/')[-3:-1]) + '_prot.pickle', 'wb') as f:
            pickle.dump(u, f)

        return self.filename + '_'.join(trajectory.split('/')[-3:-1]) + '_prot.pickle'

    def load_system(self, trajectory):
        #    print(trajectory)
        with open(self.filename + '/log.log', 'a') as f:
            f.write(trajectory + '\n')
        u = mda.Universe(trajectory + '/ca.pdb',
                        trajectory + '/system.xtc')

        with open(self.filename + '/log.log', 'a') as f:
            f.write(trajectory + ' done' + '\n')
    # return u
        with open(self.filename + '_'.join(trajectory.split('/')[-3:-1]) + '_sys.pickle', 'wb') as f:
            pickle.dump(u, f)

        return self.filename + '_'.join(trajectory.split('/')[-3:-1]) + '_sys.pickle'

    @property
    def filename(self):
        return self.pwd + '/' + self.ensemble_name