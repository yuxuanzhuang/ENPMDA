import numpy as np
import itertools

from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.lib.distances import calc_bonds
from MDAnalysis.analysis.distances import self_distance_array
from MDAnalysis.analysis.dihedrals import Dihedral as mda_dihedral

from .base import DaskChunkMdanalysis


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


class get_atomic_position(DaskChunkMdanalysis):
    name = 'at_position'

    def run_analysis(self, universe, start, stop, step):
        backbone_atoms = universe.select_atoms('backbone')

        self._feature_info = ['_'.join([str(residex), feat]) for residex, feat in itertools.product(
            backbone_atoms.residues.resindices, ['_bb_pos_x', '_bb_pos_y', '_bb_pos_z'])]

        result = []
        for ts in universe.trajectory[start:stop:step]:
            result.append(
                (backbone_atoms.positions -
                 backbone_atoms.center_of_geometry()).ravel())

        return result


class get_rmsd_init(DaskChunkMdanalysis):
    name = 'rmsd_to_init'

    def run_analysis(self, universe, start, stop, step):
        self._feature_info = ['rmsd_to_init']

        name_backbone = universe.select_atoms('backbone')
        rmsd = RMSD(name_backbone, name_backbone).run(start, stop, step)
        return rmsd.results['rmsd'].T[2]


class get_protein_hydration(DaskChunkMdanalysis):
    name = 'protein_hydration'
    universe_file = 'system'

    def run_analysis(self, universe, start, stop, step):
        self._feature_info = ['protein_hydration']

        prot_hydration = universe.select_atoms(
            "around 5 backbone", updating=True)

        n_hydration = []
        for ts in universe.trajectory[start:stop:step]:
            n_hydration.append(prot_hydration.n_atoms)
        return n_hydration
