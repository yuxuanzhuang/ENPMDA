import numpy as np
import itertools

from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.dihedrals import Dihedral as mda_dihedral

from .base import DaskChunkMdanalysis


class get_backbonetorsion(DaskChunkMdanalysis):
    name = "torsion"

    def set_feature_info(self, universe):
        prot_ag = universe.select_atoms("protein")
        feature_info = []
        res_phi_list = [
            x.resindices[0] for x in prot_ag.residues.phi_selections() if x is not None
        ]

        res_psi_list = [
            x.resindices[0] for x in prot_ag.residues.psi_selections() if x is not None
        ]

        feature_info.extend(
            [
                "_".join([str(residex), ang])
                for residex, ang in itertools.product(
                    res_phi_list, ["phi_cos", "phi_sin"]
                )
            ]
        )
        feature_info.extend(
            [
                "_".join([str(residex), ang])
                for residex, ang in itertools.product(
                    res_psi_list, ["psi_cos", "psi_sin"]
                )
            ]
        )

        return feature_info

    def run_analysis(self, universe, start, stop, step):
        prot_ag = universe.select_atoms("protein")
        phi_ag = [x for x in prot_ag.residues.phi_selections() if x is not None]
        phi_angle_ana = mda_dihedral(phi_ag).run(start, stop, step)
        phi_angle = np.concatenate(
            [
                [np.cos(np.deg2rad(phi_angle_ana.results["angles"]))],
                [np.sin(np.deg2rad(phi_angle_ana.results["angles"]))],
            ]
        ).transpose(1, 0, 2)

        psi_ag = [x for x in universe.atoms.residues.psi_selections() if x is not None]
        psi_angle_ana = mda_dihedral(psi_ag).run(start, stop, step)
        psi_angle = np.concatenate(
            [
                [np.cos(np.deg2rad(psi_angle_ana.results["angles"]))],
                [np.sin(np.deg2rad(psi_angle_ana.results["angles"]))],
            ]
        ).transpose(1, 0, 2)

        #TODO: psi_angle and phi_angle may have different shape
        try:
            torsion_angle = np.concatenate([psi_angle, phi_angle], axis=1)
        except ValueError:
            torsion_angle = np.stack([phi_angle, psi_angle], axis=1)
        torsion_angle = torsion_angle.reshape(
            torsion_angle.shape[0], torsion_angle.shape[1] * torsion_angle.shape[2]
        )

        return torsion_angle


class get_atomic_position(DaskChunkMdanalysis):
    name = "ca_position"

    def set_feature_info(self, universe):
        ca_atoms = universe.select_atoms("protein and name CA")
        return [
            "_".join([str(residex), feat])
            for residex, feat in itertools.product(
                ca_atoms.resindices, ["_ca_pos_x", "_ca_pos_y", "_ca_pos_z"]
            )
        ]

    def run_analysis(self, universe, start, stop, step):
        ca_atoms = universe.select_atoms("protein and name CA")
        result = []
        for _ in universe.trajectory[start:stop:step]:
            result.append(
                (ca_atoms.positions - ca_atoms.center_of_geometry()).ravel()
            )

        return result


class get_rmsd_init(DaskChunkMdanalysis):
    name = "rmsd_2_init"

    def set_feature_info(self, universe):
        return ["frame0"]

    def run_analysis(self, universe, start, stop, step):
        name_backbone = universe.select_atoms("protein and backbone")
        rmsd = RMSD(name_backbone, name_backbone).run(start, stop, step)
        return rmsd.results["rmsd"].T[2]


class get_protein_hydration(DaskChunkMdanalysis):
    name = "protein_hydration"
    universe_file = "system"

    def set_feature_info(self, universe):
        return ["protein_hydration"]

    def run_analysis(self, universe, start, stop, step):
        prot_hydration = universe.select_atoms("around 5 backbone", updating=True)

        n_hydration = []
        for _ in universe.trajectory[start:stop:step]:
            n_hydration.append(prot_hydration.n_atoms)
        return n_hydration
