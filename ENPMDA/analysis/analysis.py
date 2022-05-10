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
                for i, j in self.arg_comp_candidates:
                    filtered_distance_matrix.append(distance_matrix[i, j])

                r_ts.extend(filtered_distance_matrix)
            result.append(r_ts)

        self._feature_info = []
        for selection1, selection2 in selection_comb:
            for feature_chain in self.feature_list:
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
                for i, j in self.filtered_candidate:
                    filtered_distance_matrix.append(
                        1.0 / distance_matrix[i, j])

                r_ts.extend(filtered_distance_matrix)
            result.append(r_ts)

        self._feature_info = []
        for selection1, selection2 in selection_comb:
            for feature_chain in self.feature_list:
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
            self.u_ref.select_atoms('name CA'),
            ref_frame=0).run(
            start,
            stop,
            step)
        rmsd_epj = RMSD(
            name_ca,
            self.u_ref.select_atoms('name CA'),
            ref_frame=1).run(
            start,
            stop,
            step)
        rmsd_pnu = RMSD(
            name_ca,
            self.u_ref.select_atoms('name CA'),
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