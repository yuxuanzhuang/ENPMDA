from pkg_resources import resource_filename

ensemble_ala_tpr = []
ensemble_ala_traj = []
ensemble_ala_top = []

for rep in range(1, 9):
    ensemble_ala_tpr.append(resource_filename(__name__,
                                              f'datafile/ensemble_AlaDipeptide/rep{rep}/md.tpr'))
    ensemble_ala_traj.append(resource_filename(__name__,
                                               f'datafile/ensemble_AlaDipeptide/rep{rep}/md.xtc'))
    ensemble_ala_top.append(resource_filename(__name__,
                                              f'datafile/ensemble_AlaDipeptide/rep{rep}/start.pdb'))
