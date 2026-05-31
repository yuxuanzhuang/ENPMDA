#!/usr/bin/env python

"""Tests for `ENPMDA` package."""
import os
import pickle
import shutil
import tempfile
import pytest
import numpy as np
from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_equal,
    assert_array_equal,
)

from ENPMDA import MDDataFrame
from ENPMDA.preprocessing import TrajectoryEnsemble
from ENPMDATests.datafiles import (
    ensemble_ala_bonded,
    ensemble_ala_traj,
    ensemble_ala_top
)


class TestDDataFrameCreation(object):
    @pytest.fixture
    def tempdir(self):
        tempdir = tempfile.mkdtemp()
        return tempdir

    def test_inititialize_trajectoryensemble(self, tempdir):
        traj_ensemble = TrajectoryEnsemble(ensemble_name= './' + tempdir + '/test_traj_ensemble',
                                            topology_list=ensemble_ala_top,
                                            trajectory_list=ensemble_ala_traj,
                                            bonded_topology_list=ensemble_ala_bonded,
                                            updating=False,
                                            only_raw=False)
        traj_ensemble.load_ensemble()
        assert_equal(traj_ensemble.filename,
                     os.getcwd() + tempdir + '/' + 'test_traj_ensemble' + '/',
                     "ensemble name is not set correctly")

    def test_initialize_mddataframe(self, tempdir):
        md_dataframe = MDDataFrame(dataframe_name='./' + tempdir + '/test_init_dataframe')

        assert_equal(md_dataframe.filename,
                     os.getcwd() + tempdir + '/' + 'test_init_dataframe' + '/',
                     "working dir is not set correctly")

    def test_init_trajectoryensemble_abs_path(self, tempdir):
        traj_ensemble = TrajectoryEnsemble(ensemble_name=tempdir + '/test_traj_ensemble',
                                            topology_list=ensemble_ala_top,
                                            trajectory_list=ensemble_ala_traj,
                                            bonded_topology_list=ensemble_ala_bonded,
                                            updating=False,
                                            only_raw=False)
        traj_ensemble.load_ensemble()
        assert_equal(traj_ensemble.filename,
                     tempdir + '/' + 'test_traj_ensemble' + '/',
                     "ensemble name is not set correctly")

    def test_init_mddataframe_abs_path(self, tempdir):
        md_dataframe = MDDataFrame(dataframe_name=tempdir + '/test_init_dataframe')

        assert_equal(md_dataframe.filename,
                     tempdir + '/' + 'test_init_dataframe' + '/',
                     "working dir is not set correctly")


class TestAddTrajEnsemble(object):

    @pytest.fixture
    def tempdir(self):
        tempdir = tempfile.mkdtemp()
        return tempdir

    @pytest.fixture()
    def md_dataframe(self, tempdir):
        return MDDataFrame(dataframe_name=tempdir + '/test_init_dataframe')

    @pytest.fixture()
    def traj_ensemble(self, tempdir):
        traj_ensemble = TrajectoryEnsemble(ensemble_name=tempdir + '/test_ensemble',
                                  topology_list=ensemble_ala_top,
                                  trajectory_list=ensemble_ala_traj,
                                  bonded_topology_list=ensemble_ala_bonded,
                                  updating=False,
                                  only_raw=False)
        traj_ensemble.load_ensemble()
        return traj_ensemble

    @pytest.fixture()
    def raw_traj_ensemble(self, tempdir):
        traj_ensemble = TrajectoryEnsemble(ensemble_name=tempdir + '/test_ensemble',
                                  topology_list=ensemble_ala_top,
                                  trajectory_list=ensemble_ala_traj,
                                  bonded_topology_list=ensemble_ala_bonded,
                                  updating=False,
                                  only_raw=True)
        traj_ensemble.load_ensemble()
        return traj_ensemble

    def test_add_trajectory_ensemble(self, md_dataframe, traj_ensemble):
        md_dataframe.add_traj_ensemble(traj_ensemble, npartitions=10)

        assert md_dataframe.dd_dataframe is not None
        assert_array_equal(md_dataframe.dataframe.shape,
                           (168, 7),
                           "Dataframe shape is not correct")
        assert_equal(md_dataframe.dataframe.columns.tolist(),
                     ['universe_protein', 'universe_system', 'system',
                         'traj_name', 'frame', 'traj_time', 'stride'],
                     "Dataframe columns are not correct")

        assert_equal(md_dataframe.npartitions, 10,
                     "npartitions not set correctly")
        assert md_dataframe.analysis_results is not None

    def test_add_raw_trajectory_ensemble(self, md_dataframe, raw_traj_ensemble):
        md_dataframe.add_traj_ensemble(raw_traj_ensemble, npartitions=10)

        assert md_dataframe.dd_dataframe is not None
        assert_array_equal(md_dataframe.dataframe.shape,
                           (168, 7),
                           "Dataframe shape is not correct")
        assert_equal(md_dataframe.dataframe.columns.tolist(),
                     ['universe_protein', 'universe_system', 'system',
                         'traj_name', 'frame', 'traj_time', 'stride'],
                     "Dataframe columns are not correct")

        assert_equal(md_dataframe.npartitions, 10,
                     "npartitions not set correctly")
        assert md_dataframe.analysis_results is not None

    def test_mixed_chained_trajectory_list(self, tempdir):
        topology1 = os.path.join(tempdir, 'start1.pdb')
        topology2 = os.path.join(tempdir, 'start2.pdb')
        topology3 = os.path.join(tempdir, 'start3.pdb')
        traj1 = os.path.join(tempdir, 'md1.xtc')
        traj2a = os.path.join(tempdir, 'md2a.xtc')
        traj2b = os.path.join(tempdir, 'md2b.xtc')
        traj3 = os.path.join(tempdir, 'md3.xtc')
        shutil.copyfile(ensemble_ala_top[0], topology1)
        shutil.copyfile(ensemble_ala_top[1], topology2)
        shutil.copyfile(ensemble_ala_top[2], topology3)
        shutil.copyfile(ensemble_ala_traj[0], traj1)
        shutil.copyfile(ensemble_ala_traj[1], traj2a)
        shutil.copyfile(ensemble_ala_traj[2], traj2b)
        shutil.copyfile(ensemble_ala_traj[3], traj3)

        traj_ensemble = TrajectoryEnsemble(
            ensemble_name=tempdir + '/test_chain_ensemble',
            topology_list=[topology1, topology2, topology3],
            trajectory_list=[traj1, [traj2a, traj2b], traj3],
            updating=True,
            only_raw=True,
        )
        traj_ensemble.load_ensemble()

        assert_equal(len(traj_ensemble.trajectory_files), 3)
        with open(traj_ensemble.trajectory_files[0], 'rb') as f:
            universe1 = pickle.load(f)
        with open(traj_ensemble.trajectory_files[1], 'rb') as f:
            universe2 = pickle.load(f)
        with open(traj_ensemble.trajectory_files[2], 'rb') as f:
            universe3 = pickle.load(f)
        assert_equal(universe1.trajectory.n_frames, 21)
        assert_equal(universe2.trajectory.n_frames, 42)
        assert_equal(len(universe2.trajectory.filenames), 2)
        assert_equal(universe3.trajectory.n_frames, 21)
