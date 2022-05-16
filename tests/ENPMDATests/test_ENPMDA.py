#!/usr/bin/env python

"""Tests for `ENPMDA` package."""
import os
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
    ensemble_ala_tpr,
    ensemble_ala_traj,
    ensemble_ala_top
)


class TestDDataFrameCreation(object):
    @pytest.fixture
    def tempdir(self):
        tempdir = tempfile.mkdtemp()
        return tempdir

    def test_inititialize_trajectoryensemble(self, tempdir):
        traj_ensembles = TrajectoryEnsemble(ensemble_name=tempdir + '/test_traj_ensemble',
                                            topology_list=ensemble_ala_top,
                                            trajectory_list=ensemble_ala_traj,
                                            tpr_list=ensemble_ala_tpr,
                                            updating=False,
                                            only_raw=False)
        assert_equal(traj_ensembles.filename,
                     os.getcwd() + '//' + tempdir + '/' + 'test_traj_ensemble' + '/',
                     "ensemble name is not set correctly")

    def test_initialize_mddataframe(self, tempdir):
        md_dataframe = MDDataFrame(tempdir + '/test_init_dataframe')

        assert_equal(md_dataframe.working_dir,
                     os.getcwd() + '/' + tempdir + '/' + 'test_init_dataframe' + '/',
                     "working dir is not set correctly")


class TestAddTrajEnsemble(object):

    @pytest.fixture
    def tempdir(self):
        tempdir = tempfile.mkdtemp()
        return tempdir

    @pytest.fixture()
    def md_dataframe(self, tempdir):
        return MDDataFrame(tempdir + '/test_init_dataframe')

    @pytest.fixture()
    def traj_ensemble(self, tempdir):
        return TrajectoryEnsemble(ensemble_name=tempdir + '/test_ensemble',
                                  topology_list=ensemble_ala_top,
                                  trajectory_list=ensemble_ala_traj,
                                  tpr_list=ensemble_ala_tpr,
                                  updating=False,
                                  only_raw=False)

    def test_add_trajectory_ensemble(self, md_dataframe, traj_ensemble):
        md_dataframe.add_traj_ensemble(traj_ensemble, npartitions=10)

        assert md_dataframe.dd_dataframe is not None
        assert_array_equal(md_dataframe.dataframe.shape,
                           (168, 7),
                           "Dataframe shape is not correct")
        assert_equal(md_dataframe.dataframe.columns.tolist(),
                     ['universe', 'universe_system', 'system',
                         'md_name', 'frame', 'traj_time', 'stride'],
                     "Dataframe columns are not correct")

        assert_equal(md_dataframe.npartitions, 10,
                     "npartitions not set correctly")
        assert md_dataframe.analysis_results is not None
