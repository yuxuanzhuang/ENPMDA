#!/usr/bin/env python

"""Tests for `ENPMDA` package."""
import os
from pathlib.Path import tmp_path
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
from ENPMDATests.datafiles import *

topology_list = [CLIMBER_BGT_EPJ, CLIMBER_BGT_EPJPNU,
                 CLIMBER_EPJPNU_BGT, CLIMBER_EPJPNU_EPJ,
                 CLIMBER_EPJ_BGT, CLIMBER_EPJ_EPJPNU]

trajectory_list = [CLIMBER_BGT_EPJ_TRANSITION, CLIMBER_BGT_EPJPNU_TRANSITION,
                 CLIMBER_EPJPNU_BGT_TRANSITION, CLIMBER_EPJPNU_EPJ_TRANSITION,
                 CLIMBER_EPJ_BGT_TRANSITION, CLIMBER_EPJ_EPJPNU_TRANSITION]

class TestDDataFrameCreation(object):

    @pytest.fixture
    def my_filepath(self, tmpdir):
        return tmpdir.mkdir("sub").join("testCurrentTicketCount.txt")


    def test_inititialize_trajectoryensemble(self):
        traj_ensembles =  TrajectoryEnsemble(ensemble_name='test_traj_ensemble',
                                    topology_list=topology_list,
                                    trajectory_list=trajectory_list,
                                    updating=True)

        assert_equal(traj_ensembles.ensemble_name, 'test_traj_ensemble', "ensemble name is not set correctly")


    def test_initialize_mddataframe(self):
        # Universe(top, trj)
        tmp_dir = tmp_path
        tmp_dir.mkdir()
        md_dataframe = MDDataFrame(tmp_path'test_init_dataframe')

        assert_equal(md_dataframe.working_dir, os.getcwd() + '/' + 'test_init_dataframe' + '/', "working dir is not set correctly")


class TestAddTrajEnsemble(object):
    @staticmethod
    @pytest.fixture()
    def md_dataframe():
        return MDDataFrame('test_init_dataframe')

    @staticmethod
    @pytest.fixture()
    def traj_ensemble():
        return TrajectoryEnsemble(ensemble_name='test_traj_ensemble',
                                    topology_list=topology_list,
                                    trajectory_list=trajectory_list,
                                    updating=False)

    def test_add_trajectory_ensemble(self, md_dataframe, traj_ensemble):
        md_dataframe.add_traj_ensemble(traj_ensemble, npartitions=10)

        assert md_dataframe.dd_dataframe is not None
        assert_array_equal(md_dataframe.dataframe.shape, (1806,7), "Dataframe shape is not correct")
        assert_equal(md_dataframe.dataframe.columns.tolist(),
        ['universe', 'universe_system', 'system', 'md_name', 'frame', 'traj_time', 'stride'],
        "Dataframe columns are not correct")

        assert_equal(md_dataframe.npartitions, 10, "npartitions not set correctly")
        assert md_dataframe.analysis_results is not None
