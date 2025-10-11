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
    ensemble_ala_bonded,
    ensemble_ala_traj,
    ensemble_ala_top,
)

def _normpath(p: str) -> str:
    """Normalize for assertions (resolves /private/var -> /var on mac)."""
    # preserve trailing slash behavior used in the code/tests
    trailing = p.endswith(os.sep)
    rp = os.path.realpath(os.path.normpath(p))
    return rp + (os.sep if trailing and not rp.endswith(os.sep) else "")


class TestDDataFrameCreation(object):
    @pytest.fixture
    def tempdir(self):
        return tempfile.mkdtemp()

    def test_inititialize_trajectoryensemble(self, tempdir):
        # Use absolute path; don't prefix with './'
        ensemble_path = os.path.join(tempdir, 'test_traj_ensemble')
        traj_ensemble = TrajectoryEnsemble(
            ensemble_name=ensemble_path,
            topology_list=ensemble_ala_top,
            trajectory_list=ensemble_ala_traj,
            bonded_topology_list=ensemble_ala_bonded,
            updating=False,
            only_raw=False,
        )
        traj_ensemble.load_ensemble()

        # TrajectoryEnsemble.filename includes /skip<unique_skips>/
        expected = os.path.join(ensemble_path, 'skip1') + '/'
        assert_equal(
            traj_ensemble.filename,
            expected,
            "ensemble filename path is not set correctly",
        )

    def test_initialize_mddataframe(self, tempdir):
        df_path = os.path.join(tempdir, 'test_init_dataframe')
        md_dataframe = MDDataFrame(dataframe_name=df_path)

        # filename = <abs_dir>/<name>/
        expected = df_path + '/'
        assert_equal(
            md_dataframe.filename,
            expected,
            "MDDataFrame filename (save dir) is not set correctly",
        )

    def test_init_trajectoryensemble_abs_path(self, tempdir):
        ensemble_path = os.path.join(tempdir, 'test_traj_ensemble')
        traj_ensemble = TrajectoryEnsemble(
            ensemble_name=ensemble_path,
            topology_list=ensemble_ala_top,
            trajectory_list=ensemble_ala_traj,
            bonded_topology_list=ensemble_ala_bonded,
            updating=False,
            only_raw=False,
        )
        traj_ensemble.load_ensemble()

        expected = os.path.join(ensemble_path, 'skip1') + '/'
        assert_equal(
            traj_ensemble.filename,
            expected,
            "ensemble filename path is not set correctly (abs path case)",
        )

    def test_init_mddataframe_abs_path(self, tempdir):
        df_path = os.path.join(tempdir, 'test_init_dataframe')
        md_dataframe = MDDataFrame(dataframe_name=df_path)

        assert_equal(
            md_dataframe.filename,
            df_path + '/',
            "MDDataFrame filename is not set correctly (abs path case)",
        )


class TestAddTrajEnsemble(object):
    @pytest.fixture
    def tempdir(self):
        return tempfile.mkdtemp()

    @pytest.fixture()
    def md_dataframe(self, tempdir):
        return MDDataFrame(dataframe_name=os.path.join(tempdir, 'test_init_dataframe'))

    @pytest.fixture()
    def traj_ensemble(self, tempdir):
        ensemble_path = os.path.join(tempdir, 'test_ensemble')
        traj_ensemble = TrajectoryEnsemble(
            ensemble_name=ensemble_path,
            topology_list=ensemble_ala_top,
            trajectory_list=ensemble_ala_traj,
            bonded_topology_list=ensemble_ala_bonded,
            updating=False,
            only_raw=False,
        )
        traj_ensemble.load_ensemble()
        return traj_ensemble

    @pytest.fixture()
    def raw_traj_ensemble(self, tempdir):
        ensemble_path = os.path.join(tempdir, 'test_ensemble')
        traj_ensemble = TrajectoryEnsemble(
            ensemble_name=ensemble_path,
            topology_list=ensemble_ala_top,
            trajectory_list=ensemble_ala_traj,
            bonded_topology_list=ensemble_ala_bonded,
            updating=False,
            only_raw=True,
        )
        traj_ensemble.load_ensemble()
        return traj_ensemble

    def test_add_trajectory_ensemble(self, md_dataframe, traj_ensemble):
        md_dataframe.add_traj_ensemble(traj_ensemble, npartitions=10)

        assert md_dataframe.dd_dataframe is not None
        assert_array_equal(
            md_dataframe.dataframe.shape,
            (168, 7),
            "Dataframe shape is not correct",
        )
        assert_equal(
            md_dataframe.dataframe.columns.tolist(),
            [
                'universe_protein',
                'universe_system',
                'system',
                'traj_name',
                'frame',
                'traj_time',
                'stride',
            ],
            "Dataframe columns are not correct",
        )

        assert_equal(md_dataframe.npartitions, 10, "npartitions not set correctly")
        assert md_dataframe.analysis_results is not None

    def test_add_raw_trajectory_ensemble(self, md_dataframe, raw_traj_ensemble):
        md_dataframe.add_traj_ensemble(raw_traj_ensemble, npartitions=10)

        assert md_dataframe.dd_dataframe is not None
        assert_array_equal(
            md_dataframe.dataframe.shape,
            (168, 7),
            "Dataframe shape is not correct",
        )
        assert_equal(
            md_dataframe.dataframe.columns.tolist(),
            [
                'universe_protein',
                'universe_system',
                'system',
                'traj_name',
                'frame',
                'traj_time',
                'stride',
            ],
            "Dataframe columns are not correct",
        )

        assert_equal(md_dataframe.npartitions, 10, "npartitions not set correctly")
        assert md_dataframe.analysis_results is not None


class TestPathNormalizationEdgeCases:
    @pytest.fixture
    def tempdir(self):
        return tempfile.mkdtemp()

    def test_traj_ensemble_with_dot_slash_absolute(self, tempdir):
        # Even if user passes './' + ABS, we normalize to ABS
        abs_path = os.path.join(tempdir, 'ensemble_dot_abs')
        weird_input = '.' + os.sep + abs_path  # './<abs>'
        te = TrajectoryEnsemble(
            ensemble_name=weird_input,
            topology_list=ensemble_ala_top,
            trajectory_list=ensemble_ala_traj,
            bonded_topology_list=ensemble_ala_bonded,
            updating=False,
            only_raw=True,
        )
        te.load_ensemble()
        expected = os.path.join(abs_path, 'skip1') + '/'
        assert_equal(_normpath(te.filename),
                    _normpath(expected), "Dot-slash + absolute was not normalized")

    def test_traj_ensemble_multi_skip_unique_sort(self, tempdir):
        # When multiple skip values are present, folder should be 'skip1_2' (sorted unique)
        abs_path = os.path.join(tempdir, 'ensemble_multi_skip')
        te = TrajectoryEnsemble(
            ensemble_name=abs_path,
            topology_list=ensemble_ala_top,
            trajectory_list=ensemble_ala_traj,
            bonded_topology_list=ensemble_ala_bonded,
            skip=[2 if i % 2 else 1 for i in range(len(ensemble_ala_traj))],
            updating=False,
            only_raw=True,
        )
        te.load_ensemble()
        expected = os.path.join(abs_path, 'skip1_2') + '/'
        assert_equal(_normpath(te.filename),
                    _normpath(expected), "Multi-skip folder name is incorrect")

    def test_mddataframe_with_relative_input(self, tempdir, monkeypatch):
        # Change CWD to tempdir and pass a relative name; filename should resolve to ABS
        subdir_name = 'relative_df'
        monkeypatch.chdir(tempdir)
        md = MDDataFrame(dataframe_name=subdir_name)
        expected = os.path.join(tempdir, subdir_name) + '/'
        assert_equal(_normpath(md.filename),
                     _normpath(expected),
                     "Relative MDDataFrame path did not resolve to ABS")

    def test_mddataframe_dot_slash_absolute(self, tempdir):
        abs_path = os.path.join(tempdir, 'df_dot_abs')
        weird_input = '.' + os.sep + abs_path
        md = MDDataFrame(dataframe_name=weird_input)
        expected = abs_path + '/'
        assert_equal(
                _normpath(md.filename),
                _normpath(expected),
                "Dot-slash + absolute MDDataFrame not normalized")

    def test_mddataframe_filename_trailing_sep(self, tempdir):
        df_path = os.path.join(tempdir, 'df_trailing')
        md = MDDataFrame(dataframe_name=df_path)
        assert md.filename.endswith('/'), "MDDataFrame.filename must end with '/'"

    def test_handle_same_folder_trajectory_symlinks(self, tempdir):
        """
        If multiple trajectories live in the same folder, _handle_same_folder_trajectory
        should create per-file subfolders and symlink the trajectory there.
        """
        # Make a fake folder with two dummy traj files and matching tops
        traj_dir = os.path.join(tempdir, 'sys')
        os.makedirs(traj_dir, exist_ok=True)
        t1 = os.path.join(traj_dir, 'a.xtc')
        t2 = os.path.join(traj_dir, 'b.xtc')
        top1 = os.path.join(traj_dir, 'a.pdb')
        top2 = os.path.join(traj_dir, 'b.pdb')
        # create dummy files
        for p in [t1, t2, top1, top2]:
            with open(p, 'wb') as fh:
                fh.write(b'\x00')

        abs_path = os.path.join(tempdir, 'ensemble_symlink_test')
        te = TrajectoryEnsemble(
            ensemble_name=abs_path,
            topology_list=[top1, top2],
            trajectory_list=[t1, t2],
            bonded_topology_list=None,
            updating=False,
            only_raw=True,
        )

        # Call the internal routine directly to avoid MDAnalysis processing
        te._handle_same_folder_trajectory()

        # Expect the two new subfolders and symlinks
        for src in [t1, t2]:
            base = os.path.splitext(os.path.basename(src))[0]
            new_dir = os.path.join(traj_dir, base)
            new_link = os.path.join(new_dir, os.path.basename(src))
            assert os.path.isdir(new_dir), "Subfolder for same-folder trajectory not created"
            assert os.path.islink(new_link) or os.path.exists(new_link), "Symlink for trajectory not created"

        # Also expect TrajectoryEnsemble.trajectory_list to be updated to the symlinked paths
        for new_path in te.trajectory_list:
            assert traj_dir in new_path and os.path.exists(new_path), "Updated trajectory path is invalid"


class TestMDDataFrameSaveLoad:
    @pytest.fixture
    def tempdir(self):
        return tempfile.mkdtemp()

    def test_save_and_load_restores_paths(self, tempdir):
        df_path = os.path.join(tempdir, 'df_saveload')
        md = MDDataFrame(dataframe_name=df_path)

        # minimal data to allow save()
        md.dataframe.loc[0, :] = [None, None, 0, 'traj', 0, 0.0, 1]
        md.computed = True

        # Save with base name == directory name ('df_saveload')
        base = 'df_saveload'
        md.save(name=base, overwrite=True)

        # Candidates:
        # 1) base form: <dir>/<base>  → expects <dir>/<base>/<base>_md_dataframe.pickle
        by_dir_and_base = os.path.join(tempdir, base)           # /tmp/.../df_saveload
        # 2) explicit *_md_dataframe.pickle
        explicit_obj_pickle = os.path.join(df_path, f'{base}_md_dataframe.pickle')
        # 3) plain .pickle (fallback)
        explicit_df_pickle = os.path.join(df_path, f'{base}.pickle')

        for candidate in (by_dir_and_base, explicit_obj_pickle, explicit_df_pickle):
            loaded = MDDataFrame.load_dataframe(candidate)
            assert isinstance(loaded, MDDataFrame)
            assert_equal(
                _normpath(loaded.filename),
                _normpath(md.filename),
                f"Loaded MDDataFrame filename mismatch for candidate: {candidate}",
            )