"""\
===========
MDDataFrame
===========
The :class:`~ENPMDA.MDDataFrame` class both store
the metadata of simulations in the ensemble and functions as
a dask dataframe to add, compute, and store analysis.

A ``MDDataFrame`` is created from files::

    from ENPMDA import MDDataFrame
    md_dataframe = MDDataFrame()
    md_dataframe.add_traj_ensemble(traj_ensemble, npartitions=16)


Classes
=======
.. autoclass:: MDDataFrame
   :members:
"""


from datetime import datetime
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

import numpy as np

# import awkward as ak

import dask.dataframe as dd
import dask
import pandas as pd
import MDAnalysis as mda
import os
import pickle
import shutil
import gc
from tqdm import tqdm
from sklearn import preprocessing


from ENPMDA.analysis.base import AnalysisResult
from ENPMDA.preprocessing import TrajectoryEnsemble

timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
meta_data_list = [
    "universe_protein",
    "universe_system",
    "system",
    "traj_name",
    "frame",
    "traj_time",
    "stride",
]


class MDDataFrame(object):
    r"""
    Class to store the metadata and analysis results
    of the ensemble simulations.

    It uses pandas.DataFrame to store metadata
    and dask.DataFrame to distribute computation jobs
    so that the parallel analysis can be performed not
    only for one trajectory but also across simulations
    and analyses.
    """

    def __init__(
        self, dataframe_name, meta_data_list=meta_data_list, timestamp=timestamp
    ):
        """
        Parameters
        ----------
        dataframe_name: str
            The name of the dataframe
            It will be used as the folder to save
            all the analysis results.
            It can also be the absolute path to the folder.

        meta_data_list: list, optional
            List of metadata in the dataframe.
            In default, the locations of pickled universes
            of protein and system, the system index, the
            trajectory filename, the frame index, the
            trajectory time, and the stride are stored.

        timestamp: str, optional
            The timestamp of creating the ensemble
            It will be set to the current time if not provided.
        """
        self.dataframe_name = dataframe_name
        self.dataframe = pd.DataFrame(columns=meta_data_list)
        self.computed = False
        self.sorted = False

        # set working dir to absolute directory
        if not os.path.isabs(self.dataframe_name):
            self.working_dir = os.getcwd() + "/"
        else:
            self.working_dir = os.path.relpath(self.dataframe, os.getcwd())

        # the directory used for the first time
        self.init_dir = self.working_dir
        self.timestamp = timestamp
        self.trajectory_ensemble = None
        self.analysis_list = []

    def add_traj_ensemble(
        self, trajectory_ensemble: TrajectoryEnsemble, npartitions, stride=1
    ):
        """
        Parameters
        ----------
        trajectory_ensemble: ENPMDA.TrajectoryEnsemble
            The trajectory ensemble to be added to the
            dataframe.

        npartitions: int
            The number of partitions to be used in
            the dask dataframe.

        stride: int, optional
            The stride to be used in the dask dataframe.
            It is used to skip frames in the trajectory.
        """
        if self.trajectory_ensemble is not None:
            raise ValueError("Trajectory ensemble already added")

        self.trajectory_ensemble = trajectory_ensemble

        if trajectory_ensemble.protein_trajectory_files is None:
            warnings.warn(
                "The provided trajectory ensemble "
                "only contain raw trajectories "
                "all analysis will be performed on the raw trajectories",
                stacklevel=2,
            )
            self.trajectory_files = trajectory_ensemble.trajectory_files
            self.protein_trajectory_files = trajectory_ensemble.trajectory_files
            self.system_trajectory_files = trajectory_ensemble.trajectory_files
        else:
            self.trajectory_files = trajectory_ensemble.trajectory_files
            self.protein_trajectory_files = trajectory_ensemble.protein_trajectory_files
            self.system_trajectory_files = trajectory_ensemble.system_trajectory_files

        self.npartitions = npartitions
        self.stride = stride

        meta_data_jobs = []
        for ind, trajectory in enumerate(self.protein_trajectory_files):
            meta_data_jobs.append(
                dask.delayed(self._append_metadata)(trajectory, system=ind)
            )

        meta_data = dask.compute(meta_data_jobs)[0]

        for i, trajectory in enumerate(self.protein_trajectory_files):
            self.dataframe = pd.concat(
                [
                    self.dataframe,
                    pd.DataFrame(meta_data[i], columns=self.dataframe.columns),
                ],
                ignore_index=True,
            )

        self.dataframe.frame = self.dataframe.frame.apply(int)
        self.dataframe.traj_time = self.dataframe.traj_time.apply(float)

        self._init_dd_dataframe()

    def _append_metadata(self, universe, system):
        universe_system = self.system_trajectory_files[system]

        u = pickle.load(open(universe, "rb"))
        u_sys = pickle.load(open(universe_system, "rb"))
        if u.trajectory.n_frames != u_sys.trajectory.n_frames:
            raise ValueError(
                f"In system {system}, number of frames in protein and system trajectories are different!"
            )
        rep_data = []

        md_name = u.trajectory.filename
        timestep = u.trajectory.dt

        for i in range(0, u.trajectory.n_frames, self.stride):
            rep_data.append(
                [
                    universe,
                    universe_system,
                    system,
                    md_name,
                    i,
                    i * timestep,
                    self.stride,
                ]
            )
        del u
        return rep_data

    def _init_dd_dataframe(self):
        self.dd_dataframe = dd.from_pandas(self.dataframe, npartitions=self.npartitions)
        print("Requested number of partitions: ", self.npartitions)
        if self.dd_dataframe.npartitions != self.npartitions:
            print("Actual {} partitions".format(self.dd_dataframe.npartitions))
            self.npartitions = self.dd_dataframe.npartitions
        self.analysis_results = AnalysisResult(
            self.dd_dataframe,
            self.dataframe,
            working_dir=self.filename,
            timestamp=self.timestamp,
        )

    def add_analysis(self, analysis, overwrite=False, **kwargs):
        """
        Add an analysis to the dataframe.

        Parameters
        ----------
        analysis: ENPMDA.analysis.base.DaskChunkMdanalysis
            The analysis to be added to the dataframe.

        overwrite: bool, optional
            Whether to overwrite the analysis if it is
            already in the dataframe.

        **kwargs: dict, optional
            Keyword arguments to be passed to the analysis.
        """
        self.computed = False
        self.sorted = False
        if analysis.name in self.analysis_list and not overwrite:
            warnings.warn(
                f"Analysis {analysis.name} already added, add overwrite=True to overwrite",
                stacklevel=2,
            )
        elif analysis.name in self.analysis_list and overwrite:
            warnings.warn(f"Analysis {analysis.name} overwrites!", stacklevel=2)
            self.analysis_results.add_column_to_results(analysis, **kwargs)
            self.analysis_list.remove(analysis.name)
            self.analysis_list.append(analysis.name)
            print(f"Analysis {analysis.name} overwritten")
        else:
            self.analysis_results.add_column_to_results(analysis, **kwargs)
            self.analysis_list.append(analysis.name)

            print(f"Analysis {analysis.name} added")

    def compute(self):
        """
        Compute the analysis results.
        It will be append the analysis results to
        the dataframe.
        """
        if not self.computed:
            self.analysis_results.compute()
        self.analysis_results.append_to_dataframe(self.dataframe)
        self.computed = True

        # reinstantiate the dask dataframe
        self._init_dd_dataframe()

    def get_feature_info(self, feature_name):
        """
        Get the information about a feature.

        Parameters
        ----------
        feature_name: str
            The name of the feature.
        """
        feat_info = np.load(
            self.analysis_results.filename + feature_name + "_feature_info.npy"
        )
        return feat_info

    def get_feature(
        self, feature_list, extra_metadata=[], in_memory=True, working_dir=None
    ):
        """
        Get the features from the dataframe.

        Parameters
        ----------
        feature_list: list of str
            The list of features to be extracted.
        extra_metadata: list of str, optional
            The list of extra metadata to be extracted.
        in_memory: bool, optional
            Whether to load the features in memory.
        """
        meta_data = ["system", "traj_name", "frame", "traj_time"] + extra_metadata
        if not self.computed:
            self.compute()

        if isinstance(feature_list, str):
            feature_list = [feature_list]
        for feature in feature_list:
            if feature not in self.analysis_list:
                raise ValueError(f"Feature {feature} not in analysis list")

        if in_memory:
            feature_dataframe = self.dataframe[meta_data].copy()
            for feature in feature_list:
                raw_data = np.concatenate(
                    [
                        np.load(
                            location.replace(self.init_dir, self.working_dir),
                            allow_pickle=True,
                        )
                        for location, df in tqdm(
                            self.dataframe.groupby(feature, sort=False),
                            desc="Loading feature {}".format(feature),
                            total=self.dataframe[feature].nunique(),
                        )
                    ]
                )
                feat_info = np.load(
                    self.analysis_results.filename + feature + "_feature_info.npy",
                    allow_pickle=True,
                )
                col_names = [feature + "_" + feat for feat in feat_info]

                if raw_data.ndim == 1 and len(feat_info) != 1:
                    raw_data_con = []
                    for raw_data_single in raw_data:
                        raw_data_con.append(list(raw_data_single))
                    raw_data_concat = pd.DataFrame(raw_data_con, columns=col_names)
                else:
                    raw_data = raw_data.reshape(raw_data.shape[0], -1)
                    raw_data_concat = pd.DataFrame(raw_data, columns=col_names)
                feature_dataframe = pd.concat(
                    [feature_dataframe, raw_data_concat], axis=1
                )
            return feature_dataframe.reset_index(drop=True)
        else:
            if not self.sorted:
                self.sort_analysis_result()
            feature_dataframe = pd.DataFrame(columns=meta_data + feature_list)

            for ind, (system, df) in tqdm(
                enumerate(self.dataframe.groupby("system", sort=False)),
                desc="Loading features",
                total=len(self.dataframe.system.unique()),
            ):
                feature_dataframe = pd.concat(
                    [
                        feature_dataframe,
                        pd.DataFrame(
                            [
                                [
                                    system,
                                    df.traj_name.values[-1],
                                    df.frame.values[-1],
                                    df.traj_time.values[-1],
                                ]
                                + [df[feat].values[-1] for feat in feature_list]
                            ],
                            columns=feature_dataframe.columns,
                        ),
                    ]
                )

            return feature_dataframe.reset_index(drop=True)

    def save(self, name="dataframe", overwrite=False):
        """
        Compute the analysis results and
        save the dataframe to a pickle file.

        Parameters
        ----------
        name: str, optional
            The name of the pickle file.
            It will be saved in the working directory.
        overwrite: bool, optional
            Whether to overwrite the file if it exists.
        """

        if not self.computed:
            self.compute()
        self.save_name = name

        if overwrite:
            self.dump(name)
            return

        if not os.path.exists(f"{self.filename}{name}.pickle"):
            self.dump(name)
        else:
            md_dataframe_old = pickle.load(
                open(f"{self.filename}{name}_md_dataframe.pickle", "rb")
            )
            md_data_old = md_dataframe_old.dataframe

            if set(md_data_old.universe_protein) != set(
                self.dataframe.universe_protein
            ):
                print("Seeds changed")
                self.dump(name)

            if md_data_old.shape[0] != self.dataframe.shape[0]:
                print("Trajectory length changed")
                self.dump(name)

            elif set(md_data_old.columns) != set(self.dataframe.columns):
                print("# features changed")

                old_cols = md_data_old.columns
                new_cols = self.dataframe.columns
                print("New: " + np.setdiff1d(new_cols, old_cols))

                old_extra_cols = np.setdiff1d(old_cols, new_cols)

                for old_extra_col in old_extra_cols:
                    self.analysis_list.append(old_extra_col)
                    shutil.copyfile(
                        f"{md_dataframe_old.analysis_results.filename}{old_extra_col}_feature_info.npy",
                        f"{self.analysis_results.filename}{old_extra_col}_feature_info.npy",
                    )

                extra_cols = np.setdiff1d(new_cols, old_cols)

                for extra_col in extra_cols:
                    md_data_old[extra_col] = self.dataframe[extra_col]

                print("Common: " + np.intersect1d(new_cols, old_cols))
                common_cols = np.intersect1d(new_cols, old_cols)

                for common_col in common_cols:
                    md_data_old[common_col] = self.dataframe[common_col]

                self.dataframe = md_data_old
                self.dump(name, backup=True)
            else:
                print("No changes")
                self.dump(name)

    def dump(self, filename, backup=False):
        if backup:
            try:
                shutil.copyfile(
                    f"{self.filename}{filename}.pickle",
                    f"{self.filename}{filename}_{self.timestamp}.pickle",
                )
            except FileNotFoundError:
                pass
            try:
                shutil.copyfile(
                    f"{self.filename}{filename}_md_dataframe.pickle",
                    f"{self.filename}{filename}_md_dataframe_{self.timestamp}.pickle",
                )
            except FileNotFoundError:
                pass

        with open(f"{self.filename}{filename}.pickle", "wb") as f:
            pickle.dump(self.dataframe, f)
        with open(f"{self.filename}{filename}_md_dataframe.pickle", "wb") as f:
            pickle.dump(self, f)

    def sort_analysis_result(self):
        if not self.computed:
            self.compute()

        if not self.sorted:
            for feature in self.analysis_list:
                if self.dataframe[feature][0].split("_")[-1] == "0.npy":
                    print(f"{feature} already sorted")
                    continue
                print(f"start to sort {feature}.")

                #                builder = ak.ArrayBuilder()
                #                for location, df in self.dataframe.groupby(feature, sort=False):
                #                    builder.append(np.load(location, allow_pickle=True))

                old_locations = [
                    location
                    for location, df in self.dataframe.groupby(feature, sort=False)
                ]
                raw_data = np.concatenate(
                    [
                        np.load(location, allow_pickle=True)
                        for location in old_locations
                    ],
                    axis=0,
                )

                reordered_feat_loc = []
                for sys, df in self.dataframe.groupby(["system"]):
                    sys_data = raw_data[df.index[0] : df.index[-1] + 1]
                    np.save(
                        f"{self.analysis_results.filename}{feature}_{sys}.npy", sys_data
                    )
                    reordered_feat_loc.append(
                        [f"{self.analysis_results.filename}{feature}_{sys}.npy"]
                        * len(df)
                    )

                self.dataframe[feature] = np.concatenate(reordered_feat_loc)
                print(f"{feature} sorted.")
                del raw_data
                gc.collect()
                _ = [os.remove(location) for location in old_locations]

            self.sorted = True

            # update the analysis results
            self._init_dd_dataframe()

            if hasattr(self, "save_name"):
                print(f"Saving sorted results to {self.save_name}")
                self.save(self.save_name, overwrite=True)
        else:
            print("Already sorted")

    def add_analysis_result_from_data(self, data, feature_name, feature_info):
        if data.shape[0] != self.dataframe.shape[0]:
            print(
                f"Data shape {data.shape[0]} does not match the dataframe shape {self.dataframe.shape[0]}."
            )
            return

        if feature_name in self.analysis_list:
            print(f"{feature_name} already exists.")
            return

        feat_locs = []
        for sys, df in tqdm(
            self.dataframe.groupby(["system"]), total=self.dataframe.system.nunique()
        ):
            sys_data = data[df.index[0] : df.index[-1] + 1]
            np.save(
                f"{self.analysis_results.filename}{feature_name}_{sys}.npy", sys_data
            )
            feat_locs.append(
                [f"{self.analysis_results.filename}{feature_name}_{sys}.npy"] * len(df)
            )

        self.dataframe[f"{feature_name}"] = np.concatenate(feat_locs)
        self.analysis_list.append(f"{feature_name}")

        np.save(
            f"{self.analysis_results.filename}{feature_name}_feature_info.npy",
            feature_info,
        )

        if hasattr(self, "save_name"):
            self.save(self.save_name, overwrite=True)

    def remove_analysis(self, feature_name):
        """
        Remove an analysis from the dataframe.
        """
        self.analysis_list.remove(feature_name)
        self.analysis_results.dataframe = self.analysis_results.dataframe.drop(
            columns=[feature_name]
        )
        self.analysis_results.dd_dataframe = self.analysis_results.dd_dataframe.drop(
            columns=[feature_name]
        )
        # remove file
        file_paths = [
            location.replace(self.init_dir, self.working_dir)
            for location, df in self.dataframe.groupby(feature_name, sort=False)
        ]
        _ = [os.remove(file_path) for file_path in file_paths]
        self.dataframe = self.dataframe.drop(columns=[feature_name])

    def transform_to_logistic(self, feature_name, logistic):
        raw_data = np.concatenate(
            [
                np.load(
                    location.replace(self.init_dir, self.working_dir),
                    allow_pickle=True,
                )
                for location, df in self.dataframe.groupby(feature_name, sort=False)
            ]
        )

        #        if raw_data.shape[1] == 1:
        #            raw_data = np.hstack(raw_data).T

        scaler = preprocessing.MinMaxScaler(feature_range=(-logistic, logistic))

        scaled_data = scaler.fit_transform(raw_data)
        log_data = 1 / (1 + np.exp(-scaled_data))

        feat_locs = []
        for sys, df in tqdm(
            self.dataframe.groupby(["system"]), total=self.dataframe.system.nunique()
        ):
            sys_data = log_data[df.index[0] : df.index[-1] + 1]
            np.save(
                f"{self.analysis_results.filename}{feature_name}_log{logistic}_{sys}.npy",
                sys_data,
            )
            feat_locs.append(
                [
                    f"{self.analysis_results.filename}{feature_name}_log{logistic}_{sys}.npy"
                ]
                * len(df)
            )

        self.dataframe[f"{feature_name}_log{logistic}"] = np.concatenate(feat_locs)
        self.analysis_list.append(f"{feature_name}_log{logistic}")
        # TODO rename features
        shutil.copyfile(
            f"{self.analysis_results.filename}{feature_name}_feature_info.npy",
            f"{self.analysis_results.filename}{feature_name}_log{logistic}_feature_info.npy",
        )
        print("Finish transforming to logistic.")
        del raw_data
        gc.collect()

        if hasattr(self, "save_name"):
            self.save(self.save_name, overwrite=True)

    def transform_to_logistic_with_minmax(
        self, feature_name, logistic, min_arr, max_arr
    ):
        raw_data = np.concatenate(
            [
                np.load(
                    location.replace(self.init_dir, self.working_dir),
                    allow_pickle=True,
                )
                for location, df in self.dataframe.groupby(feature_name, sort=False)
            ]
        )
        scaled_data = (raw_data - min_arr) / (max_arr - min_arr)
        scaled_data = scaled_data * (2 * logistic) - logistic
        log_data = 1 / (1 + np.exp(-scaled_data))

        feat_locs = []
        for sys, df in tqdm(
            self.dataframe.groupby(["system"]), total=self.dataframe.system.nunique()
        ):
            sys_data = log_data[df.index[0] : df.index[-1] + 1]
            np.save(
                f"{self.analysis_results.filename}{feature_name}_logminmax{logistic}_{sys}.npy",
                sys_data,
            )
            feat_locs.append(
                [
                    f"{self.analysis_results.filename}{feature_name}_logminmax{logistic}_{sys}.npy"
                ]
                * len(df)
            )

        self.dataframe[f"{feature_name}_logminmax{logistic}"] = np.concatenate(
            feat_locs
        )
        self.analysis_list.append(f"{feature_name}_logminmax{logistic}")
        # TODO rename features
        shutil.copyfile(
            f"{self.analysis_results.filename}{feature_name}_feature_info.npy",
            f"{self.analysis_results.filename}{feature_name}_logminmax{logistic}_feature_info.npy",
        )
        print("Finish transforming to logistic.")
        del raw_data
        gc.collect()

        if hasattr(self, "save_name"):
            self.save(self.save_name, overwrite=True)

    def transform_to_reciprocal(self, feature_name):
        raw_data = np.concatenate(
            [
                np.load(
                    location.replace(self.init_dir, self.working_dir),
                    allow_pickle=True,
                )
                for location, df in self.dataframe.groupby(feature_name, sort=False)
            ]
        )

        #        if raw_data.shape[1] == 1:
        #            raw_data = np.hstack(raw_data).T

        _ = np.reciprocal(
            raw_data.astype(np.float64), out=raw_data, where=raw_data != 0
        )

        feat_locs = []
        for sys, df in tqdm(
            self.dataframe.groupby(["system"]), total=self.dataframe.system.nunique()
        ):
            sys_data = raw_data[df.index[0] : df.index[-1] + 1]
            np.save(
                f"{self.analysis_results.filename}{feature_name}_reciprocal_{sys}.npy",
                sys_data,
            )
            feat_locs.append(
                [f"{self.analysis_results.filename}{feature_name}_reciprocal_{sys}.npy"]
                * len(df)
            )

        self.dataframe[f"{feature_name}_reciprocal"] = np.concatenate(feat_locs)
        self.analysis_list.append(f"{feature_name}_reciprocal")

        # TODO rename features
        shutil.copyfile(
            f"{self.analysis_results.filename}{feature_name}_feature_info.npy",
            f"{self.analysis_results.filename}{feature_name}_reciprocal_feature_info.npy",
        )
        print("Finish transforming to reciprocal.")
        del raw_data
        gc.collect()

        if hasattr(self, "save_name"):
            self.save(self.save_name, overwrite=True)

    @classmethod
    def load_dataframe(cls, filename) -> "MDDataFrame":
        """
        Load the dataframe from a pickle file.

        Parameters
        ----------
        filename: str, optional
            The name of the pickle file.
        """
        if os.path.isfile(filename):
            print(f"Loading {filename}")
            with open(filename, "rb") as f:
                md_data = pickle.load(f)
        elif os.path.isfile(f"{filename}.pickle"):
            print(f"Loading {filename}.pickle")
            with open(f"{filename}.pickle", "rb") as f:
                md_data = pickle.load(f)
        else:
            print(f"Loading {filename}/{filename}_md_dataframe.pickle")
            with open(f"{filename}/{filename}_md_dataframe.pickle", "rb") as f:
                md_data = pickle.load(f)
        # check if the dataframe is a MDDataFrame
        if not isinstance(md_data, cls):
            raise TypeError("The loaded dataframe is not a MDDataFrame.")
        md_data.working_dir = os.getcwd() + "/"
        # if attribute init_dir not found
        if not hasattr(md_data, "init_dir"):
            md_data.init_dir = md_data.working_dir
        md_data.analysis_results.working_dir = md_data.filename
        return md_data

    @property
    def filename(self):
        """
        The saving location of all the pickled files.
        """
        return os.path.abspath(self.working_dir + self.dataframe_name) + "/"
