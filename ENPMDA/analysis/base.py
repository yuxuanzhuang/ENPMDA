"""\
========
Analysis
========
The :class:`~ENPMDA.analysis.base.DaskChunkMdanalysis` class is
the base class to define multi-frame parallel analysis for
the MD trajectories. It functions as a wrapper of the
MDAnalysis analysis functions to map the analysis to the Dask dataframe.
This class takes care of loading the right universe and dumping the
results as a `npy` file to avoid huge memory footprint
and dask scheduler clogging. 

To define a new analysis,
:class:`~ENPMDA.analysis.base.DaskChunkMdanalysis`
needs to be subclassed.
``set_feature_info(self, universe)`` and ``run_analysis(self, universe, start, stop, step)`` need to be defined. ``set_feature_info`` should return
a list of feature name e.g. the name of each torsion angle. ``run_analysis`` should return a list of analysis results.

``name`` will be the feature name appending to the dataframe.
In default, only protein universe file will be used to run analysis.
It can be overridden by defining ``universe_file=system``::

    from ENPMDA.analysis import DaskChunkMdanalysis
    class NewAnalysis(DaskChunkMdanalysis):
        name = 'new_analysis'
        universe_file = 'protein'

        def set_feature_info(self, universe):
            return ['some_info']

        def run_analysis(self, universe, start, stop, step):
            result = []
            for ts in universe.trajectory[start:stop:step]:
                result.append(some_analysis(universe.atoms))
            return result

Classes
=======
.. autoclass:: DaskChunkMdanalysis
   :members:
"""
import os
import dask
import numpy as np
import uuid
import pickle
import itertools
import warnings
import dask.dataframe as dd


class AnalysisResult(dict):
    """
    This class is used internally to store the results of the analysis.
    """

    def __init__(self, dd_dataframe, dataframe, working_dir, timestamp):
        super().__init__()

        self.dd_dataframe = dd_dataframe
        self.dataframe = dataframe
        self.working_dir = working_dir

        self.timestamp = timestamp

        os.makedirs(self.filename, exist_ok=True)

    def compute(self, item=None):
        if item is None:
            for it, df in self.items():
                if isinstance(df, dd.DataFrame):
                    self[it] = df.compute()[["system", it]]
        elif item in self.keys():
            if isinstance(df, dd.DataFrame):
                self[item] = df.compute()[["system", item]]
        else:
            raise ValueError(item + " not in dict")

    def add_column_to_results(self, analysis_function, **kwargs):
        # sanity check and get feature info
        check_analysis_function = analysis_function(filename=self.filename, **kwargs)
        if check_analysis_function.universe_file == "protein":
            with open(self.dataframe.iloc[0].universe_protein, "rb") as f:
                universe = pickle.load(f)
        else:
            with open(self.dataframe.iloc[0].universe_system, "rb") as f:
                universe = pickle.load(f)

        feature_info = check_analysis_function.set_feature_info(universe)
        if check_analysis_function.output == "array":
            check_result = np.asarray(
                check_analysis_function.run_analysis(universe, 0, 2, 1)
            )
        else:
            check_result = np.asarray(
                check_analysis_function.run_analysis(universe, 0, 2, 1), dtype=object
            )
        if check_result.ndim > 2:
            warnings.warn(
                f"The result of the analysis function is not 2D."
                f"Make sure the shape of feature info ({len(feature_info)}) "
                f"does match the shape of analysis ({check_result.shape[1]}).",
                stacklevel=2,
            )
            check_result = check_result.reshape(
                check_result.shape[0], len(feature_info), -1
            )
        else:
            check_result = check_result.reshape(check_result.shape[0], -1)
            if check_result.shape[1] != len(feature_info):
                raise ValueError(
                    f"The shape of feature info ({len(feature_info)}) "
                    f"does not match the shape of analysis ({check_result.shape[1]})."
                )
        if check_result.shape[0] != 2:
            raise ValueError(
                "The len of the result"
                "does not match the number of the frame."
                "Hint: run_analysis should return a list"
                "with the shape as the trajectory frames."
            )
        np.save(check_analysis_function.feature_info_loc, feature_info)
        check_analysis_function._feature_info = feature_info

        item_name = analysis_function.name
        meta = dict(self.dd_dataframe.dtypes)
        meta[item_name] = "f8"

        # add analysis to dataframe
        kwargs.update(check_analysis_function.__dict__)
        self[item_name] = self.dd_dataframe.map_partitions(
            lambda df: df.assign(**{item_name: analysis_function(**kwargs)(df)}),
            meta=meta,
            token=item_name,
        ).persist()

    def append_to_dataframe(self, dataframe):
        for item in self.keys():
            dataframe[item] = self[item].iloc[:, -1]

    @property
    def filename(self):
        return self.working_dir + "/analysis_results/" + self.timestamp + "/"


class DaskChunkMdanalysis(object):
    """
    This class is the base class for all analysis classes.
    The analysis results will be dumped as a `npy` file with
    unique uuid for each partition.
    """

    name = "analysis"
    universe_file = "protein"
    output = "array"

    def __init__(self, filename, **kwargs):
        self._feature_info = []
        self.filename = filename
        self.__dict__.update(kwargs)

        # get unique uuid for each partition
        self._partition = uuid.uuid4().hex

    def __call__(self, df):
        return self.run_df(df)

    def run_df(self, df):
        result = []

        # run analysis
        for system, df_sys in df.groupby("system", sort=False):
            # if system information has to be used set `universe_file =
            # 'system'`
            if self.universe_file == "protein":
                with open(df_sys.universe_protein.iloc[0], "rb") as f:
                    universe = pickle.load(f)
            else:
                with open(df_sys.universe_system.iloc[0], "rb") as f:
                    universe = pickle.load(f)
            start = df_sys.frame.iloc[0]
            stop = df_sys.frame.iloc[-1] + 1
            step = df_sys.stride.iloc[0]
            start, stop, step = universe.trajectory.check_slice_indices(
                start, stop, step
            )

            result.append(self.run_analysis(universe, start, stop, step))
            del universe

        result_ordered = list(itertools.chain.from_iterable(result))

        # make sure the result is 2D with the shape of the trajectory frames and
        # the number of features

        if self.output == "array":
            result_ordered_arr = np.asarray(result_ordered)
        else:
            result_ordered_arr = np.empty(len(result_ordered), dtype=object)
            result_ordered_arr[...] = [
                np.asarray(x, dtype=object).reshape(len(self.feature_info), -1)
                for x in result_ordered
            ]
        #        result_ordered = result_ordered.reshape(result_ordered.shape[0], len(self.feature_info), -1)

        # store results and only return the location of the files
        np.save(self.location, result_ordered_arr)
        n_result = len(result_ordered_arr)
        del result_ordered_arr

        return n_result * [self.location]

        # not returning the actual array
        # return list(itertools.chain.from_iterable(result))

    def set_feature_info(self, universe):
        """
        This function is used to set the feature information.
        Shold return a list of features.
        """
        raise NotImplementedError("Only for inheritance.")

    def run_analysis(self, universe, start, stop, step):
        """
        The function to be overwritten by the analysis class.
        """
        raise NotImplementedError("Only for inheritance.")

    @classmethod
    def test_on_universe(cls, universe, start=0, stop=2, step=1):
        """
        This function is used to test the analysis function on a universe.
        """
        return cls(filename='test').run_analysis(universe,
                                                 start,
                                                 stop,
                                                 step)

    @property
    def partition(self):
        return str(self._partition)

    @property
    def location(self):
        return self.filename + self.name + "_" + self.partition + ".npy"

    @property
    def feature_info_loc(self):
        return self.filename + self.name + "_feature_info.npy"

    @property
    def feature_info(self):
        return self._feature_info
