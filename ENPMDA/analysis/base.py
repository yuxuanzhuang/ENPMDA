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
`run_analysis` need to be defined. ``name`` will be the feature name
appending to the dataframe. In default, only protein universe file
will be used to run analysis. It can be overridden by defining
``universe_file=system``. Some information of the feature can also be
provided and stored in ``_feature_info``::

    from ENPMDA.analysis import DaskChunkMdanalysis
    class NewAnalysis(DaskChunkMdanalysis):
        name = 'new_analysis'
        universe_file = 'protein'

        def get_feature_info(self, universe):
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


class AnalysisResult(dict):
    """
    This class is used internally to store the results of the analysis.
    """

    def __init__(self,
                 dd_dataframe,
                 dataframe,
                 working_dir,
                 timestamp):
        super().__init__()

        self.dd_dataframe = dd_dataframe
        self.dataframe = dataframe
        self.working_dir = working_dir

        self.timestamp = timestamp

        os.makedirs(self.filename, exist_ok=True)

    def compute(self, item=None):
        if item is None:
            for item, df in self.items():
                if isinstance(df, dask.dataframe.core.DataFrame):
                    self[item] = df.compute()[['system', item]]
        elif item in self.keys():
            if isinstance(df, dask.dataframe.core.DataFrame):
                self[item] = df.compute()[['system', item]]
        else:
            raise ValueError(item + ' not in dict')

    def add_column_to_results(self, analysis_function, **kwargs):
        kwargs['filename'] = self.filename

        # sanity check and get feature info
        check_analysis_function = analysis_function(**kwargs)
        if check_analysis_function.universe_file == 'protein':
            universe = pickle.load(open(self.dataframe.iloc[0].universe_protein, "rb"))
        else:
            universe = pickle.load(open(self.dataframe.iloc[0].universe_system, "rb"))
        
        feature_info = check_analysis_function.set_feature_info(universe)
        check_result = np.asarray(check_analysis_function.run_analysis(universe, 0, 2, 1))
        check_result = check_result.reshape(check_result.shape[0], -1)
        if check_result.shape[1] != len(feature_info):
            raise ValueError(f'The shape of feature info ({len(feature_info)}) '
            f'does not match the shape of analysis ({check_result.shape[1]}).')
        if check_result.shape[0] != 2:
            raise ValueError('The len of the result'
            'does not match the number of the frame.'
            'Hint: run_analysis should return a list'
            'with the shape as the trajectory frames.')
        np.save(check_analysis_function.feature_info_loc, feature_info)

        item_name = analysis_function.name
        meta = dict(self.dd_dataframe.dtypes)
        meta[item_name] = "f8"

        # add analysis to dataframe
        self[item_name] = self.dd_dataframe.map_partitions(lambda df:
                                                      df.assign(
                                                          **{item_name: analysis_function(**kwargs)(df)}),
                                                      meta=meta,
                                                      token=item_name
                                                      ).persist()


    def append_to_dataframe(self, dataframe):
        for item in self.keys():
            dataframe[item] = self[item].iloc[:, 1]

    @property
    def filename(self):
        return self.working_dir + '/analysis_results/' + self.timestamp + '/'


class DaskChunkMdanalysis(object):
    """
    This class is the base class for all analysis classes.
    The analysis results will be dumped as a `npy` file with
    unique uuid for each partition.
    """

    name = 'analysis'
    universe_file = 'protein'

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

        # get unique uuid for each partition
        self._partition = uuid.uuid4().hex
        self.feature_info = []

    def __call__(self, df):
        return self.run_df(df)

    def run_df(self, df):
        result = []

        # run analysis
        for system, df_sys in df.groupby(['system']):
            # if system information has to be used set `universe_file =
            # 'system'`
            if self.universe_file == 'protein':
                universe = pickle.load(open(df_sys.universe_protein.iloc[0], "rb"))
            else:
                universe = pickle.load(
                    open(df_sys.universe_system.iloc[0], "rb"))
            start = df_sys.frame.iloc[0]
            stop = df_sys.frame.iloc[-1] + 1
            step = df_sys.stride.iloc[0]
            start, stop, step = universe.trajectory.check_slice_indices(
                start, stop, step)

            result.append(self.run_analysis(universe, start, stop, step))
            del universe

        result_ordered = list(itertools.chain.from_iterable(result))
        # store results and only return the location of the files
        np.save(self.location, result_ordered)
        n_result = len(result_ordered)
        del result_ordered
        return n_result * [self.location]

#        return list(itertools.chain.from_iterable(result))

    def set_feature_info(self, universe):
        """
        This function is used to set the feature information.
        Shold return a list of features.
        """
        raise NotImplementedError('Only for inheritance.')

    def run_analysis(self, universe, start, stop, step):
        """
        The function to be overwritten by the analysis class.
        """
        raise NotImplementedError('Only for inheritance.')

    @property
    def name(self):
        return self._name

    @property
    def partition(self):
        return str(self._partition)

    @property
    def location(self):
        return self.filename + self.name + '_' + self.partition + '.npy'

    @property
    def feature_info_loc(self):
        return self.filename + self.name + '_feature_info.npy'