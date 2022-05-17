"""\
========
Analysis
========
The :class:`~ENPMDA.analysis.DaskChunkMdanalysis` class is
the base class to define multi-frame parallel analysis for
the MD trajectories. It functions as a wrapper of the
MDAnalysis analysis functions to map the analysis to the Dask dataframe.
This class takes care of loading the right universe and dumping the
results as a `npy` file to avoid huge memory footprint
and dask scheduler clogging. 

To define a new analysis, ``DaskChunkMdanalysis`` needs to be subclassed.
`run_analysis` need to be defined. ``name`` will be the feature name
appending to the dataframe. In default, only protein universe file
will be used to run analysis. It can be overridden by defining
``universe_file=system``. Some information of the feature can also be
provided and stored in ``_feature_info``::

    from ENPMDA.analysis import DaskChunkMdanalysis
    class NewAnalysis(DaskChunkMdanalysis):
        name = 'newanalysis'
        universe_file = 'system'

        def run_analysis(self, universe, start, stop, step):
            self._feature_info = ['new_analysis']
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
#        self.ref_info = {}
#        self.u_ref = u_ref

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
        kwargs['ref_df'] = self.dataframe.iloc[0]
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

#        self.ref_info[item_name] = analysis_function(
#            **kwargs).run_analysis(self.u_ref, 0, 3, 1)

    def save_results(self):
        for item, df in self.items():
            if isinstance(df, dask.dataframe.core.DataFrame):
                df.to_csv(
                    self.working_dir +
                    '/analysis_results/' +
                    str(item) +
                    '-*.csv')

    def filter_result(self, column, filter_threshold):
        """
        filter results based on threshold
        not used
        """
        filter_index = []
        for datafile in set(self[column].iloc[:, 1]):
            data = np.load(datafile, allow_pickle=True)
            filter_index.extend(
                np.where(
                    np.any(
                        data < filter_threshold,
                        axis=0)))
        filter_index = list(set(filter_index[0]))

        for datafile in set(self[column].iloc[:, 1]):
            data = np.load(datafile, allow_pickle=True)
            np.save(datafile +
                    '.filter_' +
                    str(filter_threshold) +
                    '.npy', data[:, filter_index])

        self[column + '.filter_' +
             str(filter_threshold)] = self[column].copy(deep=True)
        self[column +
             '.filter_' +
             str(filter_threshold)].iloc[:, 1] = self[column].iloc[:, 1].apply(lambda x: x +
                                                                               '.filter_' +
                                                                               str(filter_threshold) +
                                                                               '.npy')

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
        self._feature_info = []

        # sanity check and get feature info
        if self.universe_file == 'protein':
            universe = pickle.load(open(self.ref_df.universe_protein, "rb"))
        else:
            universe = pickle.load(open(self.ref_df.universe_system, "rb"))
        self.set_feature_info(universe)
        check_result = np.asarray(self.run_analysis(universe, 0, 2, 1))
        check_result = check_result.reshape(check_result.shape[0], -1)
        if check_result.shape[1] != len(self._feature_info):
            raise ValueError(f'The shape of feature info ({len(self._feature_info)}) '
            f'does not match the shape of analysis {check_result.shape[1]}.')
        if check_result.shape[0] != 2:
            raise ValueError('The len of the result'
            'does not match the number of the frame.'
            'Hint: run_analysis should return a list'
            'with the shape as the trajectory frames.')

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
        np.save(self.feature_info_loc, self._feature_info)
        return n_result * [self.location]

#        return list(itertools.chain.from_iterable(result))

    def set_feature_info(self, universe):
        """
        This function is used to set the feature information
        """
        self._feature_info = []
        raise NotImplementedError('Only for inheritance.')

    def run_analysis(self, universe, start, stop, step):
        """
        The function to be overwritten by the analysis class.
        """
        raise NotImplementedError('Only for inheritance.')

    @property
    def feature_info(self):
        #  contain info of n features
        return self._feature_info

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