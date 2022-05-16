import os
import dask
import numpy as np
import uuid
import pickle
import itertools


class AnalysisResult(dict):
    def __init__(self,
                 md_data,
                 working_dir,
                 timestamp):
        """
        md_data: dask.dataframe.core.DataFrame
        store all the data

        reference structures
        """
        super().__init__()

        self.md_data = md_data
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
        item_name = analysis_function().name
        meta = dict(self.md_data.dtypes)
        meta[item_name] = "f8"

        # add analysis to dataframe
        self[item_name] = self.md_data.map_partitions(lambda df:
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
    def filename(self, timestamp=None):
        return self.working_dir + '/analysis_results/' + self.timestamp + '/'


class DaskChunkMdanalysis(object):
    name = 'analysis'
    universe_file = 'protein'

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

        # get unique uuid for each partition
        self._partition = uuid.uuid4().hex
        self._feature_info = []

    def __call__(self, df):
        return self.run_df(df)

    def run_df(self, df):
        result = []
        for system, df_sys in df.groupby(['system']):
            # if system information has to be used set `universe_file =
            # 'system'`
            if self.universe_file == 'protein':
                universe = pickle.load(open(df_sys.universe.iloc[0], "rb"))
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

    def run_analysis(self, universe, start, stop, step):
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
