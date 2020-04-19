import os

import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser
from spatial_output_dataframe import construct_spatial_output_df
import glob

projectdir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects', 'Vector_genetics',
                          'Data', 'simulation_data', 'Establishment_Burkina')


class SpatialAnalyzer(BaseAnalyzer):

    def __init__(self, exp_name, spatial_channels, sweep_variables, working_dir='.'):
        super(SpatialAnalyzer, self).__init__(working_dir=working_dir,
                                             filenames=['output/SpatialReportMalariaFiltered_%s.bin' % x for x in spatial_channels]
                                           )

        self.exp_name = exp_name
        self.sweep_variables = sweep_variables
        self.spatial_channels = spatial_channels
        self.output_fname = os.path.join(self.working_dir, "%s_spatial_data.csv" % self.exp_name)

    def select_simulation_data(self, data, simulation):

        simdata = construct_spatial_output_df(data['output/SpatialReportMalariaFiltered_%s.bin' % self.spatial_channels[0]], self.spatial_channels[0])
        if len(self.spatial_channels) > 1:
            for ch in self.spatial_channels[1:]:
                simdata = pd.merge(left=simdata,
                                   right=construct_spatial_output_df(data['output/SpatialReportMalariaFiltered_%s.bin' % ch], ch),
                                   on=['time', 'node'])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            else:
                simdata[sweep_var] = 0
        return simdata
        # # ensure data directory
        # # because we are running in threads, this could happen in multiple threads
        # try:
        #     os.makedirs('intermediate_data', exist_ok=True)
        # except FileExistsError:
        #     pass
        # out_path = f'intermediate_data/{simulation.id}.pickle'
        #
        # simdata.to_pickle(out_path)
        # # return out_path
        # return out_path

    def finalize(self, all_data):
        # load all the simulations
        # all_data_dict = dict()
        # for simulation, path in all_data.items():
        #     all_data_dict[simulation] = pd.read_pickle(path)

        data_sets_per_experiment = {}

        # for simulation, associated_data in all_data.items():
        for simulation, associated_data in all_data.items():
            experiment_name = simulation.experiment.exp_name
            if experiment_name not in data_sets_per_experiment:
                data_sets_per_experiment[experiment_name] = []

            data_sets_per_experiment[experiment_name].append(associated_data)

        for experiment_name, data_sets in data_sets_per_experiment.items():
            d = pd.concat(data_sets).reset_index(drop=True)
            self.sweep_variables.remove('Run_Number')
            d = d.groupby(self.sweep_variables + ['time', 'node'])[self.spatial_channels].apply(np.mean).reset_index()
            d.to_csv(self.output_fname, index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    out_dir = os.path.join(projectdir)

    exp_id = ''

    experiments = {
                   "spatial_sim": exp_id
                   }
    sweep_vars = ['Run_Number', 'Release_Number', 'Copy_To_Likelihood',
                  'ITN_Coverage', 'Transmission_To_Human',
                  'Start_Day', 'Num_Nodes']

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                SpatialAnalyzer(exp_name='GM_release',
                                                working_dir=out_dir,
                                                sweep_variables=sweep_vars,
                                                spatial_channels=['Prevalence', 'Daily_EIR', 'Adult_Vectors'],
                                                )
                                       ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()