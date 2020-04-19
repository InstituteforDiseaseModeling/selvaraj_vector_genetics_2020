import os

import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser


class GMEstablishmentAnalyzer(BaseAnalyzer):

    def __init__(self, exp_name, sweep_variables, working_dir='.'):
        super(GMEstablishmentAnalyzer, self).__init__(working_dir=working_dir,
                                                      filenames=
                                                      ['output/ReportVectorGenetics_gambiae_Female_SPECIFIC_GENOME.csv']
                                                      )

        self.exp_name = exp_name
        self.sweep_variables = sweep_variables
        self.output_fname = os.path.join(self.working_dir, "%s_establishment_rates.csv" % self.exp_name)
        self.output_fname_full = os.path.join(self.working_dir, "%s_establishment_rates_full.csv" % self.exp_name)

    def select_simulation_data(self, data, simulation):

        simdata = data['output/ReportVectorGenetics_gambiae_Female_SPECIFIC_GENOME.csv']
        genomes = simdata['Genome'].unique()
        simdata = simdata.pivot_table('VectorPopulation', ['Time', 'NodeID'], 'Genome').reset_index()
        simdata['sum'] = [0] * len(simdata)
        for genome in genomes:
            simdata['sum'] += simdata[genome]
        for genome in genomes:
            simdata[genome] = simdata[genome] / simdata['sum']

        simdata = simdata[['Time', 'X-a0:X-a0',
                           'X-a0:X-a1', 'X-a1:X-a1']]

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            else:
                simdata[sweep_var] = 0
        return simdata

    def finalize(self, all_data):
        data_sets_per_experiment = {}

        for simulation, associated_data in all_data.items():
            experiment_name = simulation.experiment.exp_name
            if experiment_name not in data_sets_per_experiment:
                data_sets_per_experiment[experiment_name] = []

            data_sets_per_experiment[experiment_name].append(associated_data)

        for experiment_name, data_sets in data_sets_per_experiment.items():
            d = pd.concat(data_sets).reset_index(drop=True)
            d.to_csv(self.output_fname_full)
            self.sweep_variables.remove('Run_Number')
            columns = ['X-a0:X-a0', 'X-a0:X-a1', 'X-a1:X-a1']
            d_mean = d.groupby(self.sweep_variables + ['Time'])['X-a0:X-a0', 'X-a0:X-a1',
                                                                'X-a1:X-a1'].apply(np.mean).reset_index()
            d_std = d.groupby(self.sweep_variables + ['Time'])['X-a0:X-a0', 'X-a0:X-a1',
                                                               'X-a1:X-a1'].apply(np.std).reset_index()
            for channel in columns:
                d_mean[channel + '_std'] = d_std[channel]
            d_mean.to_csv(self.output_fname, index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    projectdir = '/_data'
    out_dir = os.path.join(projectdir)

    exp_id = ''

    experiments = {
        "spatial_sim": exp_id
    }

    sweep_vars = ["Run_Number", "Coverage", "Start",
                  "Mutation_Rate1", "Mutation_Rate2", "Label"
                  ]

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                GMEstablishmentAnalyzer(exp_name='GM_release',
                                                        working_dir=out_dir,
                                                        sweep_variables=sweep_vars
                                                        )
                            ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()
