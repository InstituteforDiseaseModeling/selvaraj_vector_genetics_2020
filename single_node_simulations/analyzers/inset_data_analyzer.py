import os
import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser


class InsetAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, report_names=["InsetChart"], channels=None, sweep_variables=None, working_dir="."):
        super(InsetAnalyzer, self).__init__(working_dir=working_dir, filenames=["output/{name}.json".format(name=name)
                                                                                for name in report_names]
                                            )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.channels = channels or ['Annual EIR']
        self.reports = report_names
        self.expt_name = expt_name
        self.output_fname = os.path.join(self.working_dir, "%s_inset_data.csv" % self.expt_name)
        self.output_fname_full = os.path.join(self.working_dir, "%s_inset_data_full.csv" % self.expt_name)

    def select_simulation_data(self, data, simulation):
        simdata = []

        for report in self.reports:

            datatemp = data["output/{name}.json".format(name=report)]

            prevalence = datatemp['Channels']['Blood Smear Parasite Prevalence']['Data']
            eir = datatemp['Channels']['Daily EIR']['Data']
            clinical_cases = datatemp['Channels']['New Clinical Cases']['Data']
            true_prevalence = datatemp['Channels']['True Prevalence']['Data']

            df = pd.DataFrame(list(zip(prevalence, eir, clinical_cases, true_prevalence)),
                              columns=['RDT Prevalence', 'Annual EIR', 'Clinical Cases', 'True Prevalence'])
            df['Time'] = [i for i in range(len(df))]
            simdata.append(df)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            else:
                simdata[sweep_var] = 0
        return simdata

    def finalize(self, all_data):

        data_sets_per_experiment = {}

        # for simulation, associated_data in all_data.items():
        for simulation, associated_data in all_data.items():
            experiment_name = simulation.experiment.exp_name
            if experiment_name not in data_sets_per_experiment:
                data_sets_per_experiment[experiment_name] = []

            data_sets_per_experiment[experiment_name].append(associated_data)

        for experiment_name, data_sets in data_sets_per_experiment.items():
            d = pd.concat(data_sets).reset_index(drop=True)
            d.to_csv(self.output_fname_full)
            self.sweep_variables.remove('Run_Number')
            d_mean = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(np.mean).reset_index()
            d_std = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(np.std).reset_index()
            for channel in ['RDT Prevalence', 'Annual EIR', 'Clinical Cases', 'True Prevalence']:
                d_mean[channel + '_std'] = d_std[channel]
            d_mean.to_csv(self.output_fname, index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    projectdir = '/_data'
    out_dir = os.path.join(projectdir)

    exp_id = ''

    experiments = {
        "single_node_simulations": exp_id
    }

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                InsetAnalyzer(working_dir=out_dir,
                                                expt_name=expt_name,
                                                report_names=['InsetChart'],
                                                sweep_variables=["Run_Number", "Coverage", "Start",
                                                                 "Mutation_Rate1", "Mutation_Rate2",
                                                                 ])
                                       ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()

