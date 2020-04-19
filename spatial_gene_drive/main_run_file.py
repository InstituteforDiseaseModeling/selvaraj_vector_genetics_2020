from dtk.utils.reports import BaseVectorGeneticsReport
from malaria.reports.MalariaReport import add_filtered_spatial_report
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from dtk.vector.species import set_species_genes

from spatial_gene_drive.configure_interventions import *
from spatial_gene_drive.helper_functions import *
from malaria.interventions.health_seeking import add_health_seeking


if __name__ == "__main__":

    exp_name = "spatial_vector_genetics"

    geography = 'Burkina Faso'
    prefix = "vector_genetics"

    migration_mul = 100
    num_cores = 5
    num_years = 6
    num_seeds = 50

    ############################ SETUP CONFIG FILE and SERIALIZED FILE ##################3
    cb = configure_VC_GM_intervention_system(prefix, geography,
                                             num_cores=num_cores, num_years=num_years,
                                             migration_mul=migration_mul)

    serialized_file_list = ['', '']

    cb.update_params({
        "Report_Event_Recorder": 0,
        "Listed_Events": ["Immigrating", "Emigrating"],
        "Report_Event_Recorder_Events": ["Immigrating", "Emigrating"],
        "Report_Event_Recorder_Ignore_Events_In_List": 0,
        'Serialized_Population_Filenames': serialized_file_list,
        'logLevel_JsonConfigurable': 'ERROR',
        'logLevel_VectorHabitat': 'ERROR',
        'logLevel_Simulation': 'ERROR',
        'logLevel_StandardEventCoordinator': 'ERROR',
        'logLevel_LarvalHabitatMultiplier': 'ERROR',
        'logLevel_SimulationEventContext': 'ERROR',
        'logLevel_SusceptibilityMalaria': 'ERROR'
    })

    #################################### VECTOR GENETICS ##############################

    # Add genes
    genes = {'gambiae': [{
        "Alleles": {"a0": 1.0, "a1": 0.0},
        "Mutations": {"a1:a0": 0.05}
    }]
    }
    set_species_genes(cb, genes)

    release_numbers = [1000]
    net_coverages = [0.8]
    copy_to_likelihoods = np.arange(1.0, 1.01, 0.1).tolist()
    release_nodes = [6]
    transmission_probs = [0.3]
    start_days = [180]
    VC_and_GM = [
            [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
             ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', ''),
             ModFn(add_nets, coverage=net_coverage, number=number, num_nodes=numnodes, start_day=start_day),
             ModFn(add_release, number=number, num_nodes=numnodes, start_day=start_day),
             ModFn(add_drivers, copy_to_likelihood=likelihood),
             ModFn(add_trait_modifiers, transmission_to_human=transmission_prob)
             ]
            for seed in range(num_seeds)
            for start_day in start_days
            for number in release_numbers
            for numnodes in release_nodes
            for net_coverage in net_coverages
            for likelihood in copy_to_likelihoods
            for transmission_prob in transmission_probs
        ]

    VC = [
            [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
             ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', ''),
             ModFn(add_nets, coverage=net_coverage, number=0, num_nodes=0, start_day=start_day),
             ModFn(add_drivers, copy_to_likelihood=likelihood),
             ModFn(add_trait_modifiers, transmission_to_human=transmission_prob)
             ]
            for seed in range(num_seeds)
            for start_day in [0]
            for number in [0]
            for numnodes in [0]
            for net_coverage in net_coverages
            for likelihood in [0.0]
            for transmission_prob in [1.0]
        ]

    GM = [
        [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
         ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', ''),
         ModFn(add_nets, coverage=net_coverage, number=number, num_nodes=numnodes, start_day=start_day),
         ModFn(add_release, number=number, num_nodes=numnodes, start_day=start_day),
         ModFn(add_drivers, copy_to_likelihood=likelihood),
         ModFn(add_trait_modifiers, transmission_to_human=transmission_prob)
         ]
        for seed in range(num_seeds)
        for start_day in start_days
        for number in release_numbers
        for numnodes in release_nodes
        for net_coverage in [0.0]
        for likelihood in copy_to_likelihoods
        for transmission_prob in transmission_probs
    ]

    No_interventions = [
        [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
         ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', ''),
         ModFn(add_nets, coverage=net_coverage, number=0, num_nodes=0, start_day=start_day),
         ModFn(add_drivers, copy_to_likelihood=likelihood),
         ModFn(add_trait_modifiers, transmission_to_human=transmission_prob)
         ]
        for seed in range(num_seeds)
        for start_day in [0]
        for number in [0]
        for numnodes in [0]
        for net_coverage in [0.0]
        for likelihood in [0.0]
        for transmission_prob in [1.0]
    ]

    builder = ModBuilder.from_list(VC_and_GM + VC + GM + No_interventions)
    # builder = ModBuilder.from_list(VC)

    ################################ INTERVENTIONS ###################################

    # Health seeking
    add_health_seeking(cb,
                       targets=[{'trigger': 'NewClinicalCase',
                                 'coverage': 0.5,
                                 'agemin': 0,
                                 'agemax': 100,
                                 'seek': 1,
                                 'rate': 0.3},
                                {'trigger': 'NewSevereCase',
                                 'coverage': 0.8,
                                 'agemin': 0,
                                 'agemax': 100,
                                 'seek': 1,
                                 'rate': 0.5}
                                ],
                       drug=['Artemether', 'Lumefantrine'],
                       start_day=0,
                       broadcast_event_name='Received_Treatment')

    #################################### REPORTERS #######################################
    # cb.add_reports(BaseVectorStatsReport(type='ReportVectorStats', stratify_by_species=1))
    cb.add_reports(BaseVectorGeneticsReport(type='ReportVectorGenetics',
                                            species='gambiae',
                                            gender='VECTOR_FEMALE',
                                            include_vector_state_columns=0,
                                            stratify_by='SPECIFIC_GENOME',
                                            combine_similar_genomes=1,
                                            specific_genome_combinations_for_stratification=[
                                                {
                                                    "Allele_Combination": [
                                                        ["X", "X"],
                                                        ["a1", "*"]
                                                    ]
                                                },
                                                {
                                                    "Allele_Combination": [
                                                        ["X", "X"],
                                                        ["a0", "a0"]
                                                    ]
                                                }
                                            ]
                                            ))
    add_filtered_spatial_report(cb, start=0, end=365 * num_years,
                                channels=['Population', 'Prevalence', 'New_Clinical_Cases',
                                          'Daily_EIR', 'Adult_Vectors'])

    ###################################################################################################################

    run_sim_args = {'config_builder': cb,
                    'exp_name': exp_name,
                    'exp_builder': builder}

    SetupParser.default_block = 'HPC'

    SetupParser.init('HPC')
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
