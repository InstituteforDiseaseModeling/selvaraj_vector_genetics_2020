from malaria.interventions.health_seeking import add_health_seeking

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from single_node_simulations.helper_functions import change_vector_params, update_config_params, add_reporters, \
    add_ITNs, add_insecticides

if __name__ == "__main__":
    dir = './'
    geography = 'Burkina Faso'
    prefix = 'single_node_simulations'
    exp_name = "insecticide_resistance_single_node"

    extra = 0
    num_years = 6 + extra
    num_seeds = 50

    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM',
                                        Simulation_Duration=int(365 * num_years))

    update_config_params(cb, direc=dir, geography=geography)
    add_insecticides(cb, pyrethroid_killing=0.125, carbamate_killing=0.25)

    serialized_file_list = ['', '']
    cb.update_params({
        'Serialized_Population_Filenames': serialized_file_list
    })

    ########################## VECTOR GENETICS ############################################

    three_year_gap = [
            [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
             ModFn(change_vector_params, species='gambiae', mutation_rate1=mutation_rate1,
                   mutation_rate2=mutation_rate2, serialization_seed=0),
             ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', ''),
             ModFn(add_ITNs, start=start, insecticide='pyrethroid', label='3year'),
             ModFn(add_ITNs, start=start+3*365, insecticide='carbamate', label='3year')
             ]
            for seed in range(num_seeds)
            for mutation_rate1 in [pow(10, x) for x in range(-4, -1)]
            for mutation_rate2 in [pow(10, x) for x in range(-4, -1)]
            for start in [180 + extra*365]
        ]

    no_rotation_three_year = [
        [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
         ModFn(change_vector_params, species='gambiae', mutation_rate1=mutation_rate1,
               mutation_rate2=mutation_rate2, serialization_seed=0),
         ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', ''),
         ModFn(add_ITNs, coverage=0.6, start=start, insecticide='pyrethroid', label='No_rotation_3year'),
         ModFn(add_ITNs, coverage=0.6, start=start + 3 * 365, insecticide='pyrethroid', label='No_rotation_3year')
         ]
        for seed in range(num_seeds)
        for mutation_rate1 in [0.01]
        for mutation_rate2 in [0]
        for start in [180 + extra*365]
    ]

    no_nets = [
        [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
         ModFn(change_vector_params, species='gambiae', mutation_rate1=mutation_rate1,
               mutation_rate2=mutation_rate2, serialization_seed=0),
         ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', ''),
         ModFn(add_ITNs, coverage=0, start=start, insecticide='pyrethroid', label='No_nets'),
         ]
        for seed in range(num_seeds)
        for mutation_rate1 in [0.01]
        for mutation_rate2 in [0]
        for start in [180 + extra*365]
    ]

    no_resistance = [
        [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
         ModFn(change_vector_params, species='gambiae', mutation_rate1=mutation_rate1,
               mutation_rate2=mutation_rate2, serialization_seed=0),
         ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', ''),
         ModFn(add_ITNs, coverage=0.6, start=start, insecticide='carbamate', label='No_resistance'),
         ModFn(add_ITNs, coverage=0.6, start=start + 3 * 365, insecticide='carbamate', label='No_resistance')
         ]
        for seed in range(num_seeds)
        for mutation_rate1 in [0.01]
        for mutation_rate2 in [0]
        for start in [180 + extra * 365]
    ]

    builder = ModBuilder.from_list(no_rotation_three_year + no_nets + no_resistance)

    ############################### REPORTERS ###############################################
    vector_genetics_report = 1
    vector_stats_report = 0
    malaria_summary_report = 0

    add_reporters(cb, vector_genetics_report=vector_genetics_report)

    ############################## HEALTH SEEKING ###########################################
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

    ##########################################################################################

    run_sim_args = {'config_builder': cb,
                    'exp_name': exp_name,
                    'exp_builder': builder}

    SetupParser.default_block = 'HPC'

    SetupParser.init('HPC')
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
