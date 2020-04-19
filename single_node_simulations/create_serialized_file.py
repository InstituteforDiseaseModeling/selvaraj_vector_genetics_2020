from malaria.interventions.health_seeking import add_health_seeking
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from dtk.generic.serialization import add_SerializationTimesteps

from single_node_simulations.helper_functions import change_vector_params, update_config_params, add_reporters


if __name__ == "__main__":
    dir = './'
    geography = 'Burkina Faso'
    prefix = 'single_node_simulations'
    exp = 8.60
    exp_name = "insecticide_resistance_single_node_serialization"

    num_years = 40
    num_seeds = 1

    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM',
                                        Simulation_Duration=int(365 * num_years))

    update_config_params(cb, direc=dir, geography=geography)
    add_SerializationTimesteps(cb, [num_years*365], end_at_final=True)

    ########################## VECTOR GENETICS ####################################################

    builder = ModBuilder.from_list(
        [
            [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
             ModFn(change_vector_params, species='gambiae',
                   mutation_rate1=mutation_rate, mutation_rate2=mutation_rate, serialization_seed=seed)
             ]
            for seed in range(num_seeds)
            for mutation_rate in [0]  # pow(10, x) for x in range(-5, 0)
        ]
    )

    ############################### REPORTERS ########################################
    vector_genetics_report = 1
    vector_stats_report = 0
    malaria_summary_report = 0

    add_reporters(cb, vector_genetics_report=vector_genetics_report)

    ################################ INTERVENTIONS ###################################
    #
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
                       start_day=(num_years - 10) * 365,
                       broadcast_event_name='Received_Treatment')

    #####################################################################################

    run_sim_args = {'config_builder': cb,
                    'exp_name': exp_name,
                    'exp_builder': builder}

    SetupParser.default_block = 'HPC'

    SetupParser.init('HPC')
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
