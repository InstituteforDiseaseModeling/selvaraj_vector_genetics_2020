import os

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from dtk.generic.serialization import add_SerializationTimesteps
from dtk.vector.species import set_species, set_species_param, set_species_genes

from malaria.interventions.health_seeking import add_health_seeking


if __name__ == "__main__":
    dir = './'
    geography = 'Burkina Faso'
    prefix = "vector_genetics_spatial"

    exp = 8.6
    migration_mul = 100
    exp_name = "spatial_vector_genetics_serialization"

    # should be less than or equal to number of nodes requested by summary report and should ensure each core gets a job
    num_cores = 5

    num_years = 11
    num_seeds = 1

    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM',
                                        Num_Cores=num_cores,
                                        Simulation_Duration=int(365 * num_years))

    # set demographics file name and modify the config for the geography of interest
    cb.update_params({'Demographics_Filenames': [os.path.join(dir, 'Demographics', "%s_demographics.json" % prefix)]})
    cb.update_params({'Geography': geography})
    cb.update_params({'Vector_Species_Names': ['gambiae']})
    set_species(cb, ['gambiae'])
    set_species_param(cb, "gambiae", "Anthropophily", 0.65)

    # Spatial simulation + migration settings
    cb.update_params({
        # Match demographics file for constant population size (with exponential age distribution)
        "Age_Initialization_Distribution_Type": "DISTRIBUTION_COMPLEX",
        'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',
        'Enable_Nondisease_Mortality': 1,
        'New_Diagnostic_Sensitivity': 0.025,  # 40/uL
        'Disable_IP_Whitelist': 1,
        'Disable_NP_Whitelist': 1,

        "Vector_Migration_Filename_Local": os.path.join(dir, 'Migration', '%s_local_migration.bin' % prefix),
        'Enable_Vector_Migration': 1,  # mosquito migration
        'Enable_Vector_Migration_Local': 1,
        'Vector_Sampling_Type': "TRACK_ALL_VECTORS",
        'Vector_Sugar_Feeding_Frequency': 'VECTOR_SUGAR_FEEDING_NONE',
        "Roundtrip_Waypoints": 5,
        'Local_Migration_Filename': os.path.join(dir, 'Migration', '%s_local_migration.bin' % prefix),
        'Enable_Local_Migration': 1,
        'Migration_Model': 'FIXED_RATE_MIGRATION',  # human migration
        'Enable_Spatial_Output': 1,  # spatial reporting
        'Spatial_Output_Channels': ["Adult_Vectors", 'Infectious_Vectors', 'Population', 'Prevalence',
                                    'Daily_Bites_Per_Human'],
        "x_Local_Migration": 1,
        'x_Vector_Migration_Local': migration_mul,
        'logLevel_JsonConfigurable': 'ERROR',
        'logLevel_VectorHabitat': 'ERROR',
        'logLevel_Simulation': 'ERROR',
        'logLevel_StandardEventCoordinator': 'ERROR',
        'logLevel_LarvalHabitatMultiplier': 'ERROR',
        'logLevel_SimulationEventContext': 'ERROR',
        'logLevel_SusceptibilityMalaria': 'ERROR'
    })

    # Climate settings
    cb.update_params(({
        "Air_Temperature_Filename": os.path.join(dir, 'Climate', "%s_30arcsec_air_temperature_daily.bin" % geography),
        'Climate_Model': 'CLIMATE_CONSTANT',
        "Climate_Update_Resolution": "CLIMATE_UPDATE_DAY",
        "Land_Temperature_Filename": os.path.join(dir, 'Climate', "%s_30arcsec_air_temperature_daily.bin" % geography),
        "Rainfall_Filename": os.path.join(dir, 'Climate', "%s_30arcsec_rainfall_daily.bin" % geography),
        "Relative_Humidity_Filename": os.path.join(dir, 'Climate',
                                                   "%s_30arcsec_relative_humidity_daily.bin" % geography),
    }))

    cb.update_params({"Default_Geography_Initial_Node_Population": 1000,
                      "Default_Geography_Torus_Size": 10,
                      "Enable_Vector_Migration_Human": 0,
                      "Enable_Vector_Migration_Wind": 0,
                      "Egg_Hatch_Density_Dependence": "NO_DENSITY_DEPENDENCE",
                      "Temperature_Dependent_Feeding_Cycle": "NO_TEMPERATURE_DEPENDENCE",
                      "Enable_Drought_Egg_Hatch_Delay": 0,
                      "Enable_Egg_Mortality": 0,
                      "Enable_Temperature_Dependent_Egg_Hatching": 0,
                      "Vector_Migration_Modifier_Equation": "LINEAR"
                      })

    cb.set_param("Enable_Demographics_Builtin", 0)
    cb.set_param("Valid_Intervention_States", [])

    # Vector parameters
    set_species_param(cb, 'gambiae', 'Larval_Habitat_Types',
                      {"LINEAR_SPLINE": {
                          "Capacity_Distribution_Over_Time": {
                              "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
                                        182.5, 212.917, 243.333, 273.75, 304.167, 334.583],
                              "Values": [3, 0.8, 1.25, 0.1, 2.7, 10, 6, 35, 2.8, 1.5, 1.6, 2.1]
                          },
                          "Capacity_Distribution_Number_Of_Years": 1,
                          "Max_Larval_Capacity": pow(10, exp)
                      }})
    set_species_param(cb, "gambiae", "Adult_Life_Expectancy", 20)
    set_species_param(cb, "gambiae", "Male_Life_Expectancy", 10)
    set_species_param(cb, "gambiae", "Indoor_Feeding_Fraction", 0.9)

    add_SerializationTimesteps(cb, [num_years * 365], end_at_final=True)
    cb.update_params({
        "Report_Event_Recorder": 1,
        "Listed_Events": ["Immigrating", "Emigrating"],
        "Report_Event_Recorder_Events": ["Immigrating", "Emigrating"],
        "Report_Event_Recorder_Ignore_Events_In_List": 0
    })

    ########################## VECTOR GENETICS ####################################################
    # Add genes
    genes = {'gambiae': [{
        "Alleles": {"a0": 1.0, "a1": 0.0},
        "Mutations": {}
    }]
    }
    set_species_genes(cb, genes)

    builder = ModBuilder.from_list(
        [
            [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed)
             ]
            for seed in range(num_seeds)
        ]
    )

    ################################ INTERVENTIONS ############################################

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
                       start_day=(num_years-10)*365,
                       broadcast_event_name='Received_Treatment')
    #################################################################################################

    run_sim_args = {'config_builder': cb,
                    'exp_name': exp_name,
                    'exp_builder': builder}

    SetupParser.default_block = 'HPC'

    SetupParser.init('HPC')
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
