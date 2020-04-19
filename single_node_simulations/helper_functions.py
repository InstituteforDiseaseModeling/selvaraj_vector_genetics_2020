# Functions to run insecticide resistance simulations
import os

from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.utils.Campaign.CampaignClass import WaningEffectExponential
from dtk.utils.reports import BaseVectorGeneticsReport, BaseVectorStatsReport
from dtk.vector.species import set_species


def change_vector_params(cb, species, mutation_rate1, mutation_rate2, serialization_seed,
                         initial_a0=0.99, initial_a1=0.01, initial_a2=0.0):

    set_species(cb, [species])
    Vector_Species_Params = [
        {
            "Name": species,
            "Genes": [
                {
                    "Alleles": {
                        "a0": initial_a0,
                        "a1": initial_a1,
                        "a2": initial_a2
                    },
                    "Mutations": {"a0:a1": mutation_rate1,
                                  "a1:a0": mutation_rate1,
                                  "a2:a0": mutation_rate2,
                                  "a0:a2": mutation_rate2,
                                  "a1:a2": mutation_rate2,
                                  "a2:a1": mutation_rate2
                                  }
                }
            ],
            "Acquire_Modifier": 0.8,
            "Adult_Life_Expectancy": 20,
            "Male_Life_Expectancy": 10,
            "Anthropophily": 0.65,
            "Aquatic_Arrhenius_1": 84200000000,
            "Aquatic_Arrhenius_2": 8328,
            "Aquatic_Mortality_Rate": 0.1,
            "Cycle_Arrhenius_1": 0,
            "Cycle_Arrhenius_2": 0,
            "Cycle_Arrhenius_Reduction_Factor": 0,
            "Days_Between_Feeds": 3,
            "Egg_Batch_Size": 100,
            "Immature_Duration": 2,
            "Indoor_Feeding_Fraction": 0.9,
            "Infected_Arrhenius_1": 117000000000,
            "Infected_Arrhenius_2": 8336,
            "Infected_Egg_Batch_Factor": 0.8,
            "Infectious_Human_Feed_Mortality_Factor": 1.5,
            "Larval_Habitat_Types":
                {"LINEAR_SPLINE": {
                          "Capacity_Distribution_Over_Time": {
                              "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
                                        182.5, 212.917, 243.333, 273.75, 304.167, 334.583],
                              "Values": [3, 0.8, 1.25, 0.1, 2.7, 10, 6, 35, 2.8, 1.5, 1.6, 2.1]
                          },
                          "Capacity_Distribution_Number_Of_Years": 1,
                          "Max_Larval_Capacity": pow(10, 8.6),
                    }
                },
            "Transmission_Rate": 0.9,
            "Vector_Sugar_Feeding_Frequency": "VECTOR_SUGAR_FEEDING_NONE"
        }
    ]
    cb.update_params({'Vector_Species_Params': Vector_Species_Params})

    return {'Mutation_Rate1': mutation_rate1, 'Mutation_Rate2': mutation_rate2,
            'Serialization_seed': serialization_seed, 'a0_initial': initial_a0,
            'a1_initial': initial_a1
            }


def add_insecticides(cb, pyrethroid_killing=1.0, carbamate_killing=1.0):
    Insecticides = [
        {
            "Name": "pyrethroid",
            "Resistances": [
                {
                    "Allele_Combinations": [
                        [
                            "a1",
                            "a1"
                        ]
                    ],
                    "Blocking_Modifier": 1.0,
                    "Killing_Modifier": pyrethroid_killing,
                    "Species": "gambiae"
                }
            ]
        },
        {
            "Name": "pyrethroid",
            "Resistances": [
                {
                    "Allele_Combinations": [
                        [
                            "a1",
                            "a0"
                        ]
                    ],
                    "Blocking_Modifier": 1.0,
                    "Killing_Modifier": pyrethroid_killing*2,
                    "Species": "gambiae"
                }
            ]
        },
        {
            "Name": "carbamate",
            "Resistances": [
                {
                    "Allele_Combinations": [
                        [
                            "a2",
                            "*"
                        ]
                    ],
                    "Blocking_Modifier": 1.0,
                    "Killing_Modifier": carbamate_killing,
                    "Species": "gambiae"
                }
            ]
        }
    ]

    cb.update_params({'Insecticides': Insecticides})

    return {'Pyrethroid_Killing': pyrethroid_killing, 'Carbamate_Killing': carbamate_killing}


def add_ITNs(cb, start=0, coverage=0.8, insecticide='', label='No_nets', killing=0.9):

    add_ITN_age_season(cb, start=start, demographic_coverage=coverage,
                       insecticide=insecticide,
                       discard_times={"Expiration_Distribution_Type": "EXPONENTIAL_DISTRIBUTION",
                                      "Expiration_Period_Exponential": 660},
                       blocking_config=WaningEffectExponential(Initial_Effect=0.9, Decay_Time_Constant=730),
                       killing_config=WaningEffectExponential(Initial_Effect=killing, Decay_Time_Constant=1460))

    cb.update_params({"Report_Event_Recorder": 1,
                      "Report_Event_Recorder_Events": ["Bednet_Got_New_One",
                                                       "Bednet_Discarded"],
                      "Report_Event_Recorder_Ignore_Events_In_List": 0,
                      "Report_Event_Recorder_Individual_Properties": []})

    return {'Insecticide': insecticide, 'Coverage': coverage, 'Start': start, 'Label': label}


def update_config_params(cb, direc, geography):
    # set demographics file name and modify the config for the geography of interest
    cb.update_params(
        {'Demographics_Filenames': [os.path.join(direc, 'input_files', "vector_genetics_single_demographics.json")]})
    cb.update_params({'Geography': geography})

    # Spatial simulation + migration settings
    cb.update_params({
        # Match demographics file for constant population size (with exponential age distribution)
        "Age_Initialization_Distribution_Type": "DISTRIBUTION_COMPLEX",
        'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',
        "Enable_Vital_Dynamics": 1,
        "Enable_Natural_Mortality": 1,
        'New_Diagnostic_Sensitivity': 0.025,  # 40/uL
        'Disable_IP_Whitelist': 1,
        'Disable_NP_Whitelist': 1,

        "Custom_Individual_Events": ['Bednet_Discarded', 'Bednet_Got_New_One', 'Bednet_Using',
                                     'Received_Treatment'],

        'Vector_Sampling_Type': "VECTOR_COMPARTMENTS_NUMBER",
        'Vector_Sugar_Feeding_Frequency': 'VECTOR_SUGAR_FEEDING_NONE',
    })

    # Climate settings
    cb.update_params(({
        'Climate_Model': 'CLIMATE_CONSTANT',
        "Climate_Update_Resolution": "CLIMATE_UPDATE_DAY",
    }))

    cb.update_params({
        "Default_Geography_Initial_Node_Population": 1000,
        "Default_Geography_Torus_Size": 10,
        "Egg_Hatch_Density_Dependence": "NO_DENSITY_DEPENDENCE",
        "Temperature_Dependent_Feeding_Cycle": "NO_TEMPERATURE_DEPENDENCE",
        "Enable_Drought_Egg_Hatch_Delay": 0,
        "Enable_Egg_Mortality": 0,
        "Enable_Temperature_Dependent_Egg_Hatching": 0,
        "Parasite_Smear_Sensitivity": 0.01,
        "Insecticides": []
    })

    cb.set_param("Enable_Demographics_Builtin", 0)


def add_reporters(cb, vector_genetics_report=0, vector_stats_report=0, malaria_summary_report=0):

    if vector_genetics_report:
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
                                                    },
                                                    {
                                                        "Allele_Combination": [
                                                            ["X", "X"],
                                                            ["a2", "*"]
                                                        ]
                                                    }
                                                ]
                                                ))

        if vector_stats_report:
            cb.add_reports(BaseVectorStatsReport(type='ReportVectorStats', stratify_by_species=1))
