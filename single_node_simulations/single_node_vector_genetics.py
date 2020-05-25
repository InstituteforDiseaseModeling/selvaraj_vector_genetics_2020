from dtk.interventions.mosquito_release import add_mosquito_release
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.builders.sweep import GenericSweepBuilder
from dtk.utils.reports import BaseVectorStatsReport, BaseVectorGeneticsReport
from dtk.vector.species import set_species_genes, set_species_param, set_species_trait_modifiers, set_species_drivers, \
    set_species
from dtk.generic.climate import set_climate_constant
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser


if __name__ == "__main__":

    expname = 'vector_genetics_single_node'

    sim_duration = 365 * 6
    num_seeds = 50

    cb = DTKConfigBuilder.from_defaults('VECTOR_SIM')

    builder = GenericSweepBuilder.from_dict({'Run_Number': range(num_seeds)})

    set_climate_constant(cb)
    set_species(cb, ['arabiensis'])
    set_species_param(cb, 'arabiensis', 'Larval_Habitat_Types',
                      {"TEMPORARY_RAINFALL": 11250000000})
    set_species_param(cb, 'arabiensis', 'Male_Life_Expectancy', 5)
    set_species_param(cb, 'arabiensis', 'Adult_Life_Expectancy', 10)
    set_species_param(cb, 'arabiensis', 'Transmission_Rate', 0.5)
    set_species_param(cb, 'arabiensis', 'Indoor_Feeding_Fraction', 1.0)

    ########################## VECTOR GENETICS ####################################################
    # Add genes
    genes = {'arabiensis': [{
        "Alleles": {"a0": 1.0, "a1": 0.0},
        "Mutations": {"a0:a1": 0.01, "a1:a0": 0.01}
    }]
    }
    set_species_genes(cb, genes)

    # Add gene trait modifiers
    traits = {'arabiensis': [
        {
            "Allele_Combinations": [["X", "X"],["a0", "a1"]],
            "Trait_Modifiers": {"INFECTED_BY_HUMAN": 0}
        },
        {
            "Allele_Combinations": [["X", "X"], ["a1", "a1"]],
            "Trait_Modifiers": {"INFECTED_BY_HUMAN": 0}
        },
        {
            "Allele_Combinations": [["X", "Y"], ["a0", "a1"]],
            "Trait_Modifiers": {"MORTALITY": 0.5}
        },
        {
            "Allele_Combinations": [["X", "Y"], ["a1", "a1"]],
            "Trait_Modifiers": {"MORTALITY": 0.5}
        },
    ]
    }
    set_species_trait_modifiers(cb, traits)

    # Call output reporters
    cb.add_reports(BaseVectorStatsReport(type='ReportVectorStats', stratify_by_species=1))
    cb.add_reports(BaseVectorGeneticsReport(type='ReportVectorGenetics',
                                            species='arabiensis',
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

###################################################################################################################


# Update other simulation parameters
    cb.update_params({"Demographics_Filenames": ['single/vector_genetics_single_node_demographics.json'],
                      "Enable_Natural_Mortality": 0,
                      "Birth_Rate_Dependence": "DEMOGRAPHIC_DEP_RATE",

                      "Base_Land_Temperature": 26.0,
                      "Base_Rainfall": 100,
                      "Base_Relative_Humidity": 0.5,

                      'x_Temporary_Larval_Habitat': 0.1,

                      "Simulation_Duration": sim_duration + 1,
                      "Run_Number": 0,
                      "Base_Air_Temperature": 25.0,
                      "Enable_Demographics_Birth": 0,
                      "Enable_Climate_Stochasticity": 1,

                      "Age_Dependent_Biting_Risk_Type": "OFF",

                      '.Vector_Sampling_Type': "TRACK_ALL_VECTORS",
                      'Vector_Sugar_Feeding_Frequency': 'VECTOR_SUGAR_FEEDING_NONE',
                      "Vector_Sampling_Type": "VECTOR_COMPARTMENTS_NUMBER",
                      "Vector_Species_Names": ['arabiensis'],
                      "Enable_Vector_Aging": 1
                      })

    SetupParser('LOCAL')
    run_sim_args = {'config_builder': cb,
                    'exp_name': expname,
                    'exp_builder': builder}

    exp_manager = ExperimentManagerFactory.from_setup()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
