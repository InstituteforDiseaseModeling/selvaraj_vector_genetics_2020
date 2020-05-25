from dtk.interventions.mosquito_release import add_mosquito_release
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.reports import BaseVectorGeneticsReport, BaseVectorStatsReport
from dtk.vector.species import set_species_genes, set_species_param, set_species_drivers, set_species_trait_modifiers
from dtk.generic.climate import set_climate_constant
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModFn, ModBuilder
from simtools.SetupParser import SetupParser


def add_release(cb, number=1000, label='Mendelian'):

    species = 'funestus'
    add_mosquito_release(cb, start_day=200, species=species, repetitions=1, number=number,
                         released_genome=[['X', 'Y'], ['b0', 'b0']])
    add_mosquito_release(cb, start_day=200, species=species, repetitions=1, number=number,
                         released_genome=[['X', 'X'], ['b0', 'b0']])

    return {'Number': number, 'Label': label}


def add_genes(cb):

    # Add genes
    genes = {'arabiensis': [{
        "Alleles": {"a0": 1.0},
        "Mutations": {}
    }],
        'funestus': [{
            "Alleles": {"b1": 1.0, "b0": 0.0},
            "Mutations": {"b1:b0": 0.0, "b0:b1": 0.0}
        }],
    }
    set_species_genes(cb, genes)

    # Add gene trait modifiers
    traits = {'funestus': [
        {
            "Allele_Combinations": [["b1", "b1"]],
            "Trait_Modifiers": {"MORTALITY": 100}
        },
        {
            "Allele_Combinations": [["b0", "b0"]],
            "Trait_Modifiers": {"MORTALITY": 1}
        }
    ]
    }
    set_species_trait_modifiers(cb, traits)

    return None


if __name__ == "__main__":

    expname = 'vector_genetics_single_node_no_introgression'

    sim_duration = 365 * 10
    num_seeds = 50

    SetupParser('HPC')

    cb = DTKConfigBuilder.from_defaults('VECTOR_SIM')

    set_climate_constant(cb)
    set_species_param(cb, 'arabiensis', 'Larval_Habitat_Types',
                      {"CONSTANT": 11250000000})
    set_species_param(cb, 'arabiensis', 'Male_Life_Expectancy', 5)
    set_species_param(cb, 'arabiensis', 'Adult_Life_Expectancy', 10)
    set_species_param(cb, 'arabiensis', 'Transmission_Rate', 0.5)
    set_species_param(cb, 'arabiensis', 'Indoor_Feeding_Fraction', 1.0)

    set_species_param(cb, 'funestus', 'Larval_Habitat_Types',
                      {"CONSTANT": 11250000000})
    set_species_param(cb, 'funestus', 'Male_Life_Expectancy', 5)
    set_species_param(cb, 'funestus', 'Adult_Life_Expectancy', 10)
    set_species_param(cb, 'funestus', 'Transmission_Rate', 0.5)
    set_species_param(cb, 'funestus', 'Indoor_Feeding_Fraction', 1.0)

    ########################## VECTOR GENETICS ####################################################
    add_genes(cb)

    mendelian = [
        [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
         ModFn(add_release, label=label, number=10000),
         ]
        for seed in range(num_seeds)
        for label in ['No_inrogression']
    ]

    builder = ModBuilder.from_list(mendelian)

    # Call output reporters
    cb.add_reports(BaseVectorGeneticsReport(type='ReportVectorGenetics',
                                            species='arabiensis',
                                            gender='VECTOR_FEMALE',
                                            include_vector_state_columns=1,
                                            stratify_by='SPECIFIC_GENOME',
                                            combine_similar_genomes=1,
                                            specific_genome_combinations_for_stratification=[
                                                {
                                                    "Allele_Combination": [
                                                        ["X", "X"],
                                                        ["a0", "*"]
                                                    ]
                                                }
                                            ]
                                            ))
    cb.add_reports(BaseVectorGeneticsReport(type='ReportVectorGenetics',
                                            species='funestus',
                                            gender='VECTOR_FEMALE',
                                            include_vector_state_columns=1,
                                            stratify_by='SPECIFIC_GENOME',
                                            combine_similar_genomes=1,
                                            specific_genome_combinations_for_stratification=[
                                                {
                                                    "Allele_Combination": [
                                                        ["X", "X"],
                                                        ["b0", "b0"]
                                                    ]
                                                },
                                                {
                                                    "Allele_Combination": [
                                                        ["X", "X"],
                                                        ["b1", "*"]
                                                    ]
                                                }
                                            ]
                                            ))
    cb.add_reports(BaseVectorStatsReport(type="ReportVectorStats",
                                         stratify_by_species=1))

    cb.update_params({"Demographics_Filenames": ['VectorGeneticsSingle\Demographics.json'],
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
                      "Vector_Species_Names": ['arabiensis', 'funestus'],
                      "Enable_Vector_Aging": 1
                      })

    run_sim_args = {'config_builder': cb,
                    'exp_name': expname,
                    'exp_builder': builder}

    exp_manager = ExperimentManagerFactory.from_setup()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)