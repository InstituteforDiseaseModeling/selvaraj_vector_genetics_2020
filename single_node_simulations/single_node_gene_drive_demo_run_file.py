from dtk.interventions.mosquito_release import add_mosquito_release
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.reports import BaseVectorGeneticsReport
from dtk.vector.species import set_species_genes, set_species_param, set_species_drivers, set_species
from dtk.generic.climate import set_climate_constant
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModFn, ModBuilder
from simtools.SetupParser import SetupParser

import numpy as np


def add_release(cb, label='Mendelian', number=1000):
    add_mosquito_release(cb, start_day=200, species='arabiensis', repetitions=1, number=number,
                         released_genome=[['X', 'X'], ['a1', 'a1']])

    return {'Label': label, 'Number': number}


def add_driver(cb, copy_to_likelihood=1.0):

    drivers = {'arabiensis': [{
        "Alleles_Driven": [
            {"Allele_To_Copy": "a1",
             "Allele_To_Replace": "a0",
             "Copy_To_Likelihood": {"a0": 1-copy_to_likelihood, "a1": copy_to_likelihood}
             }],
        "Driver_Type": "CLASSIC",
        "Driving_Allele": "a1"
    }]
    }

    set_species_drivers(cb, drivers)

    return {'Copy_To_Likelihood': copy_to_likelihood}


if __name__ == "__main__":

    expname = 'vector_genetics_single_node_gene_drive'

    sim_duration = 365 * 6
    num_seeds = 50

    cb = DTKConfigBuilder.from_defaults('VECTOR_SIM')

    set_climate_constant(cb)
    set_species(cb, ['arabiensis'])
    set_species_param(cb, 'arabiensis', 'Larval_Habitat_Types',
                      {"CONSTANT": 11250000000})
    set_species_param(cb, 'arabiensis', 'Male_Life_Expectancy', 5)
    set_species_param(cb, 'arabiensis', 'Adult_Life_Expectancy', 10)
    set_species_param(cb, 'arabiensis', 'Transmission_Rate', 0.5)
    set_species_param(cb, 'arabiensis', 'Indoor_Feeding_Fraction', 1.0)

    ########################## VECTOR GENETICS ####################################################
    # Add genes
    genes = {'arabiensis': [{
        "Alleles": {"a0": 1.0, "a1": 0.0},
        "Mutations": {"a0:a1": 0.0, "a1:a0": 0.0}
    }]
    }
    set_species_genes(cb, genes)

    copy_to_likelihoods = np.arange(0.5, 0.59, 0.1).tolist()

    # Add gene trait modifiers

    gene_drive = [
        [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
         ModFn(add_release, label='Gene_drive'),
         ModFn(add_driver, copy_to_likelihood=likelihood),
         ]
        for seed in range(num_seeds)
        for likelihood in copy_to_likelihoods
    ]

    mendelian = [
        [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed),
         ModFn(add_release, label='Mendelian', number=1000000)
         ]
        for seed in range(num_seeds)
    ]

    # builder = ModBuilder.from_list(gene_drive + mendelian)
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

    SetupParser('HPC')
    run_sim_args = {'config_builder': cb,
                    'exp_name': expname,
                    'exp_builder': builder}

    exp_manager = ExperimentManagerFactory.from_setup()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
