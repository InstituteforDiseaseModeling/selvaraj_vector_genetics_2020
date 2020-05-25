Single node simulation use scripts
Prashanth Selvaraj, Apr 2020

-------------
Requirements:
dtk-tools package
dtk-tools-malaria package
malaria-toolbox package
These are available upon request from support@idmod.org

COMPS system for HPC job management
Email support@idmod.org

Input files and executable are elsewhere in the Additional File.

-------------
Scripts:

analyzers/genetic_data_analyzer.py - analyze genetic data output files

analyzers/inset_data_analyzer.py - analyze malaria epi data output files

plotting/plot_gene_drive.py - generates plots for single node gene drive simulations

plotting/plot_insecticide_resistance.py - generates plots for single node insecticide resistance simulations

plotting/plot_introgression.py - generates plots for single node introgression simulations

plotting/plot_vector_genetics.py - generates plots for single node vector genetics simulations

create_serialized_file.py - creates a population with immune profile matching a high transmission setting with 
no interventions through which different interventions scenarios can be explored

helper_functions.py - helper functions to set up single node simulations

insecticide_resistance_run_file.py - creates and runs insecticide resistance simulation scenarios

single_node_gene_drive_fitness_cost.py - creates and runs single node gene drive simulations with fitness costs

single_node_gene_drive_no_fitness_cost.py - creates and runs single node gene drive simulations without fitness costs

single_node_vector_genetics.py - creates and runs single node vector genetics simulation scenarios

single_node_no_introgression.py - creates and runs single node simulation without introgression

single_node_species_introgression.py - creates and runs single node simulation with introgression