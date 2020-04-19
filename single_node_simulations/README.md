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

create_serialized_file.py - creates a population with immune profile matching a high transmission setting with 
no interventions through which different interventions scenarios can be explored

helper_functions.py - helper functions to set up single node simulations

insecticide_resistance_run_file.py - creates and runs insecticide resistance simulation scenarios

plot_insecticide_resistance.py - generates plots for single node simulation

single_node_gene_drive_demo_run_file.py - creates and runs single node gene drive simulation scenarios

single_node_genetics_demo_run_file.py - creates and runs simple single node vector genetics simulation scenarios

analyzers/genetic_data_analyzer.py - analyze genetic data output files

analyzers/inset_data_analyzer.py - analyze malaria epi data output files