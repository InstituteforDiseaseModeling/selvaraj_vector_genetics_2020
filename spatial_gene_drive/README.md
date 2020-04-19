Spatial simulation use scripts
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

configure_interventions.py - created interventions for spatial simulations

create_serialized_file.py - creates a population with immune profile matching a high transmission setting with 
no interventions through which different interventions scenarios can be explored

helper_functions.py - helper functions to set up spatial simulations

main_run_file.py - creates and runs spatial simulation scenarios

plot_spatial_sim.py - generates plots for spatial simulation

analyzers/genetic_data_analyzer.py - analyze genetic data output files

analyzers/spatial_data_analyzer.py - analyze output spatial binary files

analyzers/spatial_output_dataframe.py - convert data from spatial binary files to a pandas dataframe