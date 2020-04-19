import numpy as np
import math
import os
import json
import pandas as pd
import configparser

from dtk.utils.Campaign.CampaignClass import WaningEffectExponential
from malaria.study_sites.site_setup_functions import summary_report_fn
from dtk.interventions.mosquito_release import add_mosquito_release
from simtools.Utilities.Experiments import retrieve_experiment
from COMPS import Client
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.vector.species import set_species_trait_modifiers, set_species_drivers


def add_summary_report(cb, start_day=0, ipfilter='', description=''):
    summary_report_fn(start=start_day + 1, interval=30.0, description='Monthly_Report_%s' %description,
                      age_bins=[5.0, 15.0, 100.0],
                      parasitemia_bins=[0.0, 0.01, 0.1, 4.0, 40.0, 4000000.0],
                      infection_bins=np.arange(0.0, 101.0, 5.0).tolist(),
                      ipfilter=ipfilter,
                      nodes={'Node_List': [i for i in range(1, 89)], "class": "NodeSetNodeList"})(cb)

    return None


def add_release(cb, number=100, num_nodes=1, start_day=180):
    nodelist = find_n_largest(n=num_nodes)
    add_mosquito_release(cb, start_day=start_day, species='gambiae', repetitions=1, number=number,
                         released_genome=[['X', 'Y'], ['a1', 'a1']],
                         nodeIDs=nodelist)

    return {'Release_Number': number, 'Num_Nodes': num_nodes, 'Start_Day': start_day}


def add_release_single_node(cb, number=1000, start_day=180):
    add_mosquito_release(cb, start_day=start_day, species='gambiae', repetitions=1, number=number,
                         released_genome=[['X', 'Y'], ['a1', 'a1']])

    return {'Release_Number': number, 'Start_Day': start_day}


def add_nets(cb, coverage=0.6, number=0, num_nodes=0, start_day=180):

    # # ITN
    add_ITN_age_season(cb, 180, demographic_coverage=coverage,
                       discard_times={"Expiration_Distribution_Type": "EXPONENTIAL_DURATION",
                                      "Expiration_Period": 639},
                       blocking_config=WaningEffectExponential(Decay_Time_Constant=730, Initial_Effect=0.6),
                       killing_config=WaningEffectExponential(Decay_Time_Constant=1460, Initial_Effect=0.7)
                       )
    add_ITN_age_season(cb, 180 + 3 * 365, demographic_coverage=coverage,
                       discard_times={"Expiration_Distribution_Type": "EXPONENTIAL_DURATION",
                                      "Expiration_Period": 639},
                       blocking_config=WaningEffectExponential(Decay_Time_Constant=730, Initial_Effect=0.6),
                       killing_config=WaningEffectExponential(Decay_Time_Constant=1460, Initial_Effect=0.7)
                       )

    return {'ITN_Coverage': coverage, 'Release_Number': number, 'Num_Nodes': num_nodes, 'Start_Day': start_day}


def add_trait_modifiers(cb, transmission_to_human=0.0):

    # Add gene trait modifiers
    traits = {'gambiae': [
        {
            "Allele_Combinations": [["a0", "a1"]],
            "Trait_Modifiers": {"TRANSMISSION_TO_HUMAN": transmission_to_human}
        },
        {
            "Allele_Combinations": [["a1", "a1"]],
            "Trait_Modifiers": {"TRANSMISSION_TO_HUMAN": transmission_to_human}
        }
    ]
    }

    set_species_trait_modifiers(cb, traits)

    return {'Transmission_To_Human': transmission_to_human}


def add_drivers(cb, copy_to_likelihood=1.0):

    drivers = {'gambiae': [{
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


def find_n_largest(n=1):

    config = configparser.ConfigParser()
    config.read('simtools.ini')
    direc = config['HPC']['input_root']

    file = os.path.join(direc, 'Demographics', '%s_demographics.json' % 'vector_genetics')
    nodes = []
    pop = []
    with open(file) as f:
        data = json.load(f)
        for node in data['Nodes']:
            nodes.append(node['NodeID'])
            pop.append(node['NodeAttributes']['InitialPopulation'])
    df = pd.DataFrame({'nodes': nodes, 'pop': pop})
    df = df.sort_values(by=['pop'], ascending=False)

    node_list = list(df['nodes'][:n])

    return node_list
