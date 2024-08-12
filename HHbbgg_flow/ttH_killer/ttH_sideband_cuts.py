import argparse
import json
import os
from pathlib import Path
import re

import awkward as ak

# LPC_FILEPREFIX = "/eos/uscms/store/group/lpcdihiggsboost/tsievert/HiggsDNA_parquet/v1"
FILEPREFIX = str()
OUTPUT_DIRPATH = f'{str(Path().absolute())}/../../ttH_sideband_cuts_output'
SINGLE_B_WPS = {
    'preEE': {'L': 0.047, 'M': 0.245, 'T': 0.6734, 'XT': 0.7862, 'XXT': 0.961},
    'postEE': {'L': 0.0499, 'M': 0.2605, 'T': 0.6915, 'XT': 0.8033, 'XXT': 0.9664}
}
SIDEBAND_MASK = 'SIDEBAND_MASK'
FILL_VALUE = -999

def sideband_cuts(data_era: str, sample):
    """
    Builds the event_mask used to do data-mc comparison in a sideband.
    """
    # Require diphoton and dijet exist (should be required in preselection, and thus be all True)
    event_mask = ak.where(sample['pt'] != FILL_VALUE, True, False) & ak.where(sample['dijet_pt'] != FILL_VALUE, True, False)
    # Require btag score above Loose WP
    EE_era_2022 = 'preEE' if re.search('preEE', data_era) is not None else 'postEE'
    event_mask = event_mask & ak.where(
        sample['lead_bjet_btagPNetB'] > SINGLE_B_WPS[EE_era_2022]['L'], True, False
    ) & ak.where(
        sample['sublead_bjet_btagPNetB'] > SINGLE_B_WPS[EE_era_2022]['L'], True, False
    )
    # Require at least 3 jets (to remove bbH background), extra jets coming from Ws
    event_mask = event_mask & ak.where(sample['jet3_pt'] != FILL_VALUE, True, False)
    # Require events with diphoton mass within Higgs window
    event_mask = event_mask & (
        ak.where(sample['mass'] >= 100, True, False) & ak.where(sample['mass'] <= 150, True, False)
    )
    # Mask out events with dijet mass within Higgs window
    event_mask = event_mask & (
        ak.where(sample['dijet_mass'] <= 70, True, False) | ak.where(sample['dijet_mass'] >= 150, True, False)
    )
    sample = ak.mask(sample, event_mask)

def main(config: dict, out_pq_size: str):
    """
    Runs the script to compute the ttH-Killer variables
    """
    dir_lists = {data_era: list() for data_era in config['data_eras']}
    
    for data_era in dir_lists.keys():
        dir_lists[data_era] = list(
            (
                set(config["samples"]) & set(os.listdir(FILEPREFIX+'/'+data_era))
            )
        )

    for data_era, dir_list in dir_lists.items():
        for dir_name in dir_list:
            for sample_type in ['nominal']: # Eventually change to os.listdir(LPC_FILEPREFIX+'/'+data_era+'/'+dir_name)
                # Load all the parquets of a single sample into an ak array
                sample = ak.concatenate(
                    [ak.from_parquet(FILEPREFIX+'/'+data_era+'/'+dir_name+'/'+sample_type+'/'+file) for file in os.listdir(FILEPREFIX+'/'+data_era+'/'+dir_name+'/'+sample_type+'/')]
                )
                sideband_cuts(sample)
        
                destdir = OUTPUT_DIRPATH + ('/' if OUTPUT_DIRPATH[-1] != '/' else '')
                masked_parquet = ak.to_parquet(sample, destdir+dir_name+'_'+sample_type+'.parquet')
                
                del sample
                if DEBUG:
                    print('======================== \n', dir_name)

def check_config_json(config_json): 
    """
    Checks if the config_json is valid. Currently identical btwn ttH_vars and sideband_cuts... may not eventually be true.
      If it stays the same, could be moved to a utils.py file.
    """
    if not os.path.exists(config_json):
        print("You provided a JSON filename, but the path doesn't exist. Maybe you mistyped the path?")
        raise SystemExit(1)
    with open(config_json, 'r') as f:
        config = json.load(f)

    minimal_config_json = {
        "data_eras", "samples", "file_prefix"
        # Should we require they specify variables to compute? or we can assume all if not passed
    }
    if set(config.keys()) < minimal_config_json:
        print(f"You provided a valid JSON filepath, however its missing some of the required keys: \n{minimal_config_json - set(config.keys())}")
        raise SystemExit(1)
    elif set(config.keys()) > minimal_config_json:
        print(f"WARNING: You provided a valid JSON file, however it contains some extra (unused) keys: \n{set(config.keys()) - minimal_config_json}")

    return config

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Process the v1 HiggsDNA parquets to add the ttH-killer variables.'
    )
    parser.add_argument('config_json', dest='config_json', action='store',
        help='JSON file that defines how the computations should be performed, including what files to run over.'
    )
    parser.add_argument('--dump', dest='output_dir_path', action='store', default=f'{str(Path().absolute())}/../../ttH_sideband_cuts_output',
        help='Name of the output path in which the processed parquets will be stored.'
    )
    parser.add_argument('-d', '--debug', dest='DEBUG', action='store_true',
        help='Toggles whether or not to print debugging statements.'
    )
    parser.add_argument('--output_parquet_size', dest='out_pq_size', action='store',
        help='*NOT IMPLEMENTED YET* Specifies the approx. size (in MB) of the output parquets. Useful with large datasets that require multi-processing. Defaults to 1 parquet per sample.'
    )  # Likely using Dask - https://stackoverflow.com/questions/63768642/pandas-df-to-parquet-write-to-multiple-smaller-files

    args = parser.parse_args()
    # Logic for config_json argument
    config = check_config_json(args.config_json)
    FILEPREFIX = config['file_prefix']
    DEBUG = args.DEBUG
    # Logic for output_dir_path argument
    OUTPUT_DIRPATH = args.output_dir_path
    if not os.path.exists(OUTPUT_DIRPATH):
        os.makedirs(OUTPUT_DIRPATH)
    out_pq_size = args.out_pq_size

    main(config, out_pq_size)