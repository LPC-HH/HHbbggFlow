import argparse
import glob
import json
import math
import os
from pathlib import Path
import re
import sys
import subprocess

import awkward as ak
import duckdb
import pyarrow.parquet as pq
import vector as vec
vec.register_awkward()

# lpc_redirector = "root://cmseos.fnal.gov/"
# lxplus_redirector = "root://eosuser.cern.ch/"
# lxplus_fileprefix = "/eos/cms/store/group/phys_b2g/HHbbgg/HiggsDNA_parquet/v1"
# LPC_FILEPREFIX = "/eos/uscms/store/group/lpcdihiggsboost/tsievert/HiggsDNA_parquet/v1"
FILEPREFIX = str()
OUTPUT_DIRPATH = f'{str(Path().absolute())}/../../ttH_killer_vars_output'
FORCE_RERUN = False
DEBUG = False

def add_ttH_vars(sample):
    
    # Funcs for Abs of cos and DeltaPhi #
    def ak_sign(ak_array, inverse=False):
        if not inverse:
            return ak.where(ak_array < 0, -1, 1)
        else:
            return ak.where(ak_array < 0, 1, -1)
            
    def ak_abs(ak_array):
        valid_entry_mask = ak.where(ak_array != -999, True, False)
        abs_ak_array = ak.where(ak_array > 0, ak_array, -ak_array)
        return ak.where(valid_entry_mask, abs_ak_array, -999)
        
    def deltaPhi(phi1, phi2):
        # angle1 and angle2 are (-pi, pi]
        # Convention: clockwise is (+), anti-clockwise is (-)
        subtract_angles = phi1 - phi2
        return ak.where(ak_abs(subtract_angles) <= math.pi, subtract_angles, subtract_angles + 2*math.pi*ak_sign(subtract_angles, inverse=True))


    # Funcs for chi^2 #
    def jets_mask(sample, jet_size, i, j, t_mask, i_mask=None, j_mask=None):
        jet_i_mask = t_mask & ak.where(
                sample[f'jet{i}_4mom'].deltaR(sample[f'lead_bjet_4mom']) > jet_size, True, False
            ) & ak.where(
                sample[f'jet{i}_4mom'].deltaR(sample[f'sublead_bjet_4mom']) > jet_size, True, False
            ) & ak.where(sample[f'jet{i}_pt'] != -999, True, False)
        if i_mask is not None:
            jet_i_mask = jet_i_mask & i_mask
        
        jet_j_mask = t_mask & ak.where(
                sample[f'jet{j}_4mom'].deltaR(sample[f'lead_bjet_4mom']) > jet_size, True, False
            ) & ak.where(
                sample[f'jet{j}_4mom'].deltaR(sample[f'sublead_bjet_4mom']) > jet_size, True, False
            ) & ak.where(sample[f'jet{j}_pt'] != -999, True, False)
        if j_mask is not None:
            jet_j_mask = jet_j_mask & j_mask
    
        return jet_i_mask, jet_j_mask

    def find_wjet_topjet(sample, num_jets, jet_size, w_mass, top_mass, t_mask, chi_form='t0'):
        jet_combos = []
        for i in range(1, num_jets+1):
            jet_combos.extend(
                [(i, j) for j in range(i+1, num_jets+1)]
            )
    
        chosen_w1jets = ak.Array(
            [
                {'i': -999, 'j': -999} for _ in range(ak.num(sample['event'], axis=0))
            ]
        )
        chosen_w1jets_deltaR = ak.Array(
            [0 for _ in range(ak.num(sample['event'], axis=0))]
        )
    
        
        for i, j in jet_combos:
            # Masks for non bjets and jet exists (jet_pt != -999)
            jet_i_mask, jet_j_mask = jets_mask(sample, jet_size, i, j, t_mask)
    
            # Select w-jets by minimizing deltaR between two not b-jets
            w1_decision_mask = (
                ak.where(
                    sample[f'jet{i}_4mom'].deltaR(sample[f'jet{j}_4mom']) < 
                    chosen_w1jets_deltaR, 
                    True, False
                ) | ak.where(
                    sample['w1jet_4mom'].mass == 0, True, False
                )
            ) & jet_i_mask & jet_j_mask
            
            sample['w1jet_4mom'] = ak.where(
                w1_decision_mask,
                sample[f'jet{i}_4mom'] + sample[f'jet{j}_4mom'], sample['w1jet_4mom']
            )
    
            
            chosen_w1jets['i'] = ak.where(
                w1_decision_mask,
                i, chosen_w1jets['i']
            )
            chosen_w1jets['j'] = ak.where(
                w1_decision_mask,
                j, chosen_w1jets['j']
            )
            chosen_w1jets_deltaR = ak.where(
                w1_decision_mask,
                sample[f'jet{i}_4mom'].deltaR(sample[f'jet{j}_4mom']), chosen_w1jets_deltaR
            )
    
        # Select bjet by minimizing deltaR between W-jet and b-jet
        bjet_mass_comparison_mask = ak.where(
            sample['w1jet_4mom'].deltaR(sample[f'lead_bjet_4mom']) <
            sample['w1jet_4mom'].deltaR(sample[f'sublead_bjet_4mom']),
            True, False
        )
        sample['top1jet_4mom'] = ak.where(
            bjet_mass_comparison_mask,
            sample['w1jet_4mom'] + sample[f'lead_bjet_4mom'], 
            sample['w1jet_4mom'] + sample[f'sublead_bjet_4mom']
        )
    
        if chi_form == 't1':
            # Select other wjet by dijet of reminaing two jets 
            for k, l in jet_combos:
                # Masks for non-bjets, not choosing same jets as w1, and jet exists (jet_pt != -999)
                jet_k_mask, jet_l_mask = jets_mask(
                    sample, jet_size, k, l, t_mask, 
                    i_mask=ak.where(chosen_w1jets['i'] != k, True, False),
                    j_mask=ak.where(chosen_w1jets['j'] != l, True, False)
                )
                
                w2_decision_mask = jet_k_mask & jet_l_mask
                sample['w2jet_4mom'] = ak.where(
                    w2_decision_mask,
                    sample[f'jet{k}_4mom'] + sample[f'jet{l}_4mom'], sample['w2jet_4mom']
                )
            # Select other bjet 
            sample['top2jet_4mom'] = ak.where(
                ~bjet_mass_comparison_mask,
                sample['w1jet_4mom'] + sample[f'lead_bjet_4mom'], 
                sample['w1jet_4mom'] + sample[f'sublead_bjet_4mom']
            )

    def chi_t0(sample, num_jets, jet_size):
        w_mass = 80.377
        top_mass = 172.76
    
        # To not include events with 4 extra jets, as its covered by chi_t1
        t_mask = ak.where(
            sample['jet4_pt'] != -999, True, False
        ) & ak.where(
            sample['jet6_pt'] == -999, True, False
        )
        
        find_wjet_topjet(
            sample, num_jets, jet_size, w_mass, top_mass, t_mask, chi_form='t0'
        )
        
        term1 = ((w_mass - ak.where(sample['w1jet_4mom'].mass == 0, -999, sample['w1jet_4mom'].mass)) / (0.1 * w_mass))**2
        term2 = ((top_mass - ak.where(sample['w1jet_4mom'].mass == 0, -999, sample['top1jet_4mom'].mass)) / (0.1 * top_mass))**2
    
        return ak.where(t_mask, term1+term2, -999)
        
    def chi_t1(sample, num_jets, jet_size):
        w_mass = 80.377
        top_mass = 172.76
    
        t_mask = ak.where(
            sample['jet6_pt'] != -999, True, False
        )
        
        find_wjet_topjet(
            sample, num_jets, jet_size, w_mass, top_mass, t_mask, chi_form='t1'
        )
    
        term1_1 = ((w_mass - ak.where(sample['w1jet_4mom'].mass == 0, -999, sample['w1jet_4mom'].mass)) / (0.1 * w_mass))**2
        term1_2 = ((top_mass - ak.where(sample['w1jet_4mom'].mass == 0, -999, sample['top1jet_4mom'].mass)) / (0.1 * top_mass))**2
    
        term2_1 = ((w_mass - ak.where(sample['w1jet_4mom'].mass == 0, -999, sample['w2jet_4mom'].mass)) / (0.1 * w_mass))**2
        term2_2 = ((top_mass - ak.where(sample['w1jet_4mom'].mass == 0, -999, sample['top2jet_4mom'].mass)) / (0.1 * top_mass))**2
    
        return ak.where(t_mask, term1_1+term1_2+term2_1+term2_2, -999)

    def deltaR_bjet_lepton(sample, lepton_type='lead', bjet_type='lead'):
        return sample[f'{lepton_type}_lepton_4mom'].deltaR(sample[f'{bjet_type}_bjet_4mom'])
    
    # Abs of cos #
    sample['abs_CosThetaStar_CS'] = ak.where(sample['CosThetaStar_CS'] >= 0, sample['CosThetaStar_CS'], -1*sample['CosThetaStar_CS'])
    sample['abs_CosThetaStar_jj'] = ak.where(sample['CosThetaStar_jj'] >= 0, sample['CosThetaStar_jj'], -1*sample['CosThetaStar_jj'])

    # DeltaPhi of (j, MET) #
    sample['DeltaPhi_j1MET'] = deltaPhi(sample['lead_bjet_phi'], sample['puppiMET_phi'])
    sample['DeltaPhi_j2MET'] = deltaPhi(sample['sublead_bjet_phi'], sample['puppiMET_phi'])

    # chi^2 #
    for field in ['lead', 'sublead']:
        sample[f'{field}_bjet_4mom'] = ak.zip(
            {
                'rho': sample[f'{field}_bjet_pt'], # rho is synonym for pt
                'phi': sample[f'{field}_bjet_phi'],
                'eta': sample[f'{field}_bjet_eta'],
                'tau': sample[f'{field}_bjet_mass'], # tau is synonym for mass
            }, with_name='Momentum4D'
        )

    for i in range(1, 7): # how to not hard-code 7 jets?
        sample[f'jet{i}_4mom'] = ak.zip(
            {
                'rho': sample[f'jet{i}_pt'],
                'phi': sample[f'jet{i}_phi'],
                'eta': sample[f'jet{i}_eta'],
                'tau': sample[f'jet{i}_mass'],
            }, with_name='Momentum4D'
        )
            
    for jet_type in ['w', 'top']:
        for jet_num in range(1, 3):
            # sample[f'{jet_type}{jet_num}jet_4mom'] = ak.copy(sample['zero_vector'])
            sample[f'{jet_type}{jet_num}jet_4mom'] = ak.Array(
                [
                    {'rho': 0, 'phi': 0, 'eta': 0, 'tau': 0} for _ in range(ak.num(sample['event'], axis=0))
                ], with_name='Momentum4D'
            )

    sample['chi_t0'] = chi_t0(sample, 6, 0.4)
    sample['chi_t1'] = chi_t1(sample, 6, 0.4)
    
    # lepton angulars #
    for field in ['lead', 'sublead']:
        sample[f'{field}_lepton_4mom'] = ak.zip(
            {
                'rho': sample[f'lepton{"1" if field == "lead" else "2"}_pt'], # rho is synonym for pt
                'phi': sample[f'lepton{"1" if field == "lead" else "2"}_phi'],
                'eta': sample[f'lepton{"1" if field == "lead" else "2"}_eta'],
                'tau': sample[f'lepton{"1" if field == "lead" else "2"}_mass'], # tau is synonym for mass
            }, with_name='Momentum4D'
        )
    
    sample['leadBjet_leadLepton'] = deltaR_bjet_lepton(sample)
    sample['leadBjet_subleadLepton'] = deltaR_bjet_lepton(sample, lepton_type='sublead')
    sample['subleadBjet_leadLepton'] = deltaR_bjet_lepton(sample, bjet_type='sublead')
    sample['subleadBjet_subleadLepton'] = deltaR_bjet_lepton(sample, lepton_type='sublead', bjet_type='sublead')

    # MET vars (only for v1 parquets) #
    if 'puppiMET_eta' not in set(sample.fields):
        sample['puppiMET_eta'] = [-999 for _ in range(ak.num(sample['event'], axis=0))]

def main(config: dict, out_pq_size: float):
    """
    Runs the script to compute the ttH-Killer variables
    """
    dir_lists = {data_era: list() for data_era in config['data_eras']}
    
    for data_era in dir_lists.keys():
        run_samples = []
        if os.path.exists(FILEPREFIX+'/'+data_era+'/processing_output.json'):
            with open(FILEPREFIX+'/'+data_era+'/processing_output.json', 'r') as f:
                processing_output = json.load(f)
            if processing_output['config_json'] == config and not FORCE_RERUN:
                run_samples = processing_output['run_samples']
        
        dir_lists[data_era] = list(
            (
                set(config["samples"]) & set(os.listdir(FILEPREFIX+'/'+data_era))
            ) - set(run_samples)
        )

    # MC Era: total era luminosity [fb^-1] #
    luminosities = {'Run3_2022preEE': 7.9804, 'Run3_2022postEE': 26.6717}
    
    # Name: cross section [fb] @ sqrrt{s}=13.6 TeV & m_H=125.09 GeV #
    #   -> Do we not need to care about other HH processes? https://arxiv.org/pdf/1910.00012.pdf
    cross_sections = {
        # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGHH?redirectedfrom=LHCPhysics.LHCHXSWGHH#Current_recommendations_for_HH_c
        'GluGluToHH': 34.43*0.0026,
        # https://xsdb-temp.app.cern.ch/xsdb/?columns=37748736&currentPage=0&pageSize=10&searchQuery=DAS%3DGG-Box-3Jets_MGG-80_13p6TeV_sherpa
        'GGJets': 88750, 
        # https://xsdb-temp.app.cern.ch/xsdb/?columns=37748736&currentPage=0&pageSize=10&searchQuery=DAS%3DGJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8
        'GJetPt20To40': 242500, 
        # https://xsdb-temp.app.cern.ch/xsdb/?columns=37748736&currentPage=0&pageSize=10&searchQuery=DAS%3DGJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8
        'GJetPt40': 919100, 
        # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#gluon_gluon_Fusion_Process
        'GluGluHToGG': 48520*0.00228,
        # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#ttH_Process
        'ttHToGG': 506.5*0.00228,
        # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#VBF_Process
        'VBFHToGG': 3779*0.00228,
        # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#WH_Process + https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#ZH_Process
        'VHToGG': (1369 + 882.4)*0.00228,
    }
    for data_era in dir_lists.keys():
        for dir_name in dir_lists[data_era]:
            if dir_name in cross_sections or re.match('Data', dir_name) is not None:
                continue
            cross_sections[dir_name] = 0.001 # set to 1e-3 [fb] for now, need to find actual numbers for many of these samples


    for data_era, dir_list in dir_lists.items():
        for dir_name in dir_list:
            for sample_type in ['nominal']: # Eventually change to os.listdir(LPC_FILEPREFIX+'/'+data_era+'/'+dir_name)
                # Load all the parquets of a single sample into an ak array
                sample = ak.concatenate(
                    [ak.from_parquet(FILEPREFIX+'/'+data_era+'/'+dir_name+'/'+sample_type+'/'+file) for file in os.listdir(FILEPREFIX+'/'+data_era+'/'+dir_name+'/'+sample_type+'/')]
                )
                add_ttH_vars(sample)
        
                if re.match('Data', dir_name) is None:
                    # Compute the sum of genWeights for proper MC rescaling.
                    sample['sumGenWeights'] = sum(
                        float(pq.read_table(file).schema.metadata[b'sum_genw_presel']) for file in glob.glob(
                            FILEPREFIX+'/'+data_era+'/'+dir_name+'/'+sample_type+'/*'
                        )
                    )
        
                    # Store luminostity computed from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis
                    #   and summing over lumis of the same type (e.g. all 22EE era lumis summed).
                    sample['luminosity'] = luminosities[data_era]
            
                    # If the process has a defined cross section, use defined xs otherwise use 1e-3 [fb] for now.
                    sample['cross_section'] = cross_sections[dir_name]
        
                    # Define eventWeight array for hist plotting.
                    sample['eventWeight'] = ak.where(sample['genWeight'] < 0, -1, 1) * (sample['luminosity'] * sample['cross_section'] / sample['sumGenWeights'])
        
                destdir = OUTPUT_DIRPATH + ('/' if OUTPUT_DIRPATH[-1] != '/' else '')
                processed_parquet = ak.to_parquet(sample, destdir+dir_name+'_'+sample_type+'.parquet')
                
                del sample
                if DEBUG:
                    print('======================== \n', dir_name)
                run_samples.append(dir_name)
                with open(FILEPREFIX+'/'+data_era+'/processing_output.json', 'w') as f:
                    processing_output = {
                        'run_samples': run_samples,
                        'config_json': config
                    }
                    json.dump(processing_output, f)

def check_config_json(config_json):
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
    parser.add_argument('--dump', dest='output_dir_path', action='store', default=f'{str(Path().absolute())}/../../ttH_killer_vars_output',
        help='Name of the output path in which the processed parquets will be stored.'
    )
    parser.add_argument('-f', '--force', dest='FORCE_RERUN', action='store_true',
        help='Forces the reprocessing to happen, even if the input config_json is identical to the previous processing.'
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
    # Global vars setting
    FILEPREFIX = config['file_prefix']
    FORCE_RERUN = args.FORCE_RERUN
    DEBUG = args.DEBUG
    # Logic for output_dir_path argument
    OUTPUT_DIRPATH = args.output_dir_path
    if not os.path.exists(OUTPUT_DIRPATH):
        os.makedirs(OUTPUT_DIRPATH)
    out_pq_size = args.out_pq_size

    main(config, out_pq_size)
