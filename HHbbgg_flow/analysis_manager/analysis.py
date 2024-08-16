import os, sys
from HHbbgg_flow.ttH_killer import process_ttH_vars, process_ttH_sideband
import logging
logger = logging.getLogger(__name__)

def run_analysis(args):
    if args.ttH_vars_bool:
        process_ttH_vars(args.ttH_vars_config, args.out_pq_size)
    elif args.ttH_sideband_bool:
        process_ttH_sideband(args.ttH_sideband_config, args.out_pq_size)
