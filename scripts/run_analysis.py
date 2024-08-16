import argparse
import os, sys
from pathlib import Path
sys.path.insert(1, os.getcwd())

from HHbbgg_flow.analysis_manager import analysis
from HHbbgg_flow.utils.logger_utils import setup_logger

def main(args):
    logger = setup_logger(args.log_level, args.log_file)
    logger.debug("Running HHbbggFlow")
    analysis.run_analysis(args)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Runs the analysis script.'
    )
    parser.add_argument('--ttH_vars', dest='ttH_vars_bool', required=False, action='store_true', default =False,
        help='Toggles whether or not to process ttH variables.'
    )
    parser.add_argument('--ttH_vars_config', dest='ttH_vars_config', required=False, action='store', default =None,
        help='JSON file that defines how the ttH_vars computations should be performed, including what files to run over.'
    )
    parser.add_argument('--ttH_sideband', dest='ttH_sideband_bool', required=False, action='store_true', default = False,
        help='Toggles whether or not to print ttH sideband cuts.'
    )
    parser.add_argument('--ttH_sideband_config', dest='ttH_sideband_config', required=False, action='store',default = None,
        help='JSON file that defines how the ttH_sideband computations should be performed, including what files to run over.'
    )
    parser.add_argument('--output_parquet_size', dest='out_pq_size', required=False, action='store', default = 0,
            help='*NOT IMPLEMENTED YET* Specifies the approx. size (in MB) of the output parquets. Useful with large datasets that require multi-processing. Defaults to 1 parquet per sample.'
        )  # Likely using Dask - https://stackoverflow.com/questions/63768642/pandas-df-to-parquet-write-to-multiple-smaller-files
    parser.add_argument('--dump', dest='output_dir_path', action='store', default=f'{str(Path().absolute())}/../../analysis_output',
        help='Name of the output path in which the analysis output will be stored.'
    )
    parser.add_argument('-d', '--debug', dest='DEBUG', action='store_true',
        help='Toggles whether or not to print debugging statements.'
    )

    parser.add_argument("--condor", required=False, action='store_true', default=None,
        help="option to run on condor"
    )
    parser.add_argument(
        "--log-level",
        required=False,
        default="INFO",
        type=str,
        help="Level of information printed by the logger")

    parser.add_argument(
        "--log-file",
        required=False,
        type=str,
        help="Name of the log file")

    args = parser.parse_args()
    main(args)


