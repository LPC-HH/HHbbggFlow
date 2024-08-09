import argparse
import os
from pathlib import Path

from ttH_killer import process_ttH_vars, process_ttH_sideband

parser = argparse.ArgumentParser(
    description='Runs the analysis script.'
)
parser.add_argument('--ttH_vars', dest='ttH_vars_bool', action='store_true',
    help='Toggles whether or not to process ttH variables.'
)
parser.add_argument('--ttH_vars_config', dest='ttH_vars_config', action='store',
    help='JSON file that defines how the ttH_vars computations should be performed, including what files to run over.'
)
parser.add_argument('--ttH_sideband', dest='ttH_sideband_bool', action='store_true',
    help='Toggles whether or not to print ttH sideband cuts.'
)
parser.add_argument('--ttH_sideband_config', dest='ttH_sideband_config', action='store',
    help='JSON file that defines how the ttH_sideband computations should be performed, including what files to run over.'
)
parser.add_argument('--output_parquet_size', dest='out_pq_size', action='store',
        help='*NOT IMPLEMENTED YET* Specifies the approx. size (in MB) of the output parquets. Useful with large datasets that require multi-processing. Defaults to 1 parquet per sample.'
    )  # Likely using Dask - https://stackoverflow.com/questions/63768642/pandas-df-to-parquet-write-to-multiple-smaller-files
parser.add_argument('--dump', dest='output_dir_path', action='store', default=f'{str(Path().absolute())}/../../analysis_output',
    help='Name of the output path in which the analysis output will be stored.'
)
parser.add_argument('-d', '--debug', dest='DEBUG', action='store_true',
    help='Toggles whether or not to print debugging statements.'
)

def main():
    args = parser.parse_args()

    if args.ttH_vars_bool:
        process_ttH_vars(args.ttH_vars_config, args.out_pq_size)
    elif args.tth_sideband_bool:
        process_ttH_sideband(args.ttH_sideband_config, args.out_pq_size)

