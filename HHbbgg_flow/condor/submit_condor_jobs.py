import json
import os, sys
import math
import subprocess
import shutil, fileinput

# Define the maximum number of files to process per job (as defined in `run_analysis.py`)
FILES_PER_JOB = 5  # Adjust this value as needed


# Get user proxy
status, proxy_info = subprocess.getstatusoutput("voms-proxy-info")
proxy_info = proxy_info.split('\n')
for line in proxy_info:
    if "path" in line:
        user_proxy = line.split(":")[-1].strip()
        
        
xrd_analysis_tarfile= "/uscms/home/idutta/nobackup/HiggsDNA_central_Run3/HHbbggFlow/hhbbgg-flow.tar.gz"
python_file="/uscms/home/idutta/nobackup/HiggsDNA_central_Run3/HHbbggFlow/scripts/run_analysis.py"
executable="exe_template.sh"
req_memory= "4096"
req_disk="10000"
req_ncpus="1"

# Path to the original SMsamples.json file
samples_file = '../metadata/samples/SMsamples.json'
samples_to_run = {
    "Run3_2022postEE" : {
        "ttHToGG" : ["nominal"]
    },
    "Run3_2022preEE" : {
        "ttHToGG" : ["nominal"]
    }
}
# Load the JSON data
with open(samples_file, 'r') as f:
    samples_data = json.load(f)
    
# Create output directories for job scripts and JSON files
output_dir = 'job_samples/'
os.makedirs(output_dir, exist_ok=True)
exe_dir = 'job_executables/'
os.makedirs(exe_dir, exist_ok=True)
log_dir='job_logs/'
os.makedirs(log_dir, exist_ok=True) 

# Prepare job list
job_list = []

# Split the samples into smaller JSON files
for era, sample_names in samples_to_run.items():
    for sname, s_types in sample_names.items():
        for stype in s_types:
            files = samples_data['data_eras'][era]['samples'][sname]['files'][stype]
            num_jobs = math.ceil(len(files) / FILES_PER_JOB)
    
            for job_index in range(num_jobs):
                # Calculate file range for this job
                start_index = job_index * FILES_PER_JOB
                end_index = start_index + FILES_PER_JOB
                job_files = files[start_index:end_index]

                # Create new JSON structure
                job_sample_info = {
                    sname: {
                        'xs': samples_data['data_eras'][era]['samples'][sname]['xs'],
                        'bf': samples_data['data_eras'][era]['samples'][sname]['bf'],
                        'files': {
                            stype: job_files
                        }
                        ,                    }
                }
        
                # Write to a new JSON file
                job_json_file = os.path.join(output_dir, f'{era}_{sname}_{stype}_job_{job_index}.json')
                with open(job_json_file, 'w') as jf:
                    json.dump(job_sample_info, jf, indent=4)
        
                # Add to job list
                job_list.append(job_json_file)
                output_file = os.path.join(log_dir, f'{era}_{sname}_{stype}_job_{job_index}.out')
                error_file = os.path.join(log_dir, f'{era}_{sname}_{stype}_job_{job_index}.err')
                log_file = os.path.join(log_dir, f'{era}_{sname}_{stype}_job_{job_index}.log')

                # Create a job_exe.sh script for this job
                job_jdl_file = os.path.join(exe_dir, f'{era}_{sname}_{stype}_job_{job_index}.jdl')
                shutil.copyfile("submit_template.txt",job_jdl_file)
                for line in fileinput.FileInput(job_jdl_file, inplace=1):
                    line=line.replace("GRID_PROXY", user_proxy)
                    line=line.replace("XRD_ANALYSIS_TARFILE",xrd_analysis_tarfile)
                    line=line.replace("CONFIG_FILE",job_json_file)
                    line=line.replace("PYTHON_FILE",python_file)
                    line=line.replace("EXECUTABLE",executable)
                    line=line.replace("OUTPUT",output_file)
                    line=line.replace("ERROR",error_file)
                    line=line.replace("LOG",log_file)
                    line=line.replace("REQ_MEMORY",req_memory)
                    line=line.replace("REQ_DISK",req_disk)
                    line=line.replace("REQ_NCPUS",req_ncpus)
                    print(line.rstrip())
                        
                submitCommand = "condor_submit "+job_jdl_file
                print(submitCommand)
                os.system(submitCommand)
                
# Write job list to file
with open('sample_jobs.txt', 'w') as jf:
    for job_file in job_list:
        jf.write(f"HHbbgg_flow/condor/{job_file}\n")
