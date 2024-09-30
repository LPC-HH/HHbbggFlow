#!/bin/bash

# Set proxy
export X509_USER_PROXY=$PWD/GRID_PROXY

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/usr/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/usr/etc/profile.d/conda.sh" ]; then
        . "/usr/etc/profile.d/conda.sh"
    else
        export PATH="/usr/bin:$PATH"
    fi
fi
unset __conda_setup

mkdir -p HHbbgg # dir to unload the conda env
xrdcp root://cmseos.fnal.gov//store/user/idutta/HiggsDNA_env_vars_Run3/HHbbgg_conda_env.tar.gz .
tar -xf HHbbgg_conda_env.tar.gz -C HHbbgg
ls -lrth HHbbgg/*/*
source HHbbgg/bin/activate

tar -xf hhbbgg-flow.tar.gz # dir to unload the HHbbggFlow repo contents
ls -lrth
python run_analysis.py

#for i in *.parquet; do xrdcp -f $i root://cmseos.fnal.gov//store/user/idutta/HiggsDNA/EOS_OUTPUT_DIR/SAMPLE/JOB/$i; rm -rf *.parquet; done
#for i in *summary*; do xrdcp -f $i root://cmseos.fnal.gov//store/user/idutta/HiggsDNA/EOS_OUTPUT_DIR/SAMPLE/JOB/$i; rm -rf *summary*; done
