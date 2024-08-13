# HHbbggFlow
The installation procedure consists in the following steps:
1. Clone this repository

git clone git@github.com:LPC-HH/HHbbggFlow.git
cd HHbbggFlow


2. Install miniconda
First get miniconda (skip if you already have it)

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow the instructions to finish the installation

# Choose your desired destination for installing minicondda (choose nobackup on CMSLPC)
# >Miniconda3 will now be installed into this location:
# ><your_home_directory_location>/miniconda3

 # > - Press ENTER to confirm the location
 # > - Press CTRL-C to abort the installation
 # > - Or specify a different location below

# >[<your_home_directory_location>/miniconda3] >>> /uscms/home/<your_username>/nobackup/miniconda3

# Make sure to choose `yes` for the following one to let the installer initialize Miniconda3
# > Do you wish the installer to initialize Miniconda3
# > by running conda init? [yes|no]


When logging in after a while (lxplus or elsewhere), you might have to do this to set your conda environment

source <path_to_your_miniconda_base>/miniconda3/bin/activate


Add this line to your ~/.bashrc file (using your favorite text editor)

export PATH="/uscms/home/<your_username>/nobackup/miniconda3/bin:$PATH"
##Note - For this to take effect, you might have to exit lpc and re-login

3. Install environment dependencies
If you already have a generic conda environment or a [HiggsDNA environment](https://gitlab.cern.ch/cms-analysis/general/HiggsDNA) installed:
```
conda activate <name-of-your-env>
conda install pip
pip install -r requirements.txt
```
And your environment will automatically be updated to have all needed packages.

For a BRAND NEW environment, the necessary dependencies can be installed by running:
```
conda env create -n HHbbgg -f full-package-list.txt
```

the conda env can become pretty large (multiple GB), so you may want to specify an installation location in the above step with
```
conda env create -f environment.yml -p <path to install conda env>
```

then activate the environment with
```
conda activate HHbbgg
```

If you specified an installation path, your activation command will instead be:
```
conda activate <path to install conda env>/HHbbgg
```
This repo provides the centralized location for the HHbbgg analysis framework. To run the analysis call:
```
python HHbbgg_flow/analysis_manager/run_analysis.py
```
or
```python
import HHbbgg_flow as hbg

hbg.analysis_runner(**kwargs)
```

This repo was developed and is maintained by the Fermilab-Purdue-Caltech collaboration.
