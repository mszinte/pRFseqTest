"""
-----------------------------------------------------------------------------------------
load_dataset.py
-----------------------------------------------------------------------------------------
Goal of the script:
Load dataset and create webviewer
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject 
sys.argv[2]: acq
sys.argv[3]: model
sys.argv[4]: prf_sign
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
python load_dataset.py sub-01 acq-2mm gauss pos
python load_dataset.py sub-01 acq-2p5mm gauss pos
python load_dataset.py sub-02 acq-2mm gauss pos
python load_dataset.py sub-02 acq-2p5mm gauss pos
> to run localy
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import json
import numpy as np
import ipdb
import platform
import time
import matplotlib.pyplot as plt
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
from utils import set_pycortex_config_file

# Get inputs
# ----------
subject = sys.argv[1]
acq = sys.argv[2]
model = sys.argv[3]
prf_sign = sys.argv[4]

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define folder
# -------------
base_dir = analysis_info['base_dir_local']

# Set pycortex db and colormaps
# -----------------------------
set_pycortex_config_file(base_dir)

# Pycortex plots
dataset_file = "{base_dir}/pp_data/{subject}/{model}/pycortex_outputs/dataset/{prf_sign}/{acq}_{prf_sign}.hdf".format(
                            base_dir = base_dir, subject = subject, acq = acq, model = model, prf_sign = prf_sign)

ds = cortex.dataset.Dataset.from_file(dataset_file)
cortex.webgl.show(data=ds)
time.sleep(5)
