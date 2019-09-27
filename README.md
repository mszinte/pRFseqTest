# pRFseqTEST
# ==========
# By :      Martin SZINTE
# Projet :  pRFseqTest
# With :    Anna MONTAGNINI & GUILLAUME MASSON
# Version:  1.0

# Version description
# ===================
# Experiment in which we first use a squqre full screen 4 direction (left/right/up/down) 
# bar pass stimuli with an attention task to the bar in order to pretest two whole brain multi-band
# acquisition sequences optimal for (occipital, parietal) frontal and subcortical structures

# +---------------------------------------+
# |          ACQUISITION SEQUENCES        |
# +-------------------+-------------------+
# |    acq-2p5mm      |      acq-2mm      |
# +-------------------+-------------------+
# | 2.5 mm isotropic  | 2.0 mm isotropic  |
# |   TR 1.2 seconds  |   TR 1.2 seconds  |
# |    Multi-band 3   |    Multi-band 4   |
# |     48 slices     |     60 slices     |
# +-------------------+-------------------+


# MRI analysis
# ------------

# 1. run mriqc on mesocentre using mri_analysis pre_fit/mriqc_srun.py or mriqc_sbatch
