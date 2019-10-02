# pRFseqTEST
By :      Martin SZINTE<br/>
Projet :  pRFseqTest<br/>
With :    Anna MONTAGNINI & Guillaume MASSON<br/>
Version:  1.0<br/>

## Version description
Experiment in which we first use a square full screen 4 direction (left/right/up/down)<br/>
bar pass stimuli with an attention task to the bar in order to pretest two whole brain multi-band<br/>
acquisition sequences optimal for (occipital, parietal) frontal and subcortical structures<br/>

## Acquisition sequences
### acq-2p5mm<br/>

* 2.5 mm isotropic<br/>  
* TR 1.2 seconds<br/>
* Multi-band 3<br/>
* 48 slices<br/>

### acq-2mm<br/>

* 2.0 mm isotropic<br/>  
* TR 1.2 seconds<br/>
* Multi-band 4<br/>
* 60 slices<br/>

## MRI analysis
1. run mriqc on mesocentre using mri_analysis/pre_fit/mriqc_srun.py or mriqc_sbatch<br/>
2. run fmriprep on mesocenter using mri_analysis/pre_fit/fmriprep_sbatch - use first anat-only option
3. check/edit pial surface on freeview, see mri_analysis/pre_fit/pial_edit
