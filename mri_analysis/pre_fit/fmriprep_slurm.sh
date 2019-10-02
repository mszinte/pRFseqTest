#!/bin/bash
#!/bin/bash
#SBATCH --mail-type=ALL 			# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -p skylake
#SBATCH --mail-user=martin.szinte@univ-amu.fr	# Your email address
#SBATCH -A a163
#SBATCH --nodes=1					# OpenMP requires a single node
#SBATCH --mem=48gb
#####SBATCH --ntasks-per-node=8
####SBATCH --ntasks=1					# Run a single serial task
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00				# Time limit hh:mm:ss
#SBATCH -e ./%N.%j.%a.err			# Standard error
#SBATCH -o ./%N.%j.%a.out			# Standard output
#SBATCH -J fmriprep			# Descriptive job name
#SBATCH --mail-type=BEGIN,END
##### END OF JOB DEFINITION  #####

SUBJECT=$1

singularity run --cleanenv -B /scratch/jsein/BIDS:/work /scratch/jsein/my_images/fmriprep-1.5.0.simg --fs-license-file /work/freesurfer/license.txt /work/pRFseqTest /work/pRFseqTest/derivatives/fmriprep participant --participant-label 01 02 --use-aroma  --cifti-output
