"""
-----------------------------------------------------------------------------------------
mriqc_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run fMRIprep on mesocentre using job mode
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main data directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: bids subject num (e.g. 01)
sys.argv[4]: server nb of processor to use (e.g 32)
sys.argv[5]: server nb of hour to request (e.g 10)
sys.argv[6]: anat only (1) or not (0)
sys.argv[7]: use of aroma (1) or not (0)
sys.argv[8]: Use low-quality tools (1) or not (0)
-----------------------------------------------------------------------------------------
Output(s):
preprocessed files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd /scratch/mszinte/projects/pRFseqTest/mri_analysis/pre_fit/
2. run python command
python fmriprep_sbatch.py [main directory] [project name] [subject num] [nb proc.] 
						  [hour proc.] [anat_only] [use_aroma] [use_sloppy]
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# imports modules
import sys
import os
import time
import ipdb
opj = os.path.join

# inputs
singularity_dir = '/scratch/mszinte/softwares/fmriprep-1.5.0.simg'
main_dir = sys.argv[1]
project_dir = sys.argv[2]
sub_num = sys.argv[3]
nb_procs = int(sys.argv[4])
hour_proc = int(sys.argv[5])
anat = int(sys.argv[6])
aroma = int(sys.argv[7])
sloppy = int(sys.argv[8])

# special input
anat_only, use_aroma, use_sloppy = '','',''
if anat == 1:
	anat_only = ' --anat-only'
if aroma == 1:
	use_aroma = ' --use-aroma'
if sloppy == 1:
	use_sloppy= ' --sloppy'

# define SLURM cmd
slurm_cmd = """\
#!/bin/sh
#SBATCH -J sub-{sub_num}_fmriprep
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -t {hour_proc}:00:00
#SBATCH --ntasks-per-node={nb_procs}
#SBATCH -o %N.%j.%a.out
#SBATCH -e %N.%j.%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=martin.szinte@gmail.com\n\n""".format(hour_proc = hour_proc, nb_procs = nb_procs, sub_num = sub_num)

# define singularity cmd
singularity_cmd = "singularity run --cleanenv -B {main_dir}:/work_dir {simg} /work_dir/{project_dir}/bids_data/ /work_dir/{project_dir}/deriv_data/fmriprep/ participant --participant_label {sub_num} -w /work_dir/{project_dir}/temp_data/ --nthreads {nb_procs:.0f} --fs-license-file /work_dir/freesurfer/license.txt --output-spaces T1w MNI152NLin2009cAsym fsaverage --verbose --notrack --no-submm-recon --cifti-output{anat_only}{use_aroma}{use_sloppy}".format(
									main_dir = main_dir,
									project_dir = project_dir,
									simg = singularity_dir,
									sub_num = sub_num,
									nb_procs = nb_procs,
									anat_only = anat_only,
									use_aroma = use_aroma,
									use_sloppy = use_sloppy
									)

# create sh folder and file
sh_dir = "{main_dir}/{project_dir}/deriv_data/fmriprep/jobs/sub-{sub_num}_fmriprep.sh".format(main_dir = main_dir, sub_num = sub_num,project_dir = project_dir,)

try:
	os.makedirs(opj(main_dir,project_dir,'deriv_data','fmriprep','jobs'))
	os.makedirs(opj(main_dir,project_dir,'deriv_data','fmriprep','log_outputs'))
except:
	pass

of = open(sh_dir, 'w')
of.write("{slurm_cmd}{singularity_cmd}".format(slurm_cmd = slurm_cmd,singularity_cmd = singularity_cmd))
of.close()

# Submit jobs
print("Submitting {sh_dir} to queue".format(sh_dir = sh_dir))
os.chdir(opj(main_dir,project_dir,'deriv_data','fmriprep','log_outputs'))
os.system("sbatch {sh_dir}".format(sh_dir = sh_dir))

