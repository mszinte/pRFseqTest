import subprocess as sb
import os
import glob
import ipdb
deb = ipdb.set_trace

subject = sys.argv[1]

# Define analysis parameters
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# base_dir = '/Users/martin/disks/meso_S/data/pRFseqTest/deriv_data/fmriprep/'
base_dir = analysis_info['base_dir_local']
fs_dir = "{base_dir}/deriv_data/fmriprep/{subject}".format(base_dir = base_dir, subject = subject)

t1mgz = 
mean_epi = 

anat_cmd = '-v {t1mgz}:grayscale=10,100'.format(t1mgz = t1mgz)
mean_epi_cmd = '-v {epi}:colormap=jet:heatscale=20,200,500:opacity=0.65'
volumes_cmd = '-f {lh_wm}:color=red:edgecolor=red -f {rh_wm}:color=red:edgecolor=red -f {lh_pial}:color=white:edgecolor=white -f {rh_pial}:color=white:edgecolor=white'.format()
slice_cmd = ' -slice {xpos} 127 127 \n -ss {opfn} \n'
freeview_cmd = 'freeview -cmd {anat_cmd} {mean_epi_cmd} {volumes_cmd} -viewport sagittal {slice_cmd}'.fomat(
                            anat_cmd = anat_cmd,
                            mean_epi_cmd = mean_epi_cmd,
                            volumes_cmd = volumes_cmd,
                            slice_cmd = slice_cmd)

 """



FS_folder = os.path.join(base_dir,'freesurfer',subject)
target_directory = '/Users/martin/disks/meso_S/data/pRFseqTest/deriv_data/fmriprep/images/{subject}'.format(subject = subject)
os.makedirs(target_directory, exist_ok=True)
cmd_file = os.path.join(target_directory, 'cmd.txt')

slices = range(0, 255)  #

sj_cmd = cmd_txt.format(
    anatomy=os.path.join(FS_folder, 'mri', 'T1.mgz'),
    lh_wm=os.path.join(FS_folder, 'surf', 'lh.white'),
    lh_pial=os.path.join(FS_folder, 'surf', 'lh.pial'),
    rh_wm=os.path.join(FS_folder, 'surf', 'rh.white'),
    rh_pial=os.path.join(FS_folder, 'surf', 'rh.pial'),
    subject=subject,
)

for sag_slice in slices:

    sj_cmd += slice_addition.format(
        xpos=sag_slice,
        opfn=os.path.join(target_directory, str(
            sag_slice).zfill(3) + '.png')
    )

sj_cmd += ' -quit \n '

with open(cmd_file, 'w') as f:
    f.write(sj_cmd)

# sb.call(freeview_command.format(cmd=cmd_file), shell=True)

convert_command = 'ffmpeg -framerate 5 -pattern_type glob -i "{target_directory}/*.png" -b:v 2M -c:v mpeg4 {target_directory}/{subject}.mp4'.format(target_directory = target_directory, subject = subject)
sb.call(convert_command, shell=True)
