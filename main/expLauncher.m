%% General experimenter launcher
%  =============================
% By :      Martin SZINTE
% Projet :  pRFseqTest experiment
% With :    Anna MONTAGNINI & Guillaume MASSON
% Version:  1.0

% Version description
% ===================
% Experiment in which we first use a squqre full screen 4 direction (left/right/up/down) 
% bar pass stimuli with an attention task to the bar in order to pretest two whole brain multi-band
% acquisition sequences optimal for (occipital, parietal) frontal and subcortical structures

% +---------------------------------------+
% |          ACQUISITION SEQUENCES        |
% +-------------------+-------------------+
% |     acq2p5mm      |       acq2mm      |
% +-------------------+-------------------+
% | 2.5 mm isotropic  | 2.0 mm isotropic  |
% |   TR 1.2 seconds  |   TR 1.2 seconds  |
% |    Multi-band 3   |    Multi-band 4   |
% |     48 slices     |     60 slices     |
% +-------------------+-------------------+

% First settings
% --------------
Screen('CloseAll');clear all;clear mex;clear functions;close all;home;AssertOpenGL;

% General settings
% ----------------
const.expName           =   'pRFseqTest';   % experiment name.
const.expStart          =   1;              % Start of a recording exp                          0 = NO  , 1 = YES
const.checkTrial        =   0;              % Print trial conditions (for debugging)            0 = NO  , 1 = YES
const.writeLogTxt       =   1;              % write a log file in addition to eyelink file      0 = NO  , 1 = YES
const.genStimuli        =   0;              % Generate all stimuli                              0 = NO  , 1 = YES
const.drawStimuli       =   0;              % Draw stimuli generated                            0 = NO  , 1 = YES
const.mkVideo           =   0;              % Make a video of a run (on mac not linux)          0 = NO  , 1 = YES

% External controls
% -----------------
const.tracker           =   1;              % run with eye tracker                              0 = NO  , 1 = YES
const.scanner           =   1;              % run in MRI scanner                                0 = NO  , 1 = YES
const.scannerTest       =   0;              % run with T returned at TR time                    0 = NO  , 1 = YES
const.room              =   1;              % run in MRI or eye-tracking room                   1 = MRI , 2 = eye-tracking

% Run order
% ---------
const.cond_run_order = [1,1;1,2;...     %    run 01 - AttendStim_acq-2p5mm_run1  | run 02 - AttendStim_acq-2mm_run1
                        1,1;1,2;...     %    run 03 - AttendStim_acq-2p5mm_run2  | run 04 - AttendStim_acq-2mm_run2
                        1,1;1,2;...     %    run 05 - AttendStim_acq-2p5mm_run3  | run 06 - AttendStim_acq-2mm_run3
                        1,1;1,2;...     %    run 07 - AttendStim_acq-2p5mm_run4  | run 08 - AttendStim_acq-2mm_run4
                        1,1;1,2];       %    run 09 - AttendStim_acq-2p5mm_run5  | run 10 - AttendStim_acq-2mm_run5

% Run number per condition
% ------------------------
const.cond_run_num =  [1;1;...
                       2;2;...
                       3;3;...
                       4;4;...
                       5;5];

% Desired screen setting
% ----------------------
const.desiredFD         =   120;            % Desired refresh rate
%fprintf(1,'\n\n\tDon''t forget to change before testing\n');
const.desiredRes        =   [1920,1080];    % Desired resolution

% Path
% ----
dir                     =   (which('expLauncher'));
cd(dir(1:end-18));

% Add Matlab path
% ---------------
addpath('config','main','conversion','eyeTracking','instructions','trials','stim','stats');

% Subject configuration
% ---------------------
[const]                 =   sbjConfig(const);
                        
% Main run
% --------
main(const);

% Analyse data
% ------------
if const.runNum > 2;        res = input(sprintf('\n\tPlot results? YES (1) NO (0) : '));
elseif const.runNum == 10;  res = 1;
else                        res = 1;
end
if res; behav_results(const.sjct,const.runNum,const.tracker); end

