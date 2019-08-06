%% General experimenter launcher
%  =============================
% By :      Martin SZINTE
% Projet :  pRFseqTest experiment
% With :    Anna MONTAGNINI & Guillaume MASSON
% Version:  1.0

% Version description
% ===================
% Experiment in which we first use a full screen 4 direction (left/right/up/down) pRF
% bar pass stimuli with an attetion toward the bar to pretest 3T Prisma acquisition 
% sequences optimal for (occipital, parietal) frontal and subcortical structures.

% TO DO
% =====
% ProPixx settings and settings for Siemens Prisma
% Siemens button press values in new scanner
% Siemens triggers of first TR
% Sequence values of the TR
% Get informed consent form
% Get scanner consent form
% Collect behavioral data in testing room
% Make analysis behavior and eye tracking in python

% First settings
% --------------
Screen('CloseAll');clear all;clear mex;clear functions;close all;home;ListenChar(1);AssertOpenGL;

% General settings
% ----------------
const.expName           =   'pRFseqTest';   % experiment name.
const.sudopwd           =   'eyelink2';     % sudo pwd on machine
const.expStart          =   1;              % Start of a recording exp                          0 = NO  , 1 = YES
const.checkTrial        =   0;              % Print trial conditions (for debugging)            0 = NO  , 1 = YES
const.writeLogTxt       =   1;              % write a log file in addition to eyelink file      0 = NO  , 1 = YES
const.genStimuli        =   0;              % Generate all stimuli                              0 = NO  , 1 = YES
const.drawStimuli       =   0;              % Draw stimuli generated                            0 = NO  , 1 = YES
const.mkVideo           =   0;              % Make a video of a run (on mac not linux)          0 = NO  , 1 = YES

% External controls
% -----------------
const.tracker           =   1;              % run with eye tracker                              0 = NO  , 1 = YES
const.scanner           =   0;              % run in MRI scanner                                0 = NO  , 1 = YES
const.scannerTest       =   1;              % run with T returned at TR time                    0 = NO  , 1 = YES
const.room              =   2;              % run in MRI or eye-tracking room                   1 = MRI , 2 = eye-tracking

% Durations
% ---------
% 1 SESSION OF 35 min each of run + topup
% AttendStim: 5 x 3:00 min (150 TR)

% Run order
% ---------
const.cond_run_order = [1,1;1,2;...     %    run 01 - AttendStimSeq1_run1  | run 02 - AttendStimSeq2_run1
                        1,1;1,2;...     %    run 03 - AttendStimSeq1_run2  | run 04 - AttendStimSeq2_run2
                        1,1;1,2;...     %    run 05 - AttendStimSeq1_run3  | run 06 - AttendStimSeq2_run3
                        1,1;1,2;...     %    run 07 - AttendStimSeq1_run4  | run 08 - AttendStimSeq2_run4
                        1,1;1,2];       %    run 09 - AttendStimSeq1_run5  | run 10 - AttendStimSeq2_run5

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
addpath('config','main','conversion','eyeTracking','instructions','trials','stim');

% Subject configuration
% ---------------------
[const]                 =   sbjConfig(const);
                        
% Main run:
% ---------
main(const);