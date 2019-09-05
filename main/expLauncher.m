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
% 3T room button press - setup Data Acquisition Toolbox
% 3T room TR triggers - setup Data Acquisition Toolbox and see settings
% Set sequence values of the TR - see Julien

% First settings
% --------------
Screen('CloseAll');clear all;clear mex;clear functions;close all;home;AssertOpenGL;

% General settings
% ----------------
const.expName           =   'pRFseqTest';   % experiment name.
const.expStart          =   0;              % Start of a recording exp                          0 = NO  , 1 = YES
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
const.room              =   1;              % run in MRI or eye-tracking room                   1 = MRI , 2 = eye-tracking

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